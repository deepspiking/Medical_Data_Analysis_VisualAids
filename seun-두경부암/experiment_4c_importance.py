# -*- coding: utf-8 -*-
"""
실험 4-7b: 주요 변수 추출·비교 (permutation importance, SHAP 대체)
------------------------------------------------------------------
SHAP이 TabICL segfault로 불가 → **permutation importance** 사용:
각 (y, 시드, fold)에서 TabICL fit 후, 변수별로 test X의 해당 열을
shuffle했을 때 score와 원본 score의 **순위 상관(Spearman ρ) 감소량**을 중요도로.

출력:
  exp1_tabicl_importance.csv  : (y, seed, fold, var, importance) 전체
  exp1_tabicl_importance_summary.csv : y별 시드/폴드 간 주요 변수 순위 일치도
"""
import os
import sys
import time
import json
import base64
import subprocess
import tempfile
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_CSV = os.path.join(BASE_DIR, "preprocessed_data.csv")
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")
os.makedirs(OUT_DIR, exist_ok=True)

SEEDS = [42, 123, 2026, 777]
N_FOLDS = 5
N_PERM = 3          # permutation 반복 (shuffle 횟수)
T_HORIZON = 24.0
Y_LABELS = ["PFS", "DSS", "LRRFS"]

COVS = ["age", "성별", "tumor size (cm)", "DOI (mm)", "differentiation",
        "budding_01vs23", "PNI", "LVI", "RM", "TIL", "TSR",
        "WPOI5_2tier", "HPV/P16_1", "HPV/P16_2", "CCRT_bin"]


def load_prep():
    d = pd.read_csv(DATA_CSV)
    flag = [c for c in d.columns if c.endswith("_known")]
    oh = [c for c in d.columns if c.startswith(("subsite_", "Tx_3tier_", "HPV/P16_"))]
    X = COVS + flag + oh
    X = list(dict.fromkeys(X))  # 중복 제거
    return d, X


def fold_split(n, y24, seed, n_folds=N_FOLDS):
    """사건층화 5-fold 분할 (실험 4 체크포인트와 동일 재현)."""
    rng = np.random.RandomState(seed)
    pos = np.nonzero(y24 == 1)[0]
    neg = np.nonzero(y24 == 0)[0]
    rng.shuffle(pos)
    rng.shuffle(neg)
    pos_folds = np.array_split(pos, n_folds)
    neg_folds = np.array_split(neg, n_folds)
    folds = []
    for f in range(n_folds):
        test_idx = np.concatenate([pos_folds[f], neg_folds[f]])
        train_idx = np.concatenate([
            np.concatenate([pos_folds[j] for j in range(n_folds) if j != f]),
            np.concatenate([neg_folds[j] for j in range(n_folds) if j != f])])
        folds.append((train_idx, test_idx))
    return folds


def _perm_fit_subprocess(args):
    """fit + permutation importance를 subprocess로 실행 (segfault 격리)."""
    X_train, y_train, X_test, seed, var_idx, out_prefix = args
    script = f'''
import os, sys, json, numpy as np
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, {BASE_DIR!r})
from tabicl import TabICLClassifier
from scipy.stats import spearmanr
X_train = np.load({out_prefix!r} + "_X.npy")
y_train = np.load({out_prefix!r} + "_y.npy")
X_test = np.load({out_prefix!r} + "_test.npy")
clf = TabICLClassifier(n_estimators=1, random_state={seed})
clf.fit(X_train, y_train)
score = clf.predict(X_test).astype(np.float64)
# permutation importance: 변수별 shuffle → 원본 score와의 Spearman ρ 감소량
vi = {var_idx!r}
n_test = X_test.shape[0]
orig_rank = np.argsort(np.argsort(score))   # 원본 score 순위
imp = {{}}
rng = np.random.RandomState({seed})
for v in vi:
    drop = []
    for _ in range({N_PERM}):
        Xp = X_test.copy()
        Xp[:, v] = rng.permutation(Xp[:, v])
        s_perm = clf.predict(Xp).astype(np.float64)
        rho = spearmanr(orig_rank, s_perm)[0]
        drop.append(rho)
    imp[str(v)] = float(np.mean(drop))
json.dump(imp, open({out_prefix!r} + "_imp.json", "w"))
os._exit(0)
'''
    np.save(out_prefix + "_X.npy", X_train)
    np.save(out_prefix + "_y.npy", y_train)
    np.save(out_prefix + "_test.npy", X_test)
    env = dict(os.environ)
    env["KMP_DUPLICATE_LIB_OK"] = "TRUE"
    env["OMP_NUM_THREADS"] = "1"
    try:
        p = subprocess.Popen([sys.executable, "-c", script],
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                             start_new_session=True, env=env)
        try:
            p.communicate(timeout=1800)
        except subprocess.TimeoutExpired:
            p.kill()
            p.communicate()
        imp_path = out_prefix + "_imp.json"
        if os.path.exists(imp_path):
            imp = json.load(open(imp_path))
            return imp
        return None
    except Exception:
        return None
    finally:
        for suffix in ("_X.npy", "_y.npy", "_test.npy", "_imp.json"):
            try:
                os.remove(out_prefix + suffix)
            except OSError:
                pass


def main():
    t0 = time.time()
    d, X = load_prep()
    n = len(d)
    print(f"주요 변수 추출 | X {len(X)}개, 시드 {len(SEEDS)} × fold {N_FOLDS} × y {len(Y_LABELS)}",
          flush=True)

    tmpdir = tempfile.mkdtemp(prefix="tabicl_perm_")
    ckpt_csv = os.path.join(OUT_DIR, "exp1_tabicl_importance.csv")
    done = set()
    if os.path.exists(ckpt_csv):
        prev = pd.read_csv(ckpt_csv)
        done = set(zip(prev["y"], prev["seed"], prev["fold"]))
        print(f"체크포인트 로드: {len(done)}개 fold 완료, 이어서 진행", flush=True)

    rows = []
    if os.path.exists(ckpt_csv):
        rows = pd.read_csv(ckpt_csv).to_dict("records")

    for y in Y_LABELS:
        y_time, y_event = f"{y}_time", f"{y}_event"
        dfx = d[X + [y_time, y_event]].dropna().reset_index(drop=True)
        dfx["y24"] = ((dfx[y_event] == 1) & (dfx[y_time] <= T_HORIZON)).astype(int)
        X_all = dfx[X].values.astype(float)
        y_all = dfx["y24"].values.astype(int)
        for seed in SEEDS:
            folds = fold_split(len(dfx), y_all, seed)
            for f, (train_idx, test_idx) in enumerate(folds):
                key = (y, seed, f)
                if key in done:
                    continue
                # 확장 bootstrap (체크포인트와 동일: seed*1000 + f*10 + 0)
                rng_try = np.random.RandomState(seed * 1000 + f * 10)
                boot_idx = rng_try.choice(train_idx, size=100, replace=True)
                out_prefix = os.path.join(tmpdir, f"perm_{y}_{seed}_{f}")
                imp = _perm_fit_subprocess(
                    (X_all[boot_idx], y_all[boot_idx], X_all[test_idx],
                     seed * 100 + f * 10, list(range(len(X))), out_prefix))
                if imp is None:
                    print(f"  [{y}] seed={seed} fold{f+1}/5 실패", flush=True)
                    continue
                for v_idx, v in enumerate(X):
                    rows.append({"y": y, "seed": seed, "fold": f,
                                 "var": v, "importance": imp.get(str(v_idx), np.nan)})
                # 즉시 저장 (segfault로 부모가 죽어도 진행 보존)
                pd.DataFrame(rows).to_csv(ckpt_csv, index=False, encoding="utf-8-sig")
                done.add(key)
                print(f"  [{y}] seed={seed} fold{f+1}/5 완료 ({time.time()-t0:.0f}s)",
                      flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT_DIR, "exp1_tabicl_importance.csv"),
               index=False, encoding="utf-8-sig")

    # 요약: y별 주요 변수 순위 + 시드/폴드 간 일치도
    if len(res):
        summ_rows = []
        for y in Y_LABELS:
            sub = res[res["y"] == y]
            if sub.empty:
                continue
            # 변수별 평균 중요도 (전체 시드/폴드)
            var_mean = sub.groupby("var")["importance"].mean().sort_values(ascending=False)
            top5 = list(var_mean.head(5).index)
            # 시드 간 순위 일치 (seed별 top 변수)
            seed_top = {}
            for seed in SEEDS:
                s = sub[sub["seed"] == seed].groupby("var")["importance"].mean()
                seed_top[seed] = list(s.sort_values(ascending=False).head(5).index)
            # 폴드 간 순위 일치 (fold별 top 변수)
            fold_top = {}
            for f in range(N_FOLDS):
                ff = sub[sub["fold"] == f].groupby("var")["importance"].mean()
                fold_top[f] = list(ff.sort_values(ascending=False).head(5).index)
            # seed 쌍 간 top5 일치율
            seed_agree = []
            seeds = list(seed_top.keys())
            for a in range(len(seeds)):
                for b in range(a + 1, len(seeds)):
                    inter = len(set(seed_top[seeds[a]]) & set(seed_top[seeds[b]]))
                    seed_agree.append(inter / 5)
            # fold 쌍 간 top5 일치율
            fold_agree = []
            folds = list(fold_top.keys())
            for a in range(len(folds)):
                for b in range(a + 1, len(folds)):
                    inter = len(set(fold_top[folds[a]]) & set(fold_top[folds[b]]))
                    fold_agree.append(inter / 5)
            summ_rows.append({
                "y": y,
                "top5_vars": ", ".join(top5),
                "seed_top5_agreement_mean": round(float(np.mean(seed_agree)), 3) if seed_agree else np.nan,
                "fold_top5_agreement_mean": round(float(np.mean(fold_agree)), 3) if fold_agree else np.nan,
                "n_rows": len(sub),
            })
            print(f"\n=== {y} 주요 변수 (전체 평균 중요도 상위 5) ===")
            print(var_mean.head(8).round(4).to_string())
            print(f"  seed 간 top5 일치율 평균: {np.mean(seed_agree):.2f}" if seed_agree else "")
            print(f"  fold 간 top5 일치율 평균: {np.mean(fold_agree):.2f}" if fold_agree else "")
        pd.DataFrame(summ_rows).to_csv(
            os.path.join(OUT_DIR, "exp1_tabicl_importance_summary.csv"),
            index=False, encoding="utf-8-sig")

    print(f"\n[완료] exp1_tabicl_importance.csv (총 {time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()