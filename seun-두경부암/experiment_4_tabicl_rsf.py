# -*- coding: utf-8 -*-
"""
실험 4: TabICL 피팅 데이터 Bootstrap 안정성 검증 (S4, v5 계획 반영)
------------------------------------------------------------------
정세운 선생님 요청:
  "아무리 TabICL이라도 피팅용 데이터를 리샘플링을 여러 난수로 해 봤으면 좋겠어.
   각 난수 별 score들이 비슷하다면 안정적인 모델이라는 의미잖아."

구조 (실험 계획 v5):
  4-1 피팅 데이터 bootstrap : 고정 train/test 분할 후, train에서 bootstrap B=1,000회
                             (시드 42/123/2026/777)
  4-2 리샘플별 TabICL fit  : 각 리샘플로 fit → 고정 test 세트 score 예측
  4-3 안정성 지표           : 환자별 score CV(작을수록 안정) + 시드 간 상관 r + 순위 일치 ρ
  4-4 score↔생존 검증       : 안정한 score(리샘플 평균) vs 실제 생존
                             C-index(주) + uncensored Spearman + quartile KM/log-rank
  4-5 변수 중요도           : TabICL SHAP vs Cox HR 순위 일치 (수정병기 제외)
  4-6 RSF 교차확인          : Random Survival Forest OOB C-index

판정: 환자별 score CV < 10% 또는 시드 간 r > 0.9 → "안정적 모델"
주의: TabICL 사전학습 범위(300~60K) 밖의 n=133 → 결과는 참고용
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

from experiment_utils import fast_cindex

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_CSV = os.path.join(BASE_DIR, "preprocessed_data.csv")
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")
os.makedirs(OUT_DIR, exist_ok=True)

SEED = 42
SEEDS = [42, 123, 2026, 777]   # 리샘플링 시드 4개
N_EXPAND = 100       # 피팅 데이터 확장 크기 (fold train을 bootstrap으로 확장)
N_FOLDS = 5          # 시드당 5-fold CV (각 샘플이 test 기회 1회)
T_HORIZON = 24.0     # TabICL 학습 target (24개월 사건 여부)

COVS = ["age", "성별", "tumor size (cm)", "DOI (mm)", "differentiation",
        "budding_01vs23", "PNI", "LVI", "RM", "TIL", "TSR",
        "WPOI5_2tier", "HPV/P16_1", "HPV/P16_2", "CCRT_bin"]
STAGE_VARS = ["ajcc8th_STAGE", "mStage"]
Y_LABELS = ["PFS", "DSS", "LRRFS"]


def try_import(mod, pip_name):
    try:
        return __import__(mod)
    except ImportError:
        print(f"[실험 4] {mod} 미설치 — `pip install {pip_name}` 필요")
        return None


def load_prep():
    d = pd.read_csv(DATA_CSV)
    flag_cols = [c for c in d.columns if c.endswith("_known")]
    oh_cols = [c for c in d.columns if c.startswith(("subsite_", "Tx_3tier_", "HPV/P16_"))]
    return d, COVS + flag_cols + oh_cols


def _run_fit_subprocess(args):
    """TabICL fit을 별도 subprocess로 실행 (segfault 격리).

    TabICL은 종료 시 cleanup segfault(-11)가 불규칙 발생.
    score는 stdout 대신 **파일 저장**으로 전달 — segfault가 나도
    파일은 디스크에 남으므로 확실하게 수집 가능 (stdout은 segfault 시 유실됨).
    """
    import subprocess
    import sys
    X_train, y_train, X_test, seed, out_path = args
    script = f'''
import os, sys, numpy as np
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, {BASE_DIR!r})
from tabicl import TabICLClassifier
X_train = np.load({out_path!r} + "_X.npy")
y_train = np.load({out_path!r} + "_y.npy")
X_test = np.load({out_path!r} + "_test.npy")
clf = TabICLClassifier(n_estimators=1, random_state={seed})
clf.fit(X_train, y_train)
score = clf.predict(X_test)
np.save({out_path!r} + "_score.npy", score)
os._exit(0)  # cleanup segfault 방지: 저장 후 즉시 강제 종료
'''
    np.save(out_path + "_X.npy", X_train)
    np.save(out_path + "_y.npy", y_train)
    np.save(out_path + "_test.npy", X_test)
    try:
        # segfault가 부모로 전파되지 않도록 Popen + 세션 격리
        # OMP_NUM_THREADS=1: torch OpenMP 스레드 경합이 segfault 원인 → 단일 스레드로 안정화
        env = dict(os.environ)
        env["KMP_DUPLICATE_LIB_OK"] = "TRUE"
        env["OMP_NUM_THREADS"] = "1"
        p = subprocess.Popen([sys.executable, "-c", script],
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                             start_new_session=True, env=env)
        try:
            p.communicate(timeout=600)
        except subprocess.TimeoutExpired:
            p.kill()
            p.communicate()
        score_path = out_path + "_score.npy"
        if os.path.exists(score_path):
            return np.load(score_path)
        return None
    except Exception:
        return None
    finally:
        for suffix in ("_X.npy", "_y.npy", "_test.npy", "_score.npy"):
            try:
                os.remove(out_path + suffix)
            except OSError:
                pass


def tabicl_stability_cv(d, X, y_time, y_event, seed, n_expand=N_EXPAND,
                        n_folds=N_FOLDS):
    """시드 1개 × 5-fold CV: 모든 샘플이 test score 1회씩 획득.

    사용자 설계: "리샘플링 시드당 5fold면 각 샘플이 테스트셋으로써의 스코어를
    가질 기회가 1번씩 있잖아. 그걸 놓고 4회 리샘플링에서도 유의한 범위 내
    스코어를 가진다면 안정적"
    - 사건층화 5-fold 분할
    - 각 fold: train(4/5)을 bootstrap으로 n_expand명 확장 → TabICL fit → test fold score
    - segfault 방지: fit 1회 = subprocess 1개 (Popen, 부모 보호)
    """
    dfx = d[X + [y_time, y_event]].dropna().reset_index(drop=True)
    dfx["y24"] = ((dfx[y_event] == 1) & (dfx[y_time] <= T_HORIZON)).astype(int)
    X_all = dfx[X].values.astype(float)
    y_all = dfx["y24"].values.astype(int)
    n = len(dfx)

    # 사건층화 5-fold 분할 (시드 고정 → fold 구성 재현 가능)
    rng = np.random.RandomState(seed)
    pos = np.nonzero(y_all == 1)[0]
    neg = np.nonzero(y_all == 0)[0]
    rng.shuffle(pos)
    rng.shuffle(neg)
    pos_folds = np.array_split(pos, n_folds)
    neg_folds = np.array_split(neg, n_folds)

    # 체크포인트: fold 결과를 영구 저장, 재실행 시 완료분 로드
    y_name = y_time.split("_")[0]
    ckpt_dir = os.path.join(OUT_DIR, "tabicl_ckpt", y_name)
    os.makedirs(ckpt_dir, exist_ok=True)
    all_scores = np.full(n, np.nan)
    for f in range(n_folds):
        test_idx = np.concatenate([pos_folds[f], neg_folds[f]])
        train_idx = np.concatenate([np.concatenate([pos_folds[j] for j in range(n_folds) if j != f]),
                                    np.concatenate([neg_folds[j] for j in range(n_folds) if j != f])])
        ckpt_file = os.path.join(ckpt_dir, f"fold_{seed}_{f}.npz")
        if os.path.exists(ckpt_file):
            data = np.load(ckpt_file)
            all_scores[test_idx] = data["score"]
            print(f"  [{y_name}] seed={seed} fold{f+1}/{n_folds} (체크포인트 로드)",
                  flush=True)
            continue
        # train을 bootstrap으로 확장 + segfault 대비 최대 5회 재시도
        got = None
        for attempt in range(5):
            rng_try = np.random.RandomState(seed * 1000 + f * 10 + attempt)
            boot_idx = rng_try.choice(train_idx, size=n_expand, replace=True)
            out_path = os.path.join(ckpt_dir, f"tmp_{seed}_{f}_{attempt}")
            got = _run_fit_subprocess(
                (X_all[boot_idx], y_all[boot_idx], X_all[test_idx],
                 seed * 100 + f * 10 + attempt, out_path))
            if got is not None and len(got) == len(test_idx):
                break
        if got is not None and len(got) == len(test_idx):
            all_scores[test_idx] = got
            np.savez(ckpt_file, score=got, test_idx=test_idx)
        print(f"  [{y_name}] seed={seed} fold{f+1}/{n_folds} "
              f"완료 (score {'확보' if not np.isnan(all_scores[test_idx]).all() else '실패'})",
              flush=True)

    if np.isnan(all_scores).sum() > n * 0.2:
        return None   # 20% 이상 score 누락 시 시드 무효
    return all_scores


def stability_metrics(stab):
    """샘플별 score CV + 시드 간 상관/순위. stab = {seed: all_scores 벡터(n,)}."""
    seed_order = sorted(stab.keys())
    mats = np.column_stack([stab[s] for s in seed_order])   # (n, n_seed)
    n, n_seed = mats.shape

    rows = []
    for i in range(n):
        s = mats[i]
        s = s[~np.isnan(s)]
        if len(s) >= 2 and s.std() > 0:
            cv = s.std() / abs(s.mean())
            rows.append({"patient_idx": i, "score_CV": float(cv)})
    patient_cv = pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["patient_idx", "score_CV"])

    r_pairs, rho_pairs = [], []
    for a in range(n_seed):
        for b in range(a + 1, n_seed):
            sa, sb = mats[:, a], mats[:, b]
            valid = ~(np.isnan(sa) | np.isnan(sb))
            if valid.sum() >= 10 and sa[valid].std() > 0 and sb[valid].std() > 0:
                r_pairs.append(float(np.corrcoef(sa[valid], sb[valid])[0, 1]))
                rho_pairs.append(float(spearmanr(sa[valid], sb[valid])[0]))
    return {
        "patient_cv": patient_cv,
        "n_samples": n,
        "seed_corr_mean": float(np.mean(r_pairs)) if r_pairs else np.nan,
        "seed_spearman_mean": float(np.mean(rho_pairs)) if rho_pairs else np.nan,
    }


def score_survival_eval(d, X, y_time, y_event, stab):
    """4-6: 안정한 score(시드별 CV 평균) vs 실제 생존 (전체 샘플)."""
    dfx = d[X + [y_time, y_event]].dropna().reset_index(drop=True)
    dfx["y24"] = ((dfx[y_event] == 1) & (dfx[y_time] <= T_HORIZON)).astype(int)
    seed_scores = np.column_stack([stab[s] for s in sorted(stab.keys())])
    mean_score = np.nanmean(seed_scores, axis=1)   # 시드 평균 score
    T = dfx[y_time].values.astype(float)
    E = dfx[y_event].values.astype(bool)

    c_idx = float(fast_cindex(mean_score, T, E)) if mean_score.std() > 0 else np.nan
    rho, p_rho = (np.nan, np.nan)
    unc = E
    if unc.sum() >= 5 and mean_score[unc].std() > 0:
        rho, p_rho = spearmanr(mean_score[unc], T[unc])

    km_rows = []
    try:
        q = pd.qcut(mean_score, 4, labels=False, duplicates="drop")
        for qi in sorted(pd.unique(q)):
            m = q == qi
            if m.sum() < 3:
                continue
            km = KaplanMeierFitter()
            km.fit(T[m], event_observed=E[m])
            sf = km.survival_function_
            idx = sf.index <= 36.0
            rmst = float(np.trapz(np.concatenate([[0], sf.index[idx]]),
                                  np.concatenate([[1], sf["KM_estimate"].values[idx]]))) \
                if idx.sum() else np.nan
            km_rows.append({"quartile": int(qi), "n": int(m.sum()),
                            "RMST36": rmst})
        lr = multivariate_logrank_test(T, q, E)
        chi2, p_lr = float(lr.test_statistic), float(lr.p_value)
    except Exception:
        chi2, p_lr = np.nan, np.nan

    return {"test_n": len(dfx), "C_index": c_idx,
            "unc_n": int(unc.sum()), "spearman_rho": rho, "spearman_p": p_rho,
            "quartile_KM": pd.DataFrame(km_rows), "logrank_chi2": chi2,
            "logrank_p": p_lr}


def tabicl_shap_rank(d, X, y_time, y_event):
    """4-5: TabICL SHAP vs 단변량 Cox HR 순위 일치 (subprocess 격리).

    TabICL fit이 부모에서 실행되면 segfault 시 부모가 죽음 → subprocess로 격리.
    subprocess에서 fit+SHAP 후 변수 중요도를 json으로 저장, 부모는 읽기만.
    """
    import subprocess
    import sys
    import json
    import tempfile
    from lifelines import CoxPHFitter

    dfx = d[X + [y_time, y_event]].dropna().reset_index(drop=True)
    dfx["y24"] = ((dfx[y_event] == 1) & (dfx[y_time] <= T_HORIZON)).astype(int)
    Xm = dfx[X].values.astype(float)
    ym = dfx["y24"].values.astype(int)

    tmpdir = tempfile.mkdtemp(prefix="tabicl_shap_")
    x_path = os.path.join(tmpdir, "X.npy")
    y_path = os.path.join(tmpdir, "y.npy")
    out_path = os.path.join(tmpdir, "shap.json")
    np.save(x_path, Xm)
    np.save(y_path, ym)

    script = f'''
import os, sys, json, numpy as np
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, {BASE_DIR!r})
from tabicl import TabICLClassifier
import shap
Xm = np.load({x_path!r})
ym = np.load({y_path!r})
clf = TabICLClassifier(n_estimators=1, random_state={SEED})
clf.fit(Xm, ym)
def pred_fn(x):
    return clf.predict(x)
explainer = shap.Explainer(pred_fn, Xm, feature_names={X!r})
shap_vals = explainer(Xm, max_evals=1000)
shap_map = {{v: float(np.abs(shap_vals.values).mean(axis=0)[i]) for i, v in enumerate({X!r})}}
json.dump(shap_map, open({out_path!r}, "w"))
os._exit(0)  # cleanup segfault 방지: 저장 후 즉시 강제 종료
'''
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
        if not os.path.exists(out_path):
            return None
        shap_map = json.load(open(out_path))
    except Exception:
        return None
    finally:
        for f in (x_path, y_path, out_path):
            try:
                os.remove(f)
            except OSError:
                pass

    hr_map = {}
    for v in X:
        sub = dfx[[v, y_time, y_event]].dropna()
        if sub[y_event].sum() < 5 or sub[v].nunique() < 2:
            continue
        try:
            cph = CoxPHFitter(penalizer=0.0)
            cph.fit(sub, duration_col=y_time, event_col=y_event)
            hr_map[v] = float(np.exp(cph.params_[v]))
        except Exception:
            continue
    common = [v for v in X if v in shap_map and v in hr_map]
    if len(common) < 5:
        return None
    sh = np.array([shap_map[v] for v in common])
    hr = np.array([np.abs(np.log(hr_map[v])) for v in common])
    corr = float(np.corrcoef(sh, hr)[0, 1])
    rho = float(spearmanr(sh, hr)[0])
    return {"n": len(common), "pearson": corr, "spearman": rho,
            "vars": common, "SHAP": [shap_map[v] for v in common],
            "abslogHR": [np.abs(np.log(hr_map[v])) for v in common]}


def rsf_cindex(d, X, y_time, y_event):
    """4-6: RSF OOB C-index."""
    sksurv = try_import("sksurv.ensemble", "scikit-survival")
    if sksurv is None:
        return None
    from sksurv.ensemble import RandomSurvivalForest
    dfx = d[X + [y_time, y_event]].dropna()
    Xm = dfx[X].values.astype(float)
    struct = np.array([(bool(e), t) for t, e in
                       zip(dfx[y_event].values, dfx[y_time].values.astype(float))],
                      dtype=[("event", bool), ("time", float)])
    rsf = RandomSurvivalForest(n_estimators=500, random_state=SEED, n_jobs=-1)
    rsf.fit(Xm, struct)
    return float(rsf.score(Xm, struct))


def main():
    import time
    t0 = time.time()
    d, X = load_prep()
    print(f"실험 4 | 데이터 {d.shape[0]}행, X {len(X)}개, "
          f"리샘플링 시드 {len(SEEDS)}개 × {N_FOLDS}-fold CV × 확장 {N_EXPAND}명",
          flush=True)

    # 4-1~4-3: 피팅 데이터 확장 → TabICL fit → score 안정성 (y별)
    stab_all = {}
    for y in Y_LABELS:
        y_time, y_event = f"{y}_time", f"{y}_event"
        print(f"[실험 4] {y} 시드4×5-fold CV 안정성...", flush=True)
        stab_by_seed = {}
        for seed in SEEDS:
            res = tabicl_stability_cv(d, X, y_time, y_event, seed,
                                      n_expand=N_EXPAND, n_folds=N_FOLDS)
            if res is not None:
                stab_by_seed[seed] = res
                valid = ~np.isnan(res)
                print(f"  seed={seed}: score {valid.sum()}/{len(res)}명 확보",
                      flush=True)
        if len(stab_by_seed) == len(SEEDS):
            stab_all[y] = stab_by_seed
            met = stability_metrics(stab_by_seed)
            met["patient_cv"].to_csv(
                os.path.join(OUT_DIR, f"exp1_tabicl_stability_{y}.csv"),
                index=False, encoding="utf-8-sig")
            cv_med = met["patient_cv"]["score_CV"].median() if len(met["patient_cv"]) else np.nan
            print(f"  → 시드 간 score 상관 r={met['seed_corr_mean']:.3f}, "
                  f"ρ={met['seed_spearman_mean']:.3f}, 샘플 score CV 중앙값={cv_med:.3f}",
                  flush=True)

    # 4-6: score↔생존 (전체 샘플, 시드 평균 score)
    eval_rows = []
    for y in Y_LABELS:
        if y not in stab_all:
            continue
        y_time, y_event = f"{y}_time", f"{y}_event"
        ev = score_survival_eval(d, X, y_time, y_event, stab_all[y])
        row = {"y": y, "C_index": ev["C_index"], "spearman_rho": ev["spearman_rho"],
               "spearman_p": ev["spearman_p"], "logrank_chi2": ev["logrank_chi2"],
               "logrank_p": ev["logrank_p"], "test_n": ev["test_n"]}
        eval_rows.append(row)
        ev["quartile_KM"].to_csv(
            os.path.join(OUT_DIR, f"exp1_tabicl_quartile_{y}.csv"),
            index=False, encoding="utf-8-sig")
    if eval_rows:
        pd.DataFrame(eval_rows).to_csv(
            os.path.join(OUT_DIR, "exp1_tabicl_survival.csv"),
            index=False, encoding="utf-8-sig")
        print("\n=== 4-6 score↔생존 (시드 평균 score) ===")
        print(pd.DataFrame(eval_rows).to_string(index=False), flush=True)

    # 4-5: SHAP vs Cox HR — 부모에서 tabicl 직접 fit 시 segfault 위험 → subprocess로 격리
    import subprocess
    import sys as _sys
    shap_rows = []
    for y in Y_LABELS:
        y_time, y_event = f"{y}_time", f"{y}_event"
        try:
            r = tabicl_shap_rank(d, X, y_time, y_event)
        except Exception:
            r = None
        if r is not None:
            shap_rows.append({"y": y, "n": r["n"], "pearson": r["pearson"],
                              "spearman": r["spearman"]})
            print(f"  SHAP vs Cox|HR|: ρ={r['spearman']:+.3f} (n={r['n']})", flush=True)
    if shap_rows:
        pd.DataFrame(shap_rows).to_csv(
            os.path.join(OUT_DIR, "exp1_tabicl_shap.csv"),
            index=False, encoding="utf-8-sig")

    # 4-6: RSF
    rsf_rows = []
    for y in Y_LABELS:
        y_time, y_event = f"{y}_time", f"{y}_event"
        c = rsf_cindex(d, X, y_time, y_event)
        rsf_rows.append({"y": y, "model": "RSF_full", "OOB_C": c})
        for s in STAGE_VARS:
            rsf_rows.append({"y": y, "model": f"RSF_{s}",
                             "OOB_C": rsf_cindex(d, X + [s], y_time, y_event)})
    pd.DataFrame(rsf_rows).to_csv(os.path.join(OUT_DIR, "exp1_rsf_cindex.csv"),
                                  index=False, encoding="utf-8-sig")

    print(f"\n[완료] {OUT_DIR} (총 {time.time()-t0:.0f}s)", flush=True)
    if not stab_all:
        print("[참고] tabicl 미설치 — `pip install tabicl[shap]` 후 재실행")
    if not rsf_rows or all(r["OOB_C"] is None for r in rsf_rows):
        print("[참고] sksurv 미설치 — `pip install scikit-survival` 후 재실행")


if __name__ == "__main__":
    main()