# -*- coding: utf-8 -*-
"""
실험 4e: 생존시간 회귀 (Survival-time regression) — classification 대체
------------------------------------------------------------------
24개월 이진화의 임의성을 없애고, **생존기간(월)을 회귀로 예측**.

모델 (5-fold × 4시드):
  AFT_LogNormal / AFT_Weibull  : lifelines — 중도절단 정식 처리 (로그선형)
  Lasso / Ridge / DecisionTree : naive(보정 없음) + IPCW(1/Ĝ(T) 가중)
  TabICL_naive / TabICL_uncensored : TabICLRegressor — sample_weight 미지원
      → naive(전체 log-time) / uncensored만(IPCW 근사)

평가 (이진화 없음):
  C-index (주, fast_cindex — 중도절단 보정 순위 일치)
  Spearman ρ (uncensored만 — 예측 log-time vs 실제 log-time)
  RMSE (log-time, uncensored만)

출력:
  exp1_surv_regression.csv     : 방식×y×시드 성능
  exp1_surv_regression_summary.csv : 방식×y 평균
"""
import os
import sys
import time
import subprocess
import tempfile
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from lifelines import KaplanMeierFitter

from experiment_utils import fast_cindex

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_CSV = os.path.join(BASE_DIR, "preprocessed_data.csv")
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")

SEEDS = [42, 123, 2026, 777]
N_FOLDS = 5
N_EXPAND = 100
Y_LABELS = ["PFS", "DSS", "LRRFS"]
# 표기: {모델}_{norm}_{naive|IPCW}
#   norm: std (StandardScaler) / power (PowerTransformer, Yeo-Johnson)
#   TabICL은 내부 처리 → normalization 없이 naive/uncensored만
METHODS = ["AFT_LogNormal", "AFT_Weibull",
           "Lasso_std_naive", "Lasso_std_IPCW", "Lasso_power_naive", "Lasso_power_IPCW",
           "Ridge_std_naive", "Ridge_std_IPCW", "Ridge_power_naive", "Ridge_power_IPCW",
           "DecisionTree_std_naive", "DecisionTree_std_IPCW",
           "DecisionTree_power_naive", "DecisionTree_power_IPCW",
           "TabICL_naive", "TabICL_uncensored"]

COVS = ["age", "성별", "tumor size (cm)", "DOI (mm)", "differentiation",
        "budding_01vs23", "PNI", "LVI", "RM", "TIL", "TSR",
        "WPOI5_2tier", "HPV/P16_1", "HPV/P16_2", "CCRT_bin"]


def load_prep():
    d = pd.read_csv(DATA_CSV)
    flag = [c for c in d.columns if c.endswith("_known")]
    oh = [c for c in d.columns if c.startswith(("subsite_", "Tx_3tier_", "HPV/P16_"))]
    X = list(dict.fromkeys(COVS + flag + oh))
    return d, X


def fold_split(n, event, seed, n_folds=N_FOLDS):
    rng = np.random.RandomState(seed)
    pos = np.nonzero(event == 1)[0]
    neg = np.nonzero(event == 0)[0]
    rng.shuffle(pos)
    rng.shuffle(neg)
    pf = np.array_split(pos, n_folds)
    nf = np.array_split(neg, n_folds)
    folds = []
    for f in range(n_folds):
        test_idx = np.concatenate([pf[f], nf[f]])
        train_idx = np.concatenate([
            np.concatenate([pf[j] for j in range(n_folds) if j != f]),
            np.concatenate([nf[j] for j in range(n_folds) if j != f])])
        folds.append((train_idx, test_idx))
    return folds


def censoring_km(T, E):
    km = KaplanMeierFitter()
    km.fit(T, event_observed=(1 - E))
    return km


def _ipcw_weights(T, E):
    km_c = censoring_km(T, E)
    g = km_c.survival_function_at_times(pd.Series(T)).values.ravel()
    w = np.where(E == 1, 1.0 / np.clip(g, 1e-3, 1.0), 0.0)
    return w


def fit_predict(method, Xtr, Ttr, Etr, Xte, seed):
    """fit → test 예측 log-time. (sklearn/lifelines — 부모 프로세스 직접)"""
    from sklearn.preprocessing import StandardScaler, PowerTransformer
    from sklearn.linear_model import Lasso, Ridge
    from sklearn.tree import DecisionTreeRegressor

    parts = method.split("_")
    base = parts[0]
    norm = parts[1] if len(parts) >= 2 and parts[1] in ("std", "power") else None
    weight_tag = parts[2] if len(parts) >= 3 else parts[1]

    if norm == "power":
        sc = PowerTransformer(method="yeo-johnson")
    else:
        sc = StandardScaler()
    Xtr_s = sc.fit_transform(Xtr)
    Xte_s = sc.transform(Xte)

    y_tr = np.log(np.clip(Ttr, 0.5, None))
    w = _ipcw_weights(Ttr, Etr) if weight_tag == "IPCW" else np.ones(len(Ttr))
    use_unc = method.endswith("_uncensored")
    if use_unc:
        mask = Etr == 1
        Xtr_s, y_tr, w = Xtr_s[mask], y_tr[mask], w[mask]

    if method == "AFT_LogNormal":
        from lifelines import LogNormalAFTFitter
        df = pd.DataFrame(Xtr_s, columns=[f"x{i}" for i in range(Xtr_s.shape[1])])
        df["T"] = Ttr
        df["E"] = Etr.astype(bool)
        m = LogNormalAFTFitter(penalizer=0.1)
        m.fit(df, duration_col="T", event_col="E")
        pred = m.predict_median(pd.DataFrame(Xte_s,
                                             columns=[f"x{i}" for i in range(Xte_s.shape[1])]))
        return np.log(np.clip(pred.values.ravel(), 0.5, None))

    if method == "AFT_Weibull":
        from lifelines import WeibullAFTFitter
        df = pd.DataFrame(Xtr_s, columns=[f"x{i}" for i in range(Xtr_s.shape[1])])
        df["T"] = Ttr
        df["E"] = Etr.astype(bool)
        m = WeibullAFTFitter(penalizer=0.1)
        m.fit(df, duration_col="T", event_col="E")
        pred = m.predict_median(pd.DataFrame(Xte_s,
                                             columns=[f"x{i}" for i in range(Xte_s.shape[1])]))
        return np.log(np.clip(pred.values.ravel(), 0.5, None))

    if base == "Lasso":
        m = Lasso(alpha=0.05, max_iter=5000, random_state=seed)
    elif base == "Ridge":
        m = Ridge(alpha=1.0, random_state=seed)
    elif base == "DecisionTree":
        m = DecisionTreeRegressor(max_depth=4, min_samples_leaf=5, random_state=seed)
    else:
        raise ValueError(method)

    m.fit(Xtr_s, y_tr, sample_weight=w if w.sum() > 0 else None)
    return m.predict(Xte_s)


def _tabicl_reg_subprocess(args):
    """TabICLRegressor fit → test log-time. (subprocess 격리 — segfault 방지)"""
    Xtr, ytr, Xte, seed, out_prefix = args
    script = f'''
import os, sys, numpy as np
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, {BASE_DIR!r})
from tabicl import TabICLRegressor
Xtr = np.load({out_prefix!r} + "_X.npy")
ytr = np.load({out_prefix!r} + "_y.npy")
Xte = np.load({out_prefix!r} + "_test.npy")
m = TabICLRegressor(n_estimators=1, random_state={seed})
m.fit(Xtr, ytr)
pred = np.asarray(m.predict(Xte, output_type="mean")).ravel()
np.save({out_prefix!r} + "_pred.npy", pred)
os._exit(0)
'''
    np.save(out_prefix + "_X.npy", Xtr)
    np.save(out_prefix + "_y.npy", ytr)
    np.save(out_prefix + "_test.npy", Xte)
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
        pred_path = out_prefix + "_pred.npy"
        if os.path.exists(pred_path):
            return np.load(pred_path)
        return None
    except Exception:
        return None
    finally:
        for suffix in ("_X.npy", "_y.npy", "_test.npy", "_pred.npy"):
            try:
                os.remove(out_prefix + suffix)
            except OSError:
                pass


def fit_predict_tabicl(Xtr, Ttr, Etr, Xte, seed, use_unc):
    y_tr = np.log(np.clip(Ttr, 0.5, None))
    if use_unc:
        mask = Etr == 1
        Xtr, y_tr = Xtr[mask], y_tr[mask]
    tmpdir = tempfile.mkdtemp(prefix="tabicl_reg_")
    out_prefix = os.path.join(tmpdir, f"r_{seed}")
    pred = _tabicl_reg_subprocess((Xtr, y_tr, Xte, seed, out_prefix))
    if pred is None:
        return None
    return np.log(np.clip(pred, 0.5, None))


def main():
    t0 = time.time()
    d, X = load_prep()
    print(f"생존시간 회귀 | X {len(X)}개, 모델 {len(METHODS)}종, "
          f"시드 {len(SEEDS)} × fold {N_FOLDS} × y {len(Y_LABELS)}", flush=True)

    ckpt_csv = os.path.join(OUT_DIR, "exp1_surv_regression.csv")
    done = set()
    if os.path.exists(ckpt_csv):
        prev = pd.read_csv(ckpt_csv)
        done = set(zip(prev["method"], prev["y"], prev["seed"]))
        print(f"체크포인트: {len(done)}개 완료", flush=True)
    rows = []
    if os.path.exists(ckpt_csv):
        rows = pd.read_csv(ckpt_csv).to_dict("records")

    for y in Y_LABELS:
        y_time, y_event = f"{y}_time", f"{y}_event"
        dfx = d[X + [y_time, y_event]].dropna().reset_index(drop=True)
        X_all = dfx[X].values.astype(float)
        T_all = dfx[y_time].values.astype(float)
        E_all = dfx[y_event].values.astype(int)

        for method in METHODS:
            is_tabicl = method.startswith("TabICL")
            use_unc = method.endswith("_uncensored")
            for seed in SEEDS:
                key = (method, y, seed)
                if key in done:
                    continue
                folds = fold_split(len(dfx), E_all, seed)
                scores = np.full(len(dfx), np.nan)
                for f, (train_idx, test_idx) in enumerate(folds):
                    rng_try = np.random.RandomState(seed * 1000 + f * 10)
                    boot_idx = rng_try.choice(train_idx, size=N_EXPAND, replace=True)
                    try:
                        if is_tabicl:
                            pred = fit_predict_tabicl(
                                X_all[boot_idx], T_all[boot_idx], E_all[boot_idx],
                                X_all[test_idx], seed * 100 + f * 10, use_unc)
                        else:
                            pred = fit_predict(
                                method, X_all[boot_idx], T_all[boot_idx],
                                E_all[boot_idx], X_all[test_idx],
                                seed * 100 + f * 10)
                        if pred is not None and len(pred) == len(test_idx):
                            scores[test_idx] = pred
                    except Exception as e:
                        print(f"  [{y}] {method} seed={seed} fold{f+1} 실패: {e}",
                              flush=True)
                valid = ~np.isnan(scores)
                if valid.sum() < 20:
                    continue
                risk = -scores[valid]
                c_idx = float(fast_cindex(risk, T_all[valid],
                                          E_all[valid].astype(bool)))
                unc = valid & (E_all == 1)
                rho, rmse = np.nan, np.nan
                if unc.sum() >= 5:
                    rho = float(spearmanr(scores[unc], np.log(T_all[unc]))[0])
                    rmse = float(np.sqrt(np.mean((scores[unc] - np.log(T_all[unc])) ** 2)))
                rows.append({"method": method, "y": y, "seed": seed,
                             "n_valid": int(valid.sum()), "n_unc": int(unc.sum()),
                             "C_index": round(c_idx, 4),
                             "spearman_rho": round(rho, 4),
                             "RMSE_logtime": round(rmse, 4)})
                pd.DataFrame(rows).to_csv(ckpt_csv, index=False, encoding="utf-8-sig")
                done.add(key)
                print(f"  [{y}] {method} seed={seed}: C={c_idx:.3f} ρ={rho:+.3f} "
                      f"({time.time()-t0:.0f}s)", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(ckpt_csv, index=False, encoding="utf-8-sig")
    summ = res.groupby(["method", "y"])[["C_index", "spearman_rho", "RMSE_logtime"]].mean().round(3)
    summ.to_csv(os.path.join(OUT_DIR, "exp1_surv_regression_summary.csv"),
                encoding="utf-8-sig")

    print("\n=== 생존시간 회귀 요약 (시드 평균) ===")
    print(summ.to_string())
    print(f"\n[완료] {ckpt_csv} (총 {time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()