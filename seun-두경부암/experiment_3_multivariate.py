# -*- coding: utf-8 -*-
"""
실험 3: Multivariate Cox (S2 대체 + S3 증분) — 단독 실행 스크립트
------------------------------------------------------------------
모델: A(임상·병리 X만) / B(A+ajcc8th_STAGE) / C(A+mStage)
지표: C-index(bootstrap CI + optimism) / 반복 5-fold CV(paired B vs C) / NRI/IDI
S2 판정: CV C-index(C) > CV C-index(B), 같은 fold paired
S3 판정: NRI(C vs A) > NRI(B vs A), IDI 동일 (증분 비교)
"""
import os
import time
import multiprocessing
import numpy as np
import pandas as pd
from experiment_utils import fast_cindex as _ci

from experiment_utils import (X_MULTI_BASE, build_multi_df, fit_multi)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")
os.makedirs(OUT_DIR, exist_ok=True)

Y_LABELS = ["PFS", "DSS", "LRRFS"]
SEEDS = [42, 123, 2026, 777]
PRIMARY_SEED = 42
N_BOOT = {PRIMARY_SEED: 1_000, 123: 300, 2026: 300, 777: 300}
NRI_TIME = 24.0


def model_cindex_bootstrap(dfx, covs, n_boot, seed):
    """multivariate 모델 C-index: bootstrap CI + optimism correction."""
    T = dfx["time"].values.astype(float)
    E = dfx["event"].values.astype(bool)
    n = len(dfx)
    rng = np.random.RandomState(seed)
    cph_app = fit_multi(dfx, covs)
    if cph_app is None:
        return np.nan, np.nan, np.nan, np.nan
    C_app = _ci(cph_app.predict_partial_hazard(dfx[covs]).values.ravel(), T, E)
    c_orig_all = np.full(n_boot, np.nan)
    optim = np.full(n_boot, np.nan)
    for b in range(n_boot):
        idx = rng.randint(0, n, n)
        boot = dfx.iloc[idx]
        try:
            cph = fit_multi(boot, covs)
            if cph is None:
                continue
            risk_orig = cph.predict_partial_hazard(dfx[covs]).values.ravel()
            c_orig_all[b] = _ci(risk_orig, T, E)
            risk_boot = cph.predict_partial_hazard(boot[covs]).values.ravel()
            C_boot = _ci(risk_boot, boot["time"].values,
                         boot["event"].values.astype(bool))
            optim[b] = C_boot - c_orig_all[b]
        except Exception:
            continue
    c_orig_all = c_orig_all[~np.isnan(c_orig_all)]
    optim = optim[~np.isnan(optim)]
    lo, hi = np.nanpercentile(c_orig_all, [2.5, 97.5])
    corr = C_app - (float(np.mean(optim)) if len(optim) else 0.0)
    return C_app, lo, hi, corr


def repeated_cv_cindex(dfx, covs, k=5, repeats=10, seed=42):
    """반복 k-fold CV → held-out C-index 평균."""
    T = dfx["time"].values.astype(float)
    E = dfx["event"].values.astype(bool)
    n = len(dfx)
    rng = np.random.RandomState(seed)
    scores = []
    for rep in range(repeats):
        perm = rng.permutation(n)
        folds = np.array_split(perm, k)
        for f in range(k):
            test_idx = folds[f]
            train_idx = np.concatenate([folds[j] for j in range(k) if j != f])
            try:
                cph = fit_multi(dfx.iloc[train_idx], covs)
                if cph is None:
                    continue
                risk = cph.predict_partial_hazard(dfx.iloc[test_idx][covs]).values.ravel()
                scores.append(_ci(risk, T[test_idx], E[test_idx]))
            except Exception:
                continue
    if not scores:
        return np.nan, np.nan
    return float(np.mean(scores)), float(np.std(scores))


def repeated_cv_cindex_paired(dfx, covs_a, covs_b, k=5, repeats=10, seed=42):
    """S2: 같은 fold에서 두 모델 fit → ΔCV = C_a - C_b 분포."""
    T = dfx["time"].values.astype(float)
    E = dfx["event"].values.astype(bool)
    n = len(dfx)
    rng = np.random.RandomState(seed)
    sa, sb = [], []
    for rep in range(repeats):
        perm = rng.permutation(n)
        folds = np.array_split(perm, k)
        for f in range(k):
            test_idx = folds[f]
            train_idx = np.concatenate([folds[j] for j in range(k) if j != f])
            tr, te = dfx.iloc[train_idx], dfx.iloc[test_idx]
            try:
                cph_a = fit_multi(tr, covs_a)
                ca = _ci(cph_a.predict_partial_hazard(te[covs_a]).values.ravel(),
                         T[test_idx], E[test_idx]) if cph_a is not None else np.nan
            except Exception:
                ca = np.nan
            try:
                cph_b = fit_multi(tr, covs_b)
                cb = _ci(cph_b.predict_partial_hazard(te[covs_b]).values.ravel(),
                         T[test_idx], E[test_idx]) if cph_b is not None else np.nan
            except Exception:
                cb = np.nan
            if not (np.isnan(ca) or np.isnan(cb)):
                sa.append(ca)
                sb.append(cb)
    if not sa:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
    sa, sb = np.array(sa), np.array(sb)
    d = sb - sa  # ΔCV = C(covs_b) - C(covs_a): 양수 = covs_b 우월
    return (float(sa.mean()), float(sa.std()),
            float(sb.mean()), float(sb.std()),
            float(d.mean()), float(np.percentile(d, 2.5)),
            float(np.percentile(d, 97.5)))


def nri_idi_multi(dfx, covs_a, covs_b, t=NRI_TIME):
    """중첩 모델 NRI/IDI. 같은 데이터에서 두 모델 예측 비교."""
    cph_a = fit_multi(dfx, covs_a)
    cph_b = fit_multi(dfx, covs_b)
    if cph_a is None or cph_b is None:
        return np.nan, np.nan
    sa = cph_a.predict_survival_function(dfx[covs_a], times=[t]).values[0].ravel()
    sb = cph_b.predict_survival_function(dfx[covs_b], times=[t]).values[0].ravel()
    ra, rb = 1.0 - sa, 1.0 - sb
    Ts = dfx["time"].values.astype(float)
    Es = dfx["event"].values.astype(bool)
    case = Es & (Ts <= t)
    ctrl = Ts > t
    if case.sum() == 0 or ctrl.sum() == 0:
        return np.nan, np.nan
    up, dn = rb > ra, rb < ra
    nri = (up[case].mean() - dn[case].mean()) - (up[ctrl].mean() - dn[ctrl].mean())
    idi = (rb[case].mean() - ra[case].mean()) - (rb[ctrl].mean() - ra[ctrl].mean())
    return nri, idi


def analyze_multivariate(args):
    label, seed = args
    t0 = time.time()
    m = build_multi_df(label)
    nb = N_BOOT[seed]
    covs = X_MULTI_BASE

    rows = []
    model_defs = [
        ("A_clinical", covs),
        ("B_ajcc8th", covs + ["ajcc8th_STAGE"]),
        ("C_modified", covs + ["mStage"]),
    ]
    for tag, covs_m in model_defs:
        C_app, lo, hi, corr = model_cindex_bootstrap(m, covs_m, nb, seed)
        cv_mean, cv_sd = repeated_cv_cindex(m, covs_m, k=5, repeats=10, seed=seed)
        rows.append({"label": label, "seed": seed, "model": tag, "n": len(m),
                     "C_app": C_app, "C_lo": lo, "C_hi": hi,
                     "C_optimism_corr": corr, "CV_mean": cv_mean, "CV_sd": cv_sd})

    covs_a = covs
    covs_b = covs + ["ajcc8th_STAGE"]
    covs_c = covs + ["mStage"]
    nri_ba, idi_ba = nri_idi_multi(m, covs_a, covs_b)
    nri_ca, idi_ca = nri_idi_multi(m, covs_a, covs_c)
    nri_cb, idi_cb = nri_idi_multi(m, covs_b, covs_c)
    for tag, nri, idi in [("B_vs_A", nri_ba, idi_ba),
                          ("C_vs_A", nri_ca, idi_ca),
                          ("C_vs_B", nri_cb, idi_cb)]:
        rows.append({"label": label, "seed": seed, "model": tag, "n": len(m),
                     "C_app": np.nan, "C_lo": np.nan, "C_hi": np.nan,
                     "C_optimism_corr": np.nan, "CV_mean": np.nan, "CV_sd": np.nan,
                     "NRI": nri, "IDI": idi})

    # S2(대체): 같은 CV fold에서 B vs C paired 비교
    cvB_m, cvB_sd, cvC_m, cvC_sd, d_cv, d_lo, d_hi = \
        repeated_cv_cindex_paired(m, covs_b, covs_c, k=5, repeats=10, seed=seed)
    rows.append({"label": label, "seed": seed, "model": "S2_C_vs_B", "n": len(m),
                 "CV_B_mean": cvB_m, "CV_B_sd": cvB_sd,
                 "CV_C_mean": cvC_m, "CV_C_sd": cvC_sd,
                 "CV_d_mean": d_cv, "CV_d_lo": d_lo, "CV_d_hi": d_hi,
                 "NRI": nri_cb, "IDI": idi_cb})

    # S3(증분 비교): NRI(C vs A) > NRI(B vs A), IDI 동일
    rows.append({"label": label, "seed": seed, "model": "S3_increment_C_vs_B",
                 "n": len(m),
                 "NRI_C_vs_A": nri_ca, "IDI_C_vs_A": idi_ca,
                 "NRI_B_vs_A": nri_ba, "IDI_B_vs_A": idi_ba,
                 "NRI_diff": (nri_ca - nri_ba) if not (np.isnan(nri_ca) or np.isnan(nri_ba)) else np.nan,
                 "IDI_diff": (idi_ca - idi_ba) if not (np.isnan(idi_ca) or np.isnan(idi_ba)) else np.nan})

    # HR 테이블 (full-data fit, PRIMARY_SEED 행에만)
    if seed == PRIMARY_SEED:
        cph_full = fit_multi(m, covs)
        if cph_full is not None:
            summ = cph_full.summary
            for v in covs:
                ci = cph_full.confidence_intervals_.loc[v]
                rows.append({"label": label, "seed": seed, "model": f"HR_{v}", "n": len(m),
                             "HR": np.exp(cph_full.params_[v]),
                             "HR_lo": np.exp(ci["95% lower-bound"]),
                             "HR_hi": np.exp(ci["95% upper-bound"]),
                             "HR_p": summ.loc[v, "p"]})

    print(f"[실험3] {label} | seed={seed} | n={len(m)} | {time.time()-t0:.0f}s", flush=True)
    return rows


def main():
    start = time.time()
    try:
        multiprocessing.set_start_method("fork", force=True)
    except RuntimeError:
        pass
    tasks = [(l, s) for l in Y_LABELS for s in SEEDS]
    nproc = min(8, multiprocessing.cpu_count())
    print(f"실험 3 (S2/S3 multivariate) | {len(tasks)} 콤보 | Pool {nproc}, fork", flush=True)

    with multiprocessing.Pool(processes=nproc) as pool:
        rows_flat = pool.map(analyze_multivariate, tasks)

    rows = [r for chunk in rows_flat for r in chunk]
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT_DIR, "exp1_multivariate.csv"), index=False,
               encoding="utf-8-sig")
    print(f"[완료] results/exp1/exp1_multivariate.csv (총 {time.time()-start:.0f}s)", flush=True)


if __name__ == "__main__":
    main()