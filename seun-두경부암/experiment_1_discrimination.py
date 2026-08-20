# -*- coding: utf-8 -*-
"""
실험 1: 수정병기 vs 기존병기 예후 구분력 비교 + multivariate 모델 검증
------------------------------------------------------------------
Part A (univariate stage 비교):
  비교 쌍: mStage vs ajcc8th_STAGE | mTstage vs T stage | mNstage vs N stage
  Y label: PFS / DSS / LRRFS
  지표: C-index(bootstrap CI, paired dC), optimism(cor.), KM+log-rank,
        NRI/IDI(24개월), time-dependent AUC(IPCW)
  seed 4개(42/123/2026/777)로 결론 동일성 확인

Part B (multivariate, X변수 -> 각 y):
  모델 A: 임상/병리 X만 | B: A + ajcc8th_STAGE | C: A + mStage
  지표: HR(95%CI,p) 테이블, 모델 C-index(bootstrap CI + optimism),
        반복 5-fold CV C-index, 중첩 모델 NRI/IDI
"""
import os
import time
import multiprocessing
from functools import lru_cache

import numpy as np
import pandas as pd

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.utils import concordance_index
from scipy.stats import spearmanr

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(BASE_DIR, "정박사님께 드릴 raw data.xlsx")
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")
os.makedirs(OUT_DIR, exist_ok=True)

LABELS = ["PFS", "DSS", "LRRFS"]
PAIRS = [
    ("mStage",   "ajcc8th_STAGE", "mStage vs ajcc8th"),
    (" mTstage", "T stage",       "mTstage vs T stage"),
    ("mNstage",  "N stage",       "mNstage vs N stage"),
]
SEEDS = [42, 123, 2026, 777]
PRIMARY_SEED = 42

N_BOOT_CINDEX = {PRIMARY_SEED: 10_000, 123: 3_000, 2026: 3_000, 777: 3_000}
N_BOOT_OPTIM = {PRIMARY_SEED: 500, 123: 200, 2026: 200, 777: 200}
N_BOOT_NRI   = {PRIMARY_SEED: 1_000, 123: 300, 2026: 300, 777: 300}

NRI_TIME = 24.0
AUC_TIMES = [12, 24, 36, 60]

COVARIATES = ["age", "성별", "tumor size (cm)", "DOI (mm)", "differentiation",
              "budding_01vs23", "PNI", "LVI", "RM", "TIL", "TSR",
              "WPOI5_2tier", "HPV_P16_pos", "CCRT_bin"]


# ----------------------------------------------------------------------------
# 데이터 로드 (lru_cache: 프로세스 내 1회만 파싱)
# ----------------------------------------------------------------------------
@lru_cache(maxsize=None)
def load_data(label):
    xl = pd.ExcelFile(DATA_PATH)
    df = xl.parse("Sheet2")
    time_col = f"{label}_month"
    keep = ["연구번호", "T stage", "N stage", "ajcc8th_STAGE",
            " mTstage", "mNstage", "mStage", "time", "event"]
    df = df.rename(columns={time_col: "time", label: "event"})
    df["HPV_P16_pos"] = (df["HPV/P16"] >= 1).astype(int)
    ccrt = pd.to_datetime(df["1st OP 후 CCRT 날짜"], errors="coerce")
    df["CCRT_bin"] = ccrt.notna().astype(int)
    df = df[keep + COVARIATES + ["HPV/P16", "1st OP 후 CCRT 날짜"]].copy()
    df["event"] = df["event"].astype(int)
    df["time"] = df["time"].astype(float)
    df["N stage"] = pd.to_numeric(df["N stage"], errors="coerce")
    for c in ["T stage", "ajcc8th_STAGE", " mTstage", "mNstage", "mStage"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def build_pair_df(df, new_col, old_col):
    out = pd.DataFrame({
        "stage_new": pd.to_numeric(df[new_col], errors="coerce"),
        "stage_old": pd.to_numeric(df[old_col], errors="coerce"),
        "time": df["time"].astype(float),
        "event": df["event"].astype(int),
    })
    return out.dropna(subset=["stage_new", "stage_old"]).reset_index(drop=True)


# ----------------------------------------------------------------------------
# 완전 벡터화 C-index (Harrell). 높은 stage = 나쁜 예후.
#   broadcast: (사건행 × 전체행) 행렬, tied time 제외
# ----------------------------------------------------------------------------
def fast_cindex(S, T, E):
    S = np.asarray(S, dtype=float)
    T = np.asarray(T, dtype=float)
    E = np.asarray(E, dtype=bool)
    ev = np.nonzero(E)[0]
    if len(ev) == 0:
        return np.nan
    Ti = T[ev]
    M = T[None, :] > Ti[:, None]
    if not M.any():
        return np.nan
    Si = S[ev][:, None]
    Sj = np.where(M, S[None, :], np.nan)
    conc = np.nansum(Si > Sj) + 0.5 * np.nansum(Si == Sj)
    return conc / M.sum()


def bootstrap_combined(pair_df, n_boot, seed):
    """단일 루프에서 C_new / C_old / paired dC 동시 계산."""
    Sn = pair_df["stage_new"].values.astype(float)
    So = pair_df["stage_old"].values.astype(float)
    T = pair_df["time"].values.astype(float)
    E = pair_df["event"].values.astype(bool)
    n = len(pair_df)
    rng = np.random.RandomState(seed)
    cn = np.empty(n_boot)
    co = np.empty(n_boot)
    d = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.randint(0, n, n)
        cn[b] = fast_cindex(Sn[idx], T[idx], E[idx])
        co[b] = fast_cindex(So[idx], T[idx], E[idx])
        d[b] = cn[b] - co[b]
    c_new = fast_cindex(Sn, T, E)
    c_old = fast_cindex(So, T, E)
    pct = lambda a: np.nanpercentile(a, [2.5, 97.5])
    return dict(
        C_new=c_new, C_new_lo=pct(cn)[0], C_new_hi=pct(cn)[1],
        C_old=c_old, C_old_lo=pct(co)[0], C_old_hi=pct(co)[1],
        dC=c_new - c_old, dC_lo=pct(d)[0], dC_hi=pct(d)[1],
        dC_excl_0=bool((pct(d)[0] > 0) or (pct(d)[1] < 0)),
    )


# ----------------------------------------------------------------------------
# Optimism correction (Harrell)
#   univariate 단일 서수 변수에서는 C-index가 beta에 무관 -> optimism ~ 0.
#   multivariate(Part B)에서 실질적 의미를 가짐.
# ----------------------------------------------------------------------------
def optimism_correct(df, stage_col, n_boot, seed):
    S = df[stage_col].values.astype(float)
    T = df["time"].values.astype(float)
    E = df["event"].values.astype(bool)
    n = len(df)
    rng = np.random.RandomState(seed)
    C_app = fast_cindex(S, T, E)
    dfc = df[["time", "event", stage_col]].rename(columns={stage_col: "stage"})
    optimism = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        boot = dfc.iloc[idx]
        try:
            cph = CoxPHFitter(penalizer=0.0)
            beta = cph.fit(boot, duration_col="time", event_col="event").params_["stage"]
            sign = 1.0 if beta >= 0 else -1.0
            C_boot = fast_cindex(sign * beta * boot["stage"].values,
                                 boot["time"].values, boot["event"].values.astype(bool))
            C_orig = fast_cindex(sign * beta * S, T, E)
            optimism.append(C_boot - C_orig)
        except Exception:
            continue
    if not optimism:
        return np.nan
    return C_app - float(np.mean(optimism))


# ----------------------------------------------------------------------------
# NRI / IDI (24개월, continuous). cases=event<=t, ctrl=time>t
# ----------------------------------------------------------------------------
def _risk_at_t(df, stage_col, t):
    dfc = pd.DataFrame({"time": df["time"].values, "event": df["event"].values,
                        "stage": df[stage_col].values.astype(float)})
    cph = CoxPHFitter(penalizer=0.0)
    cph.fit(dfc, duration_col="time", event_col="event")
    sf = cph.predict_survival_function(dfc[["stage"]], times=[t])
    return (1.0 - sf.values[0]).ravel()


def nri_idi(df, new_col, old_col, t=NRI_TIME, n_boot=100, seed=42):
    T = df["time"].values.astype(float)
    E = df["event"].values.astype(bool)
    n = len(df)

    def comp(sub):
        rn = _risk_at_t(sub, new_col, t)
        ro = _risk_at_t(sub, old_col, t)
        Ts = sub["time"].values.astype(float)
        Es = sub["event"].values.astype(bool)
        case = Es & (Ts <= t)
        ctrl = Ts > t
        if case.sum() == 0 or ctrl.sum() == 0:
            return np.nan, np.nan
        up = rn > ro
        dn = rn < ro
        nri = (up[case].mean() - dn[case].mean()) - (up[ctrl].mean() - dn[ctrl].mean())
        idi = (rn[case].mean() - ro[case].mean()) - (rn[ctrl].mean() - ro[ctrl].mean())
        return nri, idi

    nri_app, idi_app = comp(df)
    rng = np.random.RandomState(seed)
    nris, idis = [], []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        try:
            v_n, v_i = comp(df.iloc[idx].reset_index(drop=True))
            if not (np.isnan(v_n) or np.isnan(v_i)):
                nris.append(v_n)
                idis.append(v_i)
        except Exception:
            continue
    def pct(a):
        return (np.nan, np.nan) if not a else tuple(np.nanpercentile(a, [2.5, 97.5]))
    return (nri_app, *pct(nris), idi_app, *pct(idis))


# ----------------------------------------------------------------------------
# time-dependent cumulative/dynamic AUC (IPCW)
# ----------------------------------------------------------------------------
def censoring_km(time, event):
    km = KaplanMeierFitter()
    km.fit(time, event_observed=(1 - event))
    return km


def dyn_auc(stage, time, event, t, km_cens):
    S = np.asarray(stage, float)
    T = np.asarray(time, float)
    E = np.asarray(event, bool)
    case_i = np.nonzero(E & (T <= t))[0]
    ctrl_i = np.nonzero(T > t)[0]
    if len(case_i) == 0 or len(ctrl_i) == 0:
        return np.nan
    g = km_cens.survival_function_at_times(pd.Series(T[case_i])).values.ravel()
    w = 1.0 / np.clip(g, 1e-3, 1.0)
    denom = w.sum() * len(ctrl_i)
    if denom == 0:
        return np.nan
    Sc = S[ctrl_i]
    num = 0.0
    for k, i in enumerate(case_i):
        Si = S[i]
        num += w[k] * (np.count_nonzero(Si > Sc) + 0.5 * np.count_nonzero(Si == Sc))
    return num / denom


# ----------------------------------------------------------------------------
# log-rank + stage별 사건율
# ----------------------------------------------------------------------------
def logrank_stats(df, stage_col):
    lr = multivariate_logrank_test(df["time"].values, df[stage_col].values,
                                   df["event"].values)
    return float(lr.test_statistic), float(lr.p_value)


def stage_event_rates(df, stage_col):
    d = df[["time", "event", stage_col]].dropna()
    rows = []
    for g, sub in d.groupby(stage_col):
        rows.append({"stage": g, "n": len(sub), "events": int(sub["event"].sum()),
                     "event_rate": sub["event"].mean(),
                     "median_fu": float(sub["time"].median())})
    r = pd.DataFrame(rows).sort_values("stage")
    if len(r) > 2:
        corr, p = spearmanr(r["stage"], r["event_rate"])
        r["monotone_spearman"] = corr
        r["monotone_p"] = p
    return r


# ----------------------------------------------------------------------------
# Part A: univariate stage 비교 콤보
# ----------------------------------------------------------------------------
def analyze_combo(args):
    label, (new_col, old_col, pair_name), seed = args
    t0 = time.time()
    df = load_data(label)
    pair_df = build_pair_df(df, new_col, old_col)
    n = len(pair_df)
    nb_c, nb_o, nb_n = (N_BOOT_CINDEX[seed], N_BOOT_OPTIM[seed], N_BOOT_NRI[seed])

    row = {"label": label, "pair": pair_name, "seed": seed, "n": n}
    row.update(bootstrap_combined(pair_df, nb_c, seed))
    row["opt_C_new"] = optimism_correct(pair_df, "stage_new", nb_o, seed)
    row["opt_C_old"] = optimism_correct(pair_df, "stage_old", nb_o, seed)

    dfc = pair_df.rename(columns={"stage_new": "stage"})
    dfc_old = pair_df.rename(columns={"stage_old": "stage"})
    for tag, dfx in [("new", dfc), ("old", dfc_old)]:
        try:
            cph = CoxPHFitter(penalizer=0.0).fit(dfx, duration_col="time",
                                                 event_col="event")
            ci = cph.confidence_intervals_.loc["stage"]
            row[f"HR_{tag}"] = float(np.exp(cph.params_["stage"]))
            row[f"HR_{tag}_lo"] = float(np.exp(ci["95% lower-bound"]))
            row[f"HR_{tag}_hi"] = float(np.exp(ci["95% upper-bound"]))
            row[f"p_{tag}"] = float(cph.summary.loc["stage", "p"])
        except Exception:
            row[f"HR_{tag}"] = row[f"HR_{tag}_lo"] = row[f"HR_{tag}_hi"] = np.nan
            row[f"p_{tag}"] = np.nan

    nri, nri_lo, nri_hi, idi, idi_lo, idi_hi = \
        nri_idi(pair_df, "stage_new", "stage_old", NRI_TIME, nb_n, seed)
    row.update({"NRI": nri, "NRI_lo": nri_lo, "NRI_hi": nri_hi,
                "IDI": idi, "IDI_lo": idi_lo, "IDI_hi": idi_hi,
                "NRI_excl_0": bool((nri_lo > 0) or (nri_hi < 0)) if not np.isnan(nri_lo) else False})

    chi_n, p_n = logrank_stats(pair_df, "stage_new")
    chi_o, p_o = logrank_stats(pair_df, "stage_old")
    row.update({"chi2_new": chi_n, "p_lr_new": p_n,
                "chi2_old": chi_o, "p_lr_old": p_o})

    km_c = censoring_km(pair_df["time"].values, pair_df["event"].values)
    for tt in AUC_TIMES:
        row[f"AUC_new_{tt}m"] = dyn_auc(pair_df["stage_new"].values,
                                        pair_df["time"].values,
                                        pair_df["event"].values, tt, km_c)
        row[f"AUC_old_{tt}m"] = dyn_auc(pair_df["stage_old"].values,
                                        pair_df["time"].values,
                                        pair_df["event"].values, tt, km_c)

    print(f"[A] {label} | {pair_name} | seed={seed} | n={n} "
          f"| C_new={row['C_new']:.3f} C_old={row['C_old']:.3f} "
          f"dC={row['dC']:+.3f}[{row['dC_lo']:.3f},{row['dC_hi']:.3f}] "
          f"| NRI={nri:.3f} | {time.time()-t0:.0f}s", flush=True)
    return row


# ----------------------------------------------------------------------------
# Part B: multivariate (X -> 각 y) — 전처리 CSV 사용
# ----------------------------------------------------------------------------
PREP_CSV = os.path.join(BASE_DIR, "preprocessed_data.csv")

# 다변량 X: 전처리 CSV 기준 (one-hot 포함, 수정병기 구성요소 제외)
#   - 수정병기 구성요소(PD_01vs2, ENE, LN meta count, largest node,
#     LN tumor size, contra_bilateral, N stage, bone invasion)는 mStage와 공선성
#     → 모델 A/B/C 공통 X에서 제외 (실험 1 docx 재현에서만 사용)
#   - CCRT_bin 제외: Tx_3tier_2와 완전 공선성 (VIF=∞, CCRT_bin == Tx_3tier_2)
#   - 희소 범주(subsite_5=2명, subsite_8=3명, HPV/P16_1=2명)는 build_multi_df에서
#     "기타"로 병합 → 리샘플 완전분리 방지
X_MULTI_BASE = ["age", "TIL", "differentiation", "tumor size (cm)", "DOI (mm)",
                "성별", "budding_01vs23", "PNI", "LVI", "RM", "TSR",
                "WPOI5_2tier",
                "subsite_2", "subsite_3", "subsite_4", "subsite_6", "subsite_7",
                "Tx_3tier_1", "Tx_3tier_2",
                "HPV/P16_2"]

SPARSE_ONEHOT = ["subsite_5", "subsite_8", "HPV/P16_1"]


def _merge_sparse(df):
    """희소 one-hot 범주를 '기타'로 병합 (리샘플 완전분리 방지)."""
    d = df.copy()
    for c in SPARSE_ONEHOT:
        if c in d.columns:
            d[c] = 0
    return d


def build_multi_df(label):
    d = pd.read_csv(PREP_CSV)
    d = _merge_sparse(d)
    y_time, y_event = f"{label}_time", f"{label}_event"
    keep = X_MULTI_BASE + [y_time, y_event, "ajcc8th_STAGE", "mStage"]
    m = d[keep].rename(columns={y_time: "time", y_event: "event"})
    return m.dropna().reset_index(drop=True)


def _fit_multi(dfx, covs):
    cph = CoxPHFitter(penalizer=0.1)
    try:
        cph.fit(dfx[["time", "event"] + list(covs)], duration_col="time",
                event_col="event")
        return cph
    except Exception:
        return None


def _risk_multi(dfx, covs, cph=None):
    if cph is None:
        cph = _fit_multi(dfx, covs)
    return cph.predict_partial_hazard(dfx[covs]).values.ravel()


def model_cindex_bootstrap(dfx, covs, n_boot, seed):
    """multivariate 모델 C-index: bootstrap CI(C_orig 분포) + optimism correction."""
    T = dfx["time"].values.astype(float)
    E = dfx["event"].values.astype(bool)
    n = len(dfx)
    rng = np.random.RandomState(seed)
    cph_app = _fit_multi(dfx, covs)
    if cph_app is None:
        return np.nan, np.nan, np.nan, np.nan
    C_app = concordance_index(T,
                              cph_app.predict_partial_hazard(dfx[covs]).values.ravel(),
                              E)
    c_orig_all = np.full(n_boot, np.nan)
    optim = np.full(n_boot, np.nan)
    for b in range(n_boot):
        idx = rng.randint(0, n, n)
        boot = dfx.iloc[idx]
        try:
            cph = _fit_multi(boot, covs)
            if cph is None:
                continue
            risk_orig = cph.predict_partial_hazard(dfx[covs]).values.ravel()
            c_orig_all[b] = concordance_index(T, risk_orig, E)
            risk_boot = cph.predict_partial_hazard(boot[covs]).values.ravel()
            C_boot = concordance_index(boot["time"].values, risk_boot,
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
    """반복 k-fold CV: 시드별 fold 분할 -> held-out C-index 평균."""
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
                cph = _fit_multi(dfx.iloc[train_idx], covs)
                if cph is None:
                    continue
                risk = cph.predict_partial_hazard(dfx.iloc[test_idx][covs]).values.ravel()
                scores.append(concordance_index(T[test_idx], risk, E[test_idx]))
            except Exception:
                continue
    if not scores:
        return np.nan, np.nan
    return float(np.mean(scores)), float(np.std(scores))


def repeated_cv_cindex_paired(dfx, covs_a, covs_b, k=5, repeats=10, seed=42):
    """S2(대체) 검증: 같은 fold에서 두 모델을 fit해 ΔCV = C_a - C_b 분포.

    같은 train/test 분할을 두 모델에 적용 → paired 비교가 가능.
    returns: (mean_a, sd_a, mean_b, sd_b, d_mean, d_lo, d_hi)
    """
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
                cph_a = _fit_multi(tr, covs_a)
                ca = concordance_index(
                    T[test_idx],
                    cph_a.predict_partial_hazard(te[covs_a]).values.ravel(),
                    E[test_idx]) if cph_a is not None else np.nan
            except Exception:
                ca = np.nan
            try:
                cph_b = _fit_multi(tr, covs_b)
                cb = concordance_index(
                    T[test_idx],
                    cph_b.predict_partial_hazard(te[covs_b]).values.ravel(),
                    E[test_idx]) if cph_b is not None else np.nan
            except Exception:
                cb = np.nan
            if not (np.isnan(ca) or np.isnan(cb)):
                sa.append(ca)
                sb.append(cb)
    if not sa:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
    sa, sb = np.array(sa), np.array(sb)
    d = sa - sb
    return (float(sa.mean()), float(sa.std()),
            float(sb.mean()), float(sb.std()),
            float(d.mean()), float(np.percentile(d, 2.5)),
            float(np.percentile(d, 97.5)))


def _nri_idi_multi(dfx, covs_a, covs_b, t=NRI_TIME):
    """중첩 모델 NRI/IDI. 같은 데이터에서 두 모델 예측 비교."""
    cph_a = _fit_multi(dfx, covs_a)
    cph_b = _fit_multi(dfx, covs_b)
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
    nb = N_BOOT_NRI[seed]
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
    nri_ba, idi_ba = _nri_idi_multi(m, covs_a, covs_b)
    nri_ca, idi_ca = _nri_idi_multi(m, covs_a, covs_c)
    nri_cb, idi_cb = _nri_idi_multi(m, covs_b, covs_c)
    for tag, nri, idi in [("B_vs_A", nri_ba, idi_ba),
                          ("C_vs_A", nri_ca, idi_ca),
                          ("C_vs_B", nri_cb, idi_cb)]:
        rows.append({"label": label, "seed": seed, "model": tag, "n": len(m),
                     "C_app": np.nan, "C_lo": np.nan, "C_hi": np.nan,
                     "C_optimism_corr": np.nan, "CV_mean": np.nan, "CV_sd": np.nan,
                     "NRI": nri, "IDI": idi})

    # S2(대체): 같은 CV fold에서 B(기존병기) vs C(수정병기) paired 비교
    cvB_m, cvB_sd, cvC_m, cvC_sd, d_cv, d_lo, d_hi = \
        repeated_cv_cindex_paired(m, covs_b, covs_c, k=5, repeats=10, seed=seed)
    rows.append({"label": label, "seed": seed, "model": "S2_C_vs_B", "n": len(m),
                 "CV_B_mean": cvB_m, "CV_B_sd": cvB_sd,
                 "CV_C_mean": cvC_m, "CV_C_sd": cvC_sd,
                 "CV_d_mean": d_cv, "CV_d_lo": d_lo, "CV_d_hi": d_hi,
                 "NRI": nri_cb, "IDI": idi_cb})

    # S3(증분 비교, 사용자 제안): NRI(C vs A) > NRI(B vs A), IDI 동일
    rows.append({"label": label, "seed": seed, "model": "S3_increment_C_vs_B",
                 "n": len(m),
                 "NRI_C_vs_A": nri_ca, "IDI_C_vs_A": idi_ca,
                 "NRI_B_vs_A": nri_ba, "IDI_B_vs_A": idi_ba,
                 "NRI_diff": (nri_ca - nri_ba) if not (np.isnan(nri_ca) or np.isnan(nri_ba)) else np.nan,
                 "IDI_diff": (idi_ca - idi_ba) if not (np.isnan(idi_ca) or np.isnan(idi_ba)) else np.nan})

    # HR 테이블 (full-data fit, seed 무관 — PRIMARY_SEED 행에만 기록)
    if seed == PRIMARY_SEED:
        cph_full = _fit_multi(m, covs)
        if cph_full is not None:
            summ = cph_full.summary
            for v in covs:
                ci = cph_full.confidence_intervals_.loc[v]
                rows.append({"label": label, "seed": seed, "model": f"HR_{v}", "n": len(m),
                             "HR": np.exp(cph_full.params_[v]),
                             "HR_lo": np.exp(ci["95% lower-bound"]),
                             "HR_hi": np.exp(ci["95% upper-bound"]),
                             "HR_p": summ.loc[v, "p"]})

    print(f"[B] {label} | seed={seed} | n={len(m)} | {time.time()-t0:.0f}s", flush=True)
    return rows


# ----------------------------------------------------------------------------
# KM 플롯 (결정론적)
# ----------------------------------------------------------------------------
def make_km_plots():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = ["AppleGothic", "sans-serif"]
    for label in LABELS:
        df = load_data(label)
        for new_col, old_col, pair_name in PAIRS:
            pair_df = build_pair_df(df, new_col, old_col)
            fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
            for ax, stage_col, title in [(axes[0], "stage_new", "Modified (수정 병기)"),
                                         (axes[1], "stage_old", "Conventional (기존 병기)")]:
                km = KaplanMeierFitter()
                for g in sorted(pair_df[stage_col].dropna().unique()):
                    sub = pair_df[pair_df[stage_col] == g]
                    km.fit(sub["time"], event_observed=sub["event"], label=f"Stage {int(g)}")
                    km.plot_survival_function(ax=ax, ci_show=True)
                ax.set_title(f"{title} - {label}")
                ax.set_xlabel("Months (수술일 기준)")
                ax.set_ylabel("Survival probability")
                ax.grid(alpha=0.3)
            fig.tight_layout()
            fname = os.path.join(OUT_DIR,
                                 f"KM_{label}_{pair_name.replace(' ','_').replace('/','_')}.png")
            fig.savefig(fname, dpi=130)
            plt.close(fig)
    print("[km] 플롯 저장 완료", flush=True)


# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------
def main():
    start = time.time()
    try:
        multiprocessing.set_start_method("fork", force=True)
    except RuntimeError:
        pass
    tasks_a = [(l, p, s) for l in LABELS for p in PAIRS for s in SEEDS]
    tasks_b = [(l, s) for l in LABELS for s in SEEDS]
    nproc = min(8, multiprocessing.cpu_count())
    print(f"Part A {len(tasks_a)}건 + Part B {len(tasks_b)}건 "
          f"(Pool {nproc}, fork)", flush=True)

    with multiprocessing.Pool(processes=nproc) as pool:
        rows_a = pool.map(analyze_combo, tasks_a)
        rows_b_flat = pool.map(analyze_multivariate, tasks_b)

    res_a = pd.DataFrame(rows_a)
    res_a.to_csv(os.path.join(OUT_DIR, "exp1_univariate.csv"), index=False,
                 encoding="utf-8-sig")

    rows_b = [r for chunk in rows_b_flat for r in chunk]
    res_b = pd.DataFrame(rows_b)
    res_b.to_csv(os.path.join(OUT_DIR, "exp1_multivariate.csv"), index=False,
                 encoding="utf-8-sig")

    rows_rate = []
    for label in LABELS:
        df = load_data(label)
        for new_col, old_col, pair_name in PAIRS:
            pair_df = build_pair_df(df, new_col, old_col)
            for tag, c in [("new", "stage_new"), ("old", "stage_old")]:
                r = stage_event_rates(pair_df, c)
                r.insert(0, "label", label)
                r.insert(1, "pair", pair_name)
                r.insert(2, "system", tag)
                rows_rate.append(r)
    rate_df = pd.concat(rows_rate, ignore_index=True)
    rate_df.to_csv(os.path.join(OUT_DIR, "exp1_event_rates.csv"), index=False,
                   encoding="utf-8-sig")
    if "monotone_spearman" in rate_df.columns:
        print("\n=== stage 단조성 (spearman) ===")
        print(rate_df[["label", "pair", "system", "monotone_spearman",
                       "monotone_p"]].drop_duplicates().to_string(index=False),
              flush=True)

    make_km_plots()
    print(f"\n[완료] {OUT_DIR} (총 {time.time()-start:.0f}s)", flush=True)


if __name__ == "__main__":
    main()