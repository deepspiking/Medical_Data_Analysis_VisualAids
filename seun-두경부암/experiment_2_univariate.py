# -*- coding: utf-8 -*-
"""
실험 2: 수정병기 vs 기존병기 Univariate 구분력 (S1) — 단독 실행 스크립트
------------------------------------------------------------------
비교 쌍: mStage vs ajcc8th / mTstage vs T stage / mNstage vs N stage
Y: PFS/DSS/LRRFS (raw xlsx, 수술일 기준 월)
지표: C-index(bootstrap 10,000 CI, paired ΔC), AUC, KM+log-rank, NRI/IDI
"""
import os
import time
import multiprocessing
from functools import lru_cache
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from scipy.stats import spearmanr

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(BASE_DIR, "정박사님께 드릴 raw data.xlsx")
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")
os.makedirs(OUT_DIR, exist_ok=True)

LABELS = ["PFS", "DSS", "OS", "LRRFS"]
OUTCOME_SRC = {"PFS": ("PFS_month", "PFS"), "DSS": ("DSS_month", "DSS"),
               "OS": ("Death_month", "death"), "LRRFS": ("LRRFS_month", "LRRFS")}
PAIRS = [
    ("mStage",   "ajcc8th_STAGE", "mStage vs ajcc8th"),
    (" mTstage", "T stage",       "mTstage vs T stage"),
    ("mNstage",  "N stage",       "mNstage vs N stage"),
]
SEEDS = [42, 123, 2026, 777]
PRIMARY_SEED = 42

N_BOOT_CINDEX = {PRIMARY_SEED: 10_000, 123: 3_000, 2026: 3_000, 777: 3_000}
N_BOOT_OPTIM = {PRIMARY_SEED: 500, 123: 200, 2026: 200, 777: 200}
N_BOOT_NRI = {PRIMARY_SEED: 1_000, 123: 300, 2026: 300, 777: 300}
NRI_TIME = 24.0
AUC_TIMES = [12, 24, 36, 60]


@lru_cache(maxsize=None)
def load_data(label):
    xl = pd.ExcelFile(DATA_PATH)
    df = xl.parse("Sheet2")
    time_col, event_col = OUTCOME_SRC[label]
    keep = ["연구번호", "T stage", "N stage", "ajcc8th_STAGE",
            " mTstage", "mNstage", "mStage", "time", "event"]
    df = df.rename(columns={time_col: "time", event_col: "event"})
    df["event"] = df["event"].astype(int)
    df["time"] = df["time"].astype(float)
    df["N stage"] = pd.to_numeric(df["N stage"], errors="coerce")
    for c in ["T stage", "ajcc8th_STAGE", " mTstage", "mNstage", "mStage"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df[keep].copy()


def build_pair_df(df, new_col, old_col):
    out = pd.DataFrame({
        "stage_new": pd.to_numeric(df[new_col], errors="coerce"),
        "stage_old": pd.to_numeric(df[old_col], errors="coerce"),
        "time": df["time"].astype(float),
        "event": df["event"].astype(int),
    })
    return out.dropna(subset=["stage_new", "stage_old"]).reset_index(drop=True)


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
    def pct(a):
        return np.nanpercentile(a, [2.5, 97.5])
    return dict(
        C_new=c_new, C_new_lo=pct(cn)[0], C_new_hi=pct(cn)[1],
        C_old=c_old, C_old_lo=pct(co)[0], C_old_hi=pct(co)[1],
        dC=c_new - c_old, dC_lo=pct(d)[0], dC_hi=pct(d)[1],
        dC_excl_0=bool((pct(d)[0] > 0) or (pct(d)[1] < 0)),
    )


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


def logrank_stats(df, stage_col):
    lr = multivariate_logrank_test(df["time"].values, df[stage_col].values,
                                   df["event"].values)
    return float(lr.test_statistic), float(lr.p_value)


def analyze_combo(args):
    label, (new_col, old_col, pair_name), seed = args
    t0 = time.time()
    df = load_data(label)
    pair_df = build_pair_df(df, new_col, old_col)
    n = len(pair_df)
    nb_c = N_BOOT_CINDEX[seed]
    nb_o = N_BOOT_OPTIM[seed]
    nb_n = N_BOOT_NRI[seed]

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

    print(f"[실험2] {label} | {pair_name} | seed={seed} | n={n} "
          f"| C_new={row['C_new']:.3f} C_old={row['C_old']:.3f} "
          f"dC={row['dC']:+.3f}[{row['dC_lo']:.3f},{row['dC_hi']:.3f}] "
          f"| NRI={nri:.3f} | {time.time()-t0:.0f}s", flush=True)
    return row


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
    print("[실험2] KM 플롯 저장 완료", flush=True)


def main():
    import argparse
    start = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", default="PFS,DSS,OS,LRRFS",
                    help="콤마 구분 y label (기본 전체)")
    ap.add_argument("--seeds", default="42,123,2026,777", help="콤마 구분 시드")
    args = ap.parse_args()
    labels = [x for x in args.labels.split(",") if x]
    seeds = [int(x) for x in args.seeds.split(",") if x]
    try:
        multiprocessing.set_start_method("fork", force=True)
    except RuntimeError:
        pass
    tasks = [(l, p, s) for l in labels for p in PAIRS for s in seeds]
    nproc = min(8, multiprocessing.cpu_count())
    print(f"실험 2 (S1 univariate) | labels={labels} | {len(tasks)} 콤보 | Pool {nproc}, fork",
          flush=True)

    with multiprocessing.Pool(processes=nproc) as pool:
        rows = pool.map(analyze_combo, tasks)

    res = pd.DataFrame(rows)
    if labels == LABELS:
        res.to_csv(os.path.join(OUT_DIR, "exp1_univariate.csv"), index=False,
                   encoding="utf-8-sig")
    else:
        res.to_csv(os.path.join(OUT_DIR, "exp1_univariate_partial.csv"), index=False,
                   encoding="utf-8-sig")
    print(f"[완료] {time.time()-start:.0f}s — "
          f"{'results/exp1/exp1_univariate.csv' if labels == LABELS else 'exp1_univariate_partial.csv'}")
    if labels == LABELS:
        make_km_plots()


if __name__ == "__main__":
    main()