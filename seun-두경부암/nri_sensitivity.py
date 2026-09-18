# -*- coding: utf-8 -*-
"""mNstage LRRFS NRI 신호 민감도 분석 (고속판)."""
import os
import sys
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "results", "exp1", "nri_sensitivity")
os.makedirs(OUT, exist_ok=True)
DATA = os.path.join(BASE, "preprocessed_data_full.csv")
Y = ["PFS", "DSS", "OS", "LRRFS"]
PAIRS = [("mStage", "ajcc8th_STAGE", "mStage vs ajcc8th"),
         ("mTstage", "T stage", "mTstage vs T"),
         ("mNstage", "N stage_val", "mNstage vs N")]
SEED = 42


def load():
    d = pd.read_csv(DATA)
    for c in ["T stage", "ajcc8th_STAGE", "mTstage", "mNstage", "mStage", "N stage_val"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna(subset=["mStage", "mTstage", "mNstage", "N stage_val"]).reset_index(drop=True)


def fit_risk(d, col, times):
    df = d[["time", "event", col]].rename(columns={col: "s"})
    cph = CoxPHFitter(penalizer=0.0)
    cph.fit(df, duration_col="time", event_col="event")
    sf = cph.predict_survival_function(df[["s"]], times=list(times))
    return 1.0 - sf.values


def nri_parts(rn, ro, T, E, t):
    case = E & (T <= t)
    ctrl = T > t
    if case.sum() == 0 or ctrl.sum() == 0:
        return np.nan, np.nan, np.nan
    up = rn > ro
    dn = rn < ro
    ev = up[case].mean() - dn[case].mean()
    ne = up[ctrl].mean() - dn[ctrl].mean()
    return ev - ne, ev, ne


def boot_nri_fixed(rn, ro, T, E, t, b, seed):
    rng = np.random.RandomState(seed)
    n = len(T)
    vals = []
    for _ in range(b):
        i = rng.randint(0, n, n)
        v = nri_parts(rn[i], ro[i], T[i], E[i], t)[0]
        if np.isfinite(v):
            vals.append(v)
    vals = np.array(vals)
    lo, hi = np.percentile(vals, [2.5, 97.5])
    p = 2.0 * min((vals <= 0).mean(), (vals >= 0).mean())
    return lo, hi, min(p, 1.0)


def main():
    d = load()
    rows = []
    for new, old, name in PAIRS:
        for y in Y:
            dd = pd.DataFrame({"time": d[f"{y}_time"], "event": d[f"{y}_event"].astype(int),
                               new: d[new], old: d[old]}).dropna().reset_index(drop=True)
            T = dd["time"].values.astype(float)
            E = dd["event"].values.astype(bool)
            rn = fit_risk(dd, new, [24.0])[0]
            ro = fit_risk(dd, old, [24.0])[0]
            nri, ev, ne = nri_parts(rn, ro, T, E, 24.0)
            lo, hi, p = boot_nri_fixed(rn, ro, T, E, 24.0, 3000, SEED)
            rows.append({"pair": name, "y": y, "NRI": round(nri, 3), "lo": round(lo, 3),
                         "hi": round(hi, 3), "p": round(p, 4),
                         "NRI_event": round(ev, 3), "NRI_nonevent": round(ne, 3)})
    res = pd.DataFrame(rows)
    res["sig_raw"] = (res["lo"] > 0) | (res["hi"] < 0)
    order = res["p"].sort_values().index
    m = len(res)
    padj = np.empty(m)
    prev = 0.0
    for rank, idx in enumerate(order):
        val = min(1.0, res.loc[idx, "p"] * (m - rank))
        prev = max(prev, val)
        padj[idx] = prev
    res["p_holm"] = padj.round(4)
    res["sig_holm"] = res["p_holm"] < 0.05
    res.to_csv(os.path.join(OUT, "nri_all.csv"), index=False, encoding="utf-8-sig")
    print("=== NRI (t=24) 12개 비교 + Holm 보정 (fixed-score bootstrap B=3000) ===")
    print(res.to_string(index=False))

    dd = pd.DataFrame({"time": d["LRRFS_time"], "event": d["LRRFS_event"].astype(int),
                       "mNstage": d["mNstage"], "N": d["N stage_val"]}).dropna().reset_index(drop=True)
    T = dd["time"].values.astype(float)
    E = dd["event"].values.astype(bool)
    hrows = []
    for t in [12, 24, 36, 48]:
        rn = fit_risk(dd, "mNstage", [float(t)])[0]
        ro = fit_risk(dd, "N", [float(t)])[0]
        nri, ev, ne = nri_parts(rn, ro, T, E, float(t))
        lo, hi, p = boot_nri_fixed(rn, ro, T, E, float(t), 3000, SEED)
        hrows.append({"t_months": t, "NRI": round(nri, 3), "lo": round(lo, 3),
                      "hi": round(hi, 3), "p": round(p, 4),
                      "NRI_event": round(ev, 3), "NRI_nonevent": round(ne, 3)})
    hdf = pd.DataFrame(hrows)
    hdf.to_csv(os.path.join(OUT, "lrrfs_mNstage_horizons.csv"), index=False, encoding="utf-8-sig")
    print("\n=== mNstage LRRFS: 시간대별 NRI ===")
    print(hdf.to_string(index=False))

    rn24 = fit_risk(dd, "mNstage", [24.0])[0]
    ro24 = fit_risk(dd, "N", [24.0])[0]
    obs = nri_parts(rn24, ro24, T, E, 24.0)[0]
    rng = np.random.RandomState(SEED)
    null = np.array([nri_parts(rn24[rng.permutation(len(T))], ro24, T, E, 24.0)[0]
                     for _ in range(5000)])
    p_perm = float((np.abs(null) >= abs(obs)).mean())
    print("\n=== permutation test (mNstage 위험점수 재배정, B=5000) ===")
    print(f"observed NRI={obs:.3f}, permutation p={p_perm:.4f}")

    print("\n=== refit-bootstrap 재확인 (mNstage LRRFS t=24, B=200) ===")
    rng = np.random.RandomState(SEED)
    refit = []
    for _ in range(200):
        i = rng.randint(0, len(dd), len(dd))
        sub = dd.iloc[i].reset_index(drop=True)
        try:
            a = fit_risk(sub, "mNstage", [24.0])[0]
            b = fit_risk(sub, "N", [24.0])[0]
            v = nri_parts(a, b, T[i], E[i], 24.0)[0]
            if np.isfinite(v):
                refit.append(v)
        except Exception:
            pass
    refit = np.array(refit)
    print(f"refit NRI mean={refit.mean():.3f}, 95% CI "
          f"[{np.percentile(refit,2.5):.3f}, {np.percentile(refit,97.5):.3f}], "
          f"p={2*min((refit<=0).mean(),(refit>=0).mean()):.4f}")

    with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as f:
        f.write("NRI t=24 + Holm\n" + res.to_string(index=False) + "\n\n")
        f.write("mNstage LRRFS by horizon\n" + hdf.to_string(index=False) + "\n\n")
        f.write(f"permutation p={p_perm:.4f}\n")
        f.write(f"refit NRI mean={refit.mean():.3f} CI "
                f"[{np.percentile(refit,2.5):.3f},{np.percentile(refit,97.5):.3f}]\n")
    print("\n[출력]", OUT)


if __name__ == "__main__":
    main()
