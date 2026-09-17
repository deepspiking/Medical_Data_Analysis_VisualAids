# -*- coding: utf-8 -*-
"""
제안 N grouping(mN′ 4군) KM 곡선 + AJCC N 비교 (n=133)
출력: results/exp1/ngroup_proposal/KM_{y}_mNprime_vs_AJCCN.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "results", "exp1", "ngroup_proposal")
os.makedirs(OUT, exist_ok=True)

d = pd.read_csv(os.path.join(BASE_DIR, "preprocessed_data_full.csv"))
AJCC = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 6: 5}
cnt = d["LN meta count_val"].fillna(0)
dep = d["LN tumor size (mm)_val"].fillna(0)
ene = d["ENE_val"].fillna(0)
contra = d["contra_bilateral_val"].fillna(0)

mnp = np.zeros(len(d), dtype=int)
npos = cnt >= 1
high = (ene == 1) | (contra == 1) | (cnt >= 3)
mnp[npos & ~high & (dep <= 10)] = 1
mnp[npos & ~high & (dep > 10)] = 2
mnp[npos & high] = 3

plt.rcParams["font.family"] = ["AppleGothic", "sans-serif"]
LABELS = {0: "mN'0 (림프절 없음)", 1: "mN'1 (1–2, deposit≤10)",
          2: "mN'2 (1–2, deposit>10)", 3: "mN'3 (≥3 / ENE+ / 양측)"}
AJCC_LABELS = {0: "N0", 1: "N1", 2: "N2a", 3: "N2b", 4: "N2c", 5: "N3b"}

for y in ["PFS", "DSS", "OS", "LRRFS"]:
    T = d[f"{y}_time"].values.astype(float)
    E = d[f"{y}_event"].values.astype(int)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4))
    for ax, g, lab_map, title in [
            (axes[0], mnp, LABELS, "제안 mN′ 4군"),
            (axes[1], d["N stage_val"].map(AJCC).values, AJCC_LABELS, "기존 AJCC N")]:
        km = KaplanMeierFitter()
        for gg in sorted(pd.Series(g).dropna().unique()):
            m = g == gg
            if m.sum() == 0:
                continue
            km.fit(T[m], event_observed=E[m], label=f"{lab_map.get(int(gg), gg)} (n={int(m.sum())})")
            km.plot_survival_function(ax=ax, ci_show=False)
        ax.set_title(f"{title} — {y}")
        ax.set_xlabel("Months")
        ax.set_ylabel("Survival probability")
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"KM_{y}_mNprime_vs_AJCCN.png"), dpi=130)
    plt.close(fig)
    print("saved", f"KM_{y}_mNprime_vs_AJCCN.png")

# reclassification
aj = pd.Series(d["N stage_val"].map(AJCC).values)
print("\nmN′ vs AJCC N crosstab:")
print(pd.crosstab(aj, pd.Series(mnp), rownames=["AJCC N"], colnames=["mN′"]).to_string())
