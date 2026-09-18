# -*- coding: utf-8 -*-
"""
Stage를 제외한 구성요소 중요도 분석
- T/N stage(mTstage·mNstage·T stage·N stage) 제외
- 구성요소(size, DOI, PD, bone, depthbone, LN count, deposit, ENE, contra) + 임상인자
- 중요도: |β|×SD (표준화), |β|×range (nomogram 축 길이), HR
출력: results/exp1/component_importance/
"""
import os
import sys
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "results", "exp1", "component_importance")
os.makedirs(OUT, exist_ok=True)
plt.rcParams["font.family"] = ["AppleGothic", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

COMP = {
    "Tumor size (cm)": "tumor size (cm)",
    "DOI (mm)": "DOI (mm)",
    "PD cellularity": "PD_01vs2",
    "Bone invasion": "bone invasion_val",
    "Depth of bone inv. (mm)": "depth of bone invasion (mm)_val",
    "LN count": "LN meta count_val",
    "LN deposit (mm)": "LN tumor size (mm)_val",
    "ENE": "ENE_val",
    "Bilateral/contralat.": "contra_bilateral_val",
}
CLIN = {
    "Age": "age", "Male": "sex_male", "Differentiation": "differentiation",
    "PNI": "PNI", "LVI": "LVI", "WPOI5": "WPOI5_2tier",
    "CCRT": "CCRT_bin", "HPV/P16 (+)": "HPV_2",
}
Y = ["DSS", "OS", "PFS", "LRRFS"]


def load():
    d = pd.read_csv(os.path.join(BASE, "preprocessed_data_full.csv"))
    d["sex_male"] = (pd.to_numeric(d["성별"], errors="coerce") == 1).astype(int)
    d["HPV_2"] = pd.to_numeric(d["HPV/P16_2"], errors="coerce").fillna(0)
    for c in list(COMP.values()) + ["differentiation", "PNI", "LVI", "WPOI5_2tier", "CCRT_bin", "age"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    for c in COMP.values():
        d[c] = d[c].fillna(0)
    return d


def fit_importance(d, y, var_map):
    cols = list(var_map.values())
    dd = d[cols + [f"{y}_time", f"{y}_event"]].rename(
        columns={f"{y}_time": "T", f"{y}_event": "E"}).dropna().reset_index(drop=True)
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(dd[["T", "E"] + cols], duration_col="T", event_col="E")
    sd = {c: dd[c].std(ddof=1) for c in cols}
    rng = {c: dd[c].max() - dd[c].min() for c in cols}
    rows = []
    for name, c in var_map.items():
        b = float(cph.params_[c])
        rows.append({
            "variable": name, "coef": round(b, 4),
            "HR": round(float(np.exp(b)), 3),
            "p": round(float(cph.summary.loc[c, "p"]), 4),
            "std_beta": round(abs(b) * sd[c], 4),
            "range_beta": round(abs(b) * rng[c], 4),
        })
    out = pd.DataFrame(rows)
    out["rank_std"] = out["std_beta"].rank(ascending=False).astype(int)
    return out.sort_values("std_beta", ascending=False).reset_index(drop=True)


def main():
    d = load()
    for y in Y:
        a = fit_importance(d, y, COMP)
        a.to_csv(os.path.join(OUT, f"importance_components_{y}.csv"),
                 index=False, encoding="utf-8-sig")
        print(f"\n===== {y} — 구성요소만 (stage 제외) =====")
        print(a[["variable", "HR", "p", "std_beta", "range_beta", "rank_std"]].to_string(index=False))

    vmap = {**COMP, **CLIN}
    for y in Y:
        b = fit_importance(d, y, vmap)
        b.to_csv(os.path.join(OUT, f"importance_all_{y}.csv"), index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax, y in zip(axes.ravel(), Y):
        a = pd.read_csv(os.path.join(OUT, f"importance_components_{y}.csv"))
        a = a.sort_values("std_beta")
        colors = ["#c0392b" if v in ("PD cellularity", "DOI (mm)", "Tumor size (cm)")
                  else "#7f8c8d" for v in a["variable"]]
        ax.barh(a["variable"], a["std_beta"], color=colors)
        for i, (v, s) in enumerate(zip(a["variable"], a["std_beta"])):
            ax.text(s + 0.01, i, f"{s:.2f}", va="center", fontsize=8)
        ax.set_title(f"{y} — 구성요소 중요도 (|β|×SD)", fontsize=11, fontweight="bold")
        ax.set_xlabel("standardized contribution")
        ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "component_importance.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("\n[save] component_importance.png")
    print("[출력]", OUT)


if __name__ == "__main__":
    main()
