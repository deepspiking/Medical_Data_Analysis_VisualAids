# -*- coding: utf-8 -*-
"""실험1 v2 forest plot: docx 참조 vs 재현(x0) HR — 방향·유의성 일치 색상."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(BASE, "results", "exp1", "exp1_docx_reproduction.csv"))
cmp = d[d["docx_HR"].notna() & d["HR"].notna()].copy()
cmp["ok"] = cmp["direction_match"] & cmp["sig_match"]
cmp["dir_only"] = cmp["direction_match"] & ~cmp["sig_match"]
cmp = cmp.sort_values(["outcome", "type", "var"]).reset_index(drop=True)

plt.rcParams["font.family"] = ["AppleGothic", "sans-serif"]
fig, ax = plt.subplots(figsize=(9, 13))
labels, colors = [], []
y = np.arange(len(cmp))
for i, (_, r) in enumerate(cmp.iterrows()):
    color = "#2e8b57" if r["ok"] else ("#e67e22" if r["dir_only"] else "#d62728")
    ax.errorbar(np.log(r["docx_HR"]), y[i], xerr=[[np.log(r["docx_HR"]) - np.log(r["docx_lo"])],
                                                  [np.log(r["docx_hi"]) - np.log(r["docx_HR"])]],
                fmt="o", color="steelblue", ms=6, capsize=2, zorder=3)
    ax.errorbar(np.log(r["HR"]), y[i], xerr=[[np.log(r["HR"]) - np.log(r["lo"])],
                                             [np.log(r["hi"]) - np.log(r["HR"])]],
                fmt="s", color=color, ms=6, capsize=2, zorder=3)
    labels.append(f"{r['outcome']} {r['type'][:4]} {r['var']}")
    colors.append(color)
ax.axvline(0, color="grey", lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=8)
xt = [-2, -1, 0, 1, 2, 3, 4]
ax.set_xticks(xt)
ax.set_xticklabels([f"{np.exp(v):.1f}" for v in xt])
ax.set_xlabel("HR (log scale) — 파랑 원=docx, 사각형=재현('x'→0)")
from matplotlib.lines import Line2D
handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue", label="docx 참조"),
           Line2D([0], [0], marker="s", color="w", markerfacecolor="#2e8b57", label="완전 일치"),
           Line2D([0], [0], marker="s", color="w", markerfacecolor="#e67e22", label="방향만 일치"),
           Line2D([0], [0], marker="s", color="w", markerfacecolor="#d62728", label="불일치")]
ax.legend(handles=handles, loc="lower right", fontsize=8)
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
out = os.path.join(BASE, "results", "exp1", "실험1_rawdata_vs_bootstrap_forest.png")
fig.savefig(out, dpi=130)
print("저장:", out)
