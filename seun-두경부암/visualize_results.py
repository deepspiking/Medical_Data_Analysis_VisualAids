# -*- coding: utf-8 -*-
"""
결과 시각화: 수정병기 구분력 검증 논문급 차트 (gowun 브랜치 스타일 참고, dpi=300)
------------------------------------------------------------------
실험 1-5 결과 CSV를 읽어 다음 차트 생성 (모두 results/exp1/plots/):

  V1. KM 곡선 (mStage vs ajcc8th, 라벨별)        : 생존곡선 분리 + log-rank p
  V2. C-index forest plot (쌍×y)                  : modified vs conventional + 95% CI
  V3. time-dependent AUC 곡선 (t=12/24/36/60)     : 라벨별 2개 시스템 비교
  V4. NRI/IDI 막대 차트 (쌍×y)                    : 재분류 개선 + 95% CI
  V5. score↔생존 산점도 (censored 구분)           : 실험 5
  V6. RMST 4분위 차트 (score별)                    : 실험 5, t*=36
  V7. 비단조성 진단 (quintile event rate 라인)     : 실험 5f
  V8. multivariate HR forest plot (모델 C)         : 실험 3
  V9. S2 paired CV 비교 (B vs C, 라벨별)           : 실험 3
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["AppleGothic", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE_DIR, "results", "exp1")
OUT = os.path.join(RES, "plots")
os.makedirs(OUT, exist_ok=True)

DPI = 300
Y_LABELS = ["PFS", "DSS", "LRRFS"]
PAIRS = ["mStage vs ajcc8th", "mTstage vs T stage", "mNstage vs N stage"]
COL_MOD = "#d62728"   # modified (빨강)
COL_OLD = "#1f77b4"   # conventional (파랑)


def _load(name):
    p = os.path.join(RES, name)
    return pd.read_csv(p) if os.path.exists(p) else None


def _save(fig, fname):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] {fname}")


# ----------------------------------------------------------------------------
# V1. KM 곡선 — 실험 2 (raw 데이터에서 직접 계산)
# ----------------------------------------------------------------------------
def v1_km():
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import multivariate_logrank_test
    xl = pd.ExcelFile(os.path.join(BASE_DIR, "정박사님께 드릴 raw data.xlsx"))
    raw = xl.parse("Sheet2")
    for y in Y_LABELS:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
        for ax, stage_col, title, color in [
                (axes[0], "mStage", "Modified staging", COL_MOD),
                (axes[1], "ajcc8th_STAGE", "Conventional staging", COL_OLD)]:
            d = raw[[stage_col, f"{y}_month", y]].dropna()
            d.columns = ["stage", "time", "event"]
            km = KaplanMeierFitter()
            for g in sorted(d["stage"].dropna().unique()):
                sub = d[d["stage"] == g]
                if len(sub) < 2:
                    continue
                km.fit(sub["time"], event_observed=sub["event"],
                       label=f"Stage {int(g)}")
                km.plot_survival_function(ax=ax, ci_show=True,
                                          color=color, linewidth=2)
            lr = multivariate_logrank_test(d["time"], d["stage"], d["event"])
            ax.text(0.03, 0.03, f"log-rank p = {lr.p_value:.2e}",
                    transform=ax.transAxes, fontsize=11,
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="gray"))
            ax.set_title(f"{title} - {y}")
            ax.set_xlabel("Months (수술일 기준)")
            ax.set_ylabel("Survival probability")
            ax.grid(alpha=0.3)
        _save(fig, f"V1_KM_{y}.png")


# ----------------------------------------------------------------------------
# V2. C-index forest plot — 실험 2 (exp1_univariate.csv)
# ----------------------------------------------------------------------------
def v2_cindex_forest(df):
    if df is None:
        return
    d = df[df["seed"] == 42].copy()
    fig, ax = plt.subplots(figsize=(10, 7))
    ypos = 0
    ticks = []
    for pair in PAIRS:
        for y in Y_LABELS:
            r = d[(d["pair"] == pair) & (d["label"] == y)]
            if r.empty:
                continue
            r = r.iloc[0]
            ax.errorbar([r["C_old"]], [ypos], xerr=[[r["C_old"] - r["C_old_lo"]],
                                                    [r["C_old_hi"] - r["C_old"]]],
                        fmt="o", color=COL_OLD, capsize=3, label="Conventional" if ypos == 0 else None)
            ax.errorbar([r["C_new"]], [ypos + 0.35], xerr=[[r["C_new"] - r["C_new_lo"]],
                                                           [r["C_new_hi"] - r["C_new"]]],
                        fmt="s", color=COL_MOD, capsize=3, label="Modified" if ypos == 0 else None)
            ticks.append((ypos + 0.175, f"{pair}\n{y}"))
            ypos += 1
    ax.set_yticks([t for t, _ in ticks])
    ax.set_yticklabels([l for _, l in ticks], fontsize=9)
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.6)
    ax.set_xlim(0.3, 1.0)
    ax.set_xlabel("C-index (95% CI)")
    ax.set_title("C-index: Modified vs Conventional staging (seed=42)")
    ax.legend(loc="lower right")
    ax.grid(axis="x", alpha=0.3)
    _save(fig, "V2_Cindex_forest.png")


# ----------------------------------------------------------------------------
# V3. time-dependent AUC 곡선 — 실험 2
# ----------------------------------------------------------------------------
def v3_td_auc(df):
    if df is None:
        return
    d = df[df["seed"] == 42]
    ts = [12, 24, 36, 60]
    for y in Y_LABELS:
        fig, ax = plt.subplots(figsize=(8, 5.5))
        for pair in PAIRS:
            r = d[(d["pair"] == pair) & (d["label"] == y)]
            if r.empty:
                continue
            r = r.iloc[0]
            auc_new = [r[f"AUC_new_{t}m"] for t in ts]
            auc_old = [r[f"AUC_old_{t}m"] for t in ts]
            ax.plot(ts, auc_new, "o-", color=COL_MOD, linewidth=2,
                    label=f"Modified ({pair})")
            ax.plot(ts, auc_old, "s--", color=COL_OLD, linewidth=1.5,
                    label=f"Conventional ({pair})")
        ax.axhline(0.5, color="gray", linestyle=":", alpha=0.6)
        ax.set_xlabel("Time (months)")
        ax.set_ylabel("Cumulative/dynamic AUC (IPCW)")
        ax.set_title(f"Time-dependent AUC - {y}")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        _save(fig, f"V3_AUC_{y}.png")


# ----------------------------------------------------------------------------
# V4. NRI/IDI 막대 — 실험 2
# ----------------------------------------------------------------------------
def v4_nri_idi(df):
    if df is None:
        return
    d = df[df["seed"] == 42].copy()
    d = d[["pair", "label", "NRI", "NRI_lo", "NRI_hi", "IDI", "IDI_lo", "IDI_hi"]].dropna(subset=["NRI"])
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, (val, lo, hi), title, color in [
            (axes[0], ("NRI", "NRI_lo", "NRI_hi"), "NRI (24개월 재분류)", COL_MOD),
            (axes[1], ("IDI", "IDI_lo", "IDI_hi"), "IDI", COL_OLD)]:
        labels = [f"{p}\n{y}" for p, y in zip(d["pair"], d["label"])]
        x = np.arange(len(d))
        ax.bar(x, d[val], color=color, alpha=0.8)
        ax.errorbar(x, d[val], yerr=[d[val] - d[lo], d[hi] - d[val]],
                    fmt="none", ecolor="black", capsize=3)
        ax.axhline(0, color="gray", linestyle="--", alpha=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8, rotation=20, ha="right")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.3)
    _save(fig, "V4_NRI_IDI.png")


# ----------------------------------------------------------------------------
# V5. score↔생존 산점도 — 실험 5
# ----------------------------------------------------------------------------
def v5_scatter():
    xl = pd.ExcelFile(os.path.join(BASE_DIR, "정박사님께 드릴 raw data.xlsx"))
    raw = xl.parse("Sheet2")
    for y in Y_LABELS:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)
        for ax, stage_col, title, color in [
                (axes[0], "mStage", "Modified", COL_MOD),
                (axes[1], "ajcc8th_STAGE", "Conventional", COL_OLD)]:
            d = raw[[stage_col, f"{y}_month", y]].dropna()
            d.columns = ["score", "time", "event"]
            jit = np.random.RandomState(42).normal(0, 0.08, len(d))
            for ev, marker, c, lab in [(1, "o", "#d62728", "Event"),
                                       (0, "^", "#1f77b4", "Censored")]:
                m = d["event"] == ev
                ax.scatter(d.loc[m, "score"] + jit[m], d.loc[m, "time"],
                           marker=marker, color=c, alpha=0.75, s=40, label=lab)
            ax.set_xlabel(f"Stage ({stage_col})")
            ax.set_ylabel(f"{y} months")
            ax.set_title(f"{title} - {y}")
            ax.legend(fontsize=9)
            ax.grid(alpha=0.3)
        _save(fig, f"V5_scatter_{y}.png")


# ----------------------------------------------------------------------------
# V6. RMST 4분위 — 실험 5 (exp1_rmst.csv)
# ----------------------------------------------------------------------------
def v6_rmst(df):
    if df is None:
        return
    fig, ax = plt.subplots(figsize=(11, 6))
    d = df.copy()
    labels, vals, los, his, cols = [], [], [], [], []
    for y in Y_LABELS:
        for sc in ["mStage", "ajcc8th_STAGE"]:
            r = d[(d["y"] == y) & (d["score_col"] == sc)]
            if r.empty:
                continue
            for _, row in r.iterrows():
                labels.append(f"{y} | {sc} Q{int(row['quartile'])+1}")
                vals.append(row["RMST"])
                los.append(row["RMST_lo"])
                his.append(row["RMST_hi"])
                cols.append(COL_MOD if sc == "mStage" else COL_OLD)
    x = np.arange(len(labels))
    ax.bar(x, vals, color=cols, alpha=0.85)
    ax.errorbar(x, vals, yerr=[np.array(vals) - np.array(los),
                               np.array(his) - np.array(vals)],
                fmt="none", ecolor="black", capsize=2.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7, rotation=45, ha="right")
    ax.set_ylabel(f"RMST (t* = 36개월)")
    ax.set_title("RMST by score quartile (red=Modified, blue=Conventional)")
    ax.grid(axis="y", alpha=0.3)
    _save(fig, "V6_RMST_quartile.png")


# ----------------------------------------------------------------------------
# V7. 비단조성 진단 — 실험 5f (exp1_nonmonotonicity.csv)
# ----------------------------------------------------------------------------
def v7_nonmonotone(df):
    if df is None:
        return
    d = df[df["y"] == "PFS"]
    top = d.head(9)
    fig, axes = plt.subplots(3, 3, figsize=(14, 11))
    for ax, (_, r) in zip(axes.flat, top.iterrows()):
        rates = eval(r["event_rates"]) if isinstance(r["event_rates"], str) else r["event_rates"]
        ax.plot(range(1, len(rates) + 1), rates, "o-", color=COL_MOD, linewidth=2)
        ax.set_title(f"{r['var']} [{r['classification']}]", fontsize=9)
        ax.set_xlabel("Quintile")
        ax.set_ylabel("PFS event rate")
        ax.grid(alpha=0.3)
        ax.set_ylim(0, 1)
    _save(fig, "V7_nonmonotonicity.png")


# ----------------------------------------------------------------------------
# V8. multivariate HR forest plot — 실험 3 (exp1_multivariate.csv HR_ 행)
# ----------------------------------------------------------------------------
def v8_hr_forest(df):
    if df is None:
        return
    d = df[df["model"].str.startswith("HR_", na=False)].copy()
    if d.empty:
        return
    d = d[d["seed"] == 42]
    for y in Y_LABELS:
        r = d[d["label"] == y].copy()
        if r.empty:
            continue
        r = r.sort_values("HR", ascending=True).reset_index(drop=True)
        fig, ax = plt.subplots(figsize=(9, max(6, 0.45 * len(r))))
        ypos = np.arange(len(r))
        ax.errorbar(r["HR"], ypos, xerr=[r["HR"] - r["HR_lo"], r["HR_hi"] - r["HR"]],
                    fmt="o", color=COL_MOD, capsize=3, linewidth=1.5)
        ax.axvline(1, color="gray", linestyle="--", alpha=0.7)
        ax.set_yticks(ypos)
        ax.set_yticklabels(r["model"].str.replace("HR_", ""), fontsize=9)
        ax.set_xscale("log")
        ax.set_xlabel("Hazard ratio (95% CI, log scale)")
        ax.set_title(f"Multivariate Cox HR - {y} (모델 A 기준)")
        ax.grid(axis="x", alpha=0.3)
        _save(fig, f"V8_HR_forest_{y}.png")


# ----------------------------------------------------------------------------
# V9. S2 paired CV 비교 — 실험 3 (exp1_multivariate.csv S2_C_vs_B 행)
# ----------------------------------------------------------------------------
def v9_s2_cv(df):
    if df is None:
        return
    d = df[(df["model"] == "S2_C_vs_B") & (df["seed"] == 42)].copy()
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(d))
    width = 0.35
    ax.bar(x - width / 2, d["CV_B_mean"], width, color=COL_OLD, alpha=0.85,
           label="Conventional (B: 공변량+ajcc8th)")
    ax.bar(x + width / 2, d["CV_C_mean"], width, color=COL_MOD, alpha=0.85,
           label="Modified (C: 공변량+mStage)")
    for i, (_, r) in enumerate(d.iterrows()):
        ax.text(i - width / 2, r["CV_B_mean"] + 0.01,
                f"{r['CV_B_mean']:.3f}", ha="center", fontsize=9)
        ax.text(i + width / 2, r["CV_C_mean"] + 0.01,
                f"{r['CV_C_mean']:.3f}", ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(d["label"])
    ax.set_ylabel("CV C-index (5-fold × 10회 반복)")
    ax.set_title("S2 대체 검증: 같은 CV fold에서 B vs C")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    _save(fig, "V9_S2_paired_CV.png")


def main():
    uni = _load("exp1_univariate.csv")
    multi = _load("exp1_multivariate.csv")
    rmst = _load("exp1_rmst.csv")
    nonmono = _load("exp1_nonmonotonicity.csv")

    v1_km()
    v2_cindex_forest(uni)
    v3_td_auc(uni)
    v4_nri_idi(uni)
    v5_scatter()
    v6_rmst(rmst)
    v7_nonmonotone(nonmono)
    v8_hr_forest(multi)
    v9_s2_cv(multi)
    print(f"\n모든 차트 저장: {OUT}")


if __name__ == "__main__":
    main()