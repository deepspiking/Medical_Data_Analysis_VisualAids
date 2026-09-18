# -*- coding: utf-8 -*-
"""
수정병기 기반 예후 Nomogram + Calibration + DCA + C-index (정세운 선생님 요청)
================================================================================
의료계 표준 예측도구 세트를 파이썬으로 구현:
  1) Nomogram            : 다변량 Cox → 개별 환자 3년/5년 생존확률 계산 도표
  2) Calibration plot    : 예측 생존확률 vs 관측(KM) — 5분위 + bootstrap CI
  3) DCA (survival)      : IPCW 보정 net benefit vs threshold probability
  4) C-index             : Harrell C + bootstrap 95% CI (+ optimism-corrected)

입력 : preprocessed_data_full.csv (n=133, #52 정정 반영)
출력 : results/exp1/nomogram/  (그림 + 요약 txt/csv)

실행:
  python3 nomogram_analysis.py --endpoints OS,PFS --times 36,60
"""
import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from experiment_utils import fast_cindex

warnings.filterwarnings("ignore")

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "preprocessed_data_full.csv")
OUT = os.path.join(BASE, "results", "exp1", "nomogram")
os.makedirs(OUT, exist_ok=True)

plt.rcParams["font.family"] = ["AppleGothic", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

YMAP = {"OS": ("OS_time", "OS_event"), "PFS": ("PFS_time", "PFS_event"),
        "DSS": ("DSS_time", "DSS_event"), "LRRFS": ("LRRFS_time", "LRRFS_event")}

# 후보 예후인자 (병기 구성요소와 중복되는 원시 T/N 변수는 제외, 수정병기 축은 포함)
CANDIDATES = ["age", "male", "tumor size (cm)", "DOI (mm)", "differentiation",
              "budding_01vs23", "PNI", "LVI", "RM", "TIL", "TSR",
              "WPOI5_2tier", "HPV/P16_2", "CCRT_bin", "mTstage", "mNstage"]
FORCE_IN = ["mTstage", "mNstage"]   # 연구 핵심(수정병기 축) — 반드시 유지
LABELS = {
    "age": "Age (years)", "male": "Sex (male)", "tumor size (cm)": "Tumor size (cm)",
    "DOI (mm)": "DOI (mm)", "differentiation": "Differentiation", "budding_01vs23": "Budding",
    "PNI": "PNI", "LVI": "LVI", "RM": "Resection margin", "TIL": "TIL", "TSR": "TSR",
    "WPOI5_2tier": "WPOI5", "HPV/P16_2": "HPV/P16 (+)", "CCRT_bin": "CCRT",
    "mTstage": "Modified T stage", "mNstage": "Modified N stage",
}


def load():
    d = pd.read_csv(DATA)
    d = d.copy()
    d["male"] = (pd.to_numeric(d["성별"], errors="coerce") == 1).astype(int)
    return d


def prep_xy(d, y):
    tc, ec = YMAP[y]
    cols = ["연구번호"] + [c for c in CANDIDATES if c in d.columns] + [tc, ec]
    m = d[cols].rename(columns={tc: "T", ec: "E"}).dropna().reset_index(drop=True)
    m["E"] = m["E"].astype(int)
    return m


# ---------------------------------------------------------------------------
# 변수 선택 (univariate p<0.20 → multivariable backward elimination p<0.10)
# ---------------------------------------------------------------------------
def univariate_p(df, feats):
    out = {}
    for f in feats:
        try:
            cph = CoxPHFitter(penalizer=0.0)
            cph.fit(df[["T", "E", f]], duration_col="T", event_col="E")
            out[f] = float(cph.summary.loc[f, "p"])
        except Exception:
            out[f] = np.nan
    return out


def select_vars(df, y):
    base = [c for c in CANDIDATES if c in df.columns]
    uni = univariate_p(df, base)
    keep = [f for f in base if (f in FORCE_IN) or (uni.get(f, 1.0) < 0.20)]
    # backward elimination by p-value (>0.10 제거, 단 FORCE_IN 유지)
    while len(keep) > 0:
        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(df[["T", "E"] + keep], duration_col="T", event_col="E")
        p = cph.summary["p"].drop(labels=[f for f in FORCE_IN if f in keep], errors="ignore")
        if len(p) == 0:
            break
        worst = p.idxmax()
        if p.max() > 0.10:
            keep.remove(worst)
        else:
            break
    return keep, uni


def fit_final(df, feats):
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(df[["T", "E"] + feats], duration_col="T", event_col="E")
    return cph


# ---------------------------------------------------------------------------
# 예측 헬퍼
# ---------------------------------------------------------------------------
def surv_at(cph, df, t, feats):
    lp = cph.predict_log_partial_hazard(df[feats]).values.ravel()
    s0 = cph.baseline_survival_
    idx = s0.index[s0.index <= t]
    if len(idx) == 0:
        return np.full(len(df), np.nan)
    s0t = float(s0.iloc[len(idx) - 1, 0])
    return s0t ** np.exp(lp)


def check_baseline_mapping(cph, df, feats, t):
    """formula S0(t)^exp(lp) == predict_survival_function 인지 self-check."""
    a = surv_at(cph, df, t, feats)
    b = cph.predict_survival_function(df[feats], times=[t]).values.ravel()
    return float(np.nanmax(np.abs(a - b)))


# ---------------------------------------------------------------------------
# Nomogram 계산 + 그리기
# ---------------------------------------------------------------------------
def nomogram_points(cph, df, feats):
    beta = cph.params_
    # 각 변수의 선형예측 기여 범위
    contribs = {}
    for f in feats:
        vals = df[f].values.astype(float)
        b = float(beta[f])
        c = b * vals
        contribs[f] = (c.min(), c.max())
    ranges = {f: (contribs[f][1] - contribs[f][0]) for f in feats}
    R = max(max(ranges.values()), 1e-9)
    max_points = {f: 100.0 * ranges[f] / R for f in feats}

    def points_of(f, v):
        b = float(beta[f])
        return 100.0 * (b * v - contribs[f][0]) / R

    lp_shift = float(sum(contribs[f][0] for f in feats))
    return dict(contribs=contribs, max_points=max_points, R=R,
                lp_shift=lp_shift, points_of=points_of)


def nice_ticks(lo, hi, target=6):
    span = hi - lo
    if span <= 0:
        return np.array([lo])
    raw = span / max(target - 1, 1)
    mag = 10 ** np.floor(np.log10(raw))
    step = mag
    for m in (1, 2, 2.5, 5, 10):
        step = m * mag
        if span / step <= target + 1:
            break
    start = np.ceil(lo / step) * step
    return np.arange(start, hi + 1e-9, step)


def draw_nomogram(cph, df, feats, times, endpoint, fname, cindex_txt):
    nc = nomogram_points(cph, df, feats)
    beta = cph.params_
    total_max = float(sum(nc["max_points"].values()))
    # lin pred ↔ total points 관계: LP = lp_shift + (R/100)*TP  (raw βx 기준)
    # lifelines predict_log_partial_hazard = Σβx - Σβ*mean 이므로 상수만 다름 → 실측으로 보정
    lp_meas = cph.predict_log_partial_hazard(df[feats]).values.ravel()
    tp_meas = np.array([sum(nc["points_of"](f, r[f]) for f in feats)
                        for _, r in df[feats].iterrows()])
    A = np.polyfit(tp_meas, lp_meas, 1)
    slope, intercept = A[0], A[1]

    def lp_of_tp(tp):
        return intercept + slope * tp

    def tp_of_surv(t, s):
        s0 = cph.baseline_survival_
        idx = s0.index[s0.index <= t]
        if len(idx) == 0:
            return np.nan
        s0t = float(s0.iloc[len(idx) - 1, 0])
        if not (0 < s < 1) or not (0 < s0t < 1):
            return np.nan
        lp = np.log(np.log(s) / np.log(s0t))
        return (lp - intercept) / slope

    nrows = 4 + len(feats) + len(times)
    fig, ax = plt.subplots(figsize=(11, 0.62 * nrows + 1.8))
    ypos = {}
    yy = 0.0
    ypos["Points"] = yy
    yy += 1.35
    for f in feats:
        ypos[f] = yy
        yy += 1.0
    ypos["Total"] = yy + 0.25
    ypos["Risk"] = ypos["Total"] + 1.05
    ypos["LP"] = ypos["Risk"] + 1.05
    for t in times:
        ypos[f"t{t}"] = ypos["LP"] + 1.0 + 0.85 * times.index(t)
    ymax = max(ypos.values()) + 0.6
    var_palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b",
                   "#17becf", "#e377c2", "#7f7f7f", "#bcbd22", "#aec7e8"]
    var_colors = {f: var_palette[i % len(var_palette)] for i, f in enumerate(feats)}
    for f in ("mTstage", "mNstage"):
        if f in var_colors:
            var_colors[f] = "#c0392b"

    XMAX = 100.0  # 물리 좌표 0..100

    def top_x(p):      # points/variable 행
        return p

    def bot_x(tp):     # total points 이하 행
        return tp / total_max * XMAX if total_max else 0.0

    # Points ruler
    ax.plot([top_x(0), top_x(100)], [ypos["Points"]] * 2, color="black", lw=1.2)
    for p in range(0, 101, 2):
        big = (p % 10 == 0)
        ax.plot([top_x(p), top_x(p)], [ypos["Points"], ypos["Points"] + (0.16 if big else 0.09)],
                color="black", lw=1.0 if big else 0.6)
        if big:
            ax.text(top_x(p), ypos["Points"] + 0.22, str(p), ha="center", va="bottom", fontsize=8)
    ax.text(-2.5, ypos["Points"], "Points", ha="right", va="center", fontsize=10, fontweight="bold")

    for f in feats:
        yr = ypos[f]
        mp = nc["max_points"][f]
        c = var_colors[f]
        ax.plot([top_x(0), top_x(mp)], [yr] * 2, color=c, lw=2.4, solid_capstyle="round")
        vals = sorted(df[f].dropna().unique())
        if len(vals) > 8:
            ticks = nice_ticks(float(min(vals)), float(max(vals)))
        else:
            ticks = vals
        for v in ticks:
            p = nc["points_of"](f, v)
            ax.plot([top_x(p), top_x(p)], [yr, yr - 0.15], color=c, lw=1.1)
            lab = f"{v:g}" if isinstance(v, (int, float, np.floating)) else str(v)
            ax.text(top_x(p), yr - 0.26, lab, ha="center", va="top", fontsize=8.5)
        ax.text(-2.5, yr, LABELS.get(f, f), ha="right", va="center",
                fontsize=10.5, fontweight="bold", color=c)

    # Total Points
    yT = ypos["Total"]
    ax.plot([bot_x(0), bot_x(total_max)], [yT] * 2, color="black", lw=1.4)
    step = 50 if total_max > 300 else (25 if total_max > 150 else 20)
    for tp in np.arange(0, total_max + 1e-9, step):
        ax.plot([bot_x(tp), bot_x(tp)], [yT, yT - 0.16], color="black", lw=1.0)
        ax.text(bot_x(tp), yT - 0.26, f"{tp:g}", ha="center", va="top", fontsize=8)
    ax.text(-2.5, yT, "Total Points", ha="right", va="center", fontsize=10, fontweight="bold")

    from matplotlib.patches import Rectangle
    yR = ypos["Risk"]
    tp33, tp67 = np.percentile(tp_meas, [33.3, 66.7])
    bands = [(0.0, tp33, "#2ca02c", "Low-risk"),
             (tp33, tp67, "#ff7f0e", "Medium-risk"),
             (tp67, total_max, "#d62728", "High-risk")]
    for x0, x1, col, lab in bands:
        if x1 <= x0:
            continue
        xa, xb = bot_x(x0), bot_x(x1)
        ax.add_patch(Rectangle((xa, yR - 0.34), xb - xa, 0.68,
                               facecolor=col, edgecolor="white", lw=1.2, zorder=2))
        ax.text((xa + xb) / 2, yR, lab, ha="center", va="center",
                fontsize=10, fontweight="bold", color="white", zorder=3)
    ax.text(-2.5, yR, "Risk group", ha="right", va="center", fontsize=10, fontweight="bold")

    # Linear Predictor
    yL = ypos["LP"]
    lp_lo, lp_hi = lp_of_tp(0), lp_of_tp(total_max)
    ax.plot([bot_x(0), bot_x(total_max)], [yL] * 2, color="black", lw=1.2)
    for lpv in np.linspace(np.floor(lp_lo * 2) / 2, np.ceil(lp_hi * 2) / 2, 9):
        tp = (lpv - intercept) / slope
        if 0 <= tp <= total_max:
            ax.plot([bot_x(tp), bot_x(tp)], [yL, yL - 0.14], color="black", lw=0.9)
            ax.text(bot_x(tp), yL - 0.24, f"{lpv:.1f}", ha="center", va="top", fontsize=8)
    ax.text(-2.5, yL, "Linear Predictor", ha="right", va="center", fontsize=10, fontweight="bold")

    # Survival scales
    s_colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
    for ti, t in enumerate(times):
        yS = ypos[f"t{t}"]
        col = s_colors[ti % len(s_colors)]
        label = f"{t}-month survival" if t % 12 else f"{t//12}-year survival"
        ax.text(-2.5, yS, label, ha="right", va="center", fontsize=10.5,
                fontweight="bold", color=col)
        s_ticks = [0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
        tps = [(s, tp_of_surv(t, s)) for s in s_ticks]
        tps = [(s, tp) for s, tp in tps if tp is not None and np.isfinite(tp) and 0 <= tp <= total_max]
        if not tps:
            continue
        ax.plot([bot_x(tps[0][1]), bot_x(tps[-1][1])], [yS] * 2, color=col, lw=2.2)
        for s, tp in tps:
            ax.plot([bot_x(tp), bot_x(tp)], [yS, yS - 0.15], color=col, lw=1.0)
            ax.text(bot_x(tp), yS - 0.26, f"{s:g}", ha="center", va="top", fontsize=8.5)

    ax.set_xlim(-16, 104)
    ax.set_ylim(ymax, -0.7)
    ax.axis("off")
    ax.set_title(f"Prognostic Nomogram for {endpoint} (n=133)\n{cindex_txt}",
                 fontsize=12.5, fontweight="bold", pad=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {fname}")
    return dict(total_max=total_max, intercept=intercept, slope=slope,
                max_points=nc["max_points"])


# ---------------------------------------------------------------------------
# Calibration (5분위, KM observed @ t)
# ---------------------------------------------------------------------------
def calibration(cph, df, feats, y, t, fname, n_groups=5):
    S = surv_at(cph, df, t, feats)
    risk = 1 - S
    order = np.argsort(risk)
    groups = np.array_split(order, n_groups)
    pts = []
    for g in groups:
        pred = float(np.mean(risk[g]))
        sub = df.iloc[g]
        km = KaplanMeierFitter()
        km.fit(sub["T"], event_observed=sub["E"])
        obs = 1 - float(km.predict(t))
        # Greenwood 기반 CI
        ci = km.confidence_interval_survival_function_
        row = ci.index[ci.index <= t]
        if len(row):
            lo = 1 - float(ci.loc[row[-1], ci.columns[1]])
            hi = 1 - float(ci.loc[row[-1], ci.columns[0]])
        else:
            lo = hi = np.nan
        pts.append((pred, obs, lo, hi, len(g)))
    dfc = pd.DataFrame(pts, columns=["pred", "obs", "lo", "hi", "n"])

    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    ax.plot([0, 1], [0, 1], ls="--", color="gray", lw=1, label="Ideal")
    ax.errorbar(dfc["pred"], dfc["obs"],
                yerr=[dfc["obs"] - dfc["lo"], dfc["hi"] - dfc["obs"]],
                fmt="o-", color="#c0392b", ms=6, capsize=3, lw=1.5, label="Nomogram")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel(f"Predicted {t//12}-year event probability", fontsize=11)
    ax.set_ylabel(f"Observed {t//12}-year event (KM)", fontsize=11)
    ax.set_title(f"Calibration — {y}, {t//12}-year", fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {fname}")
    return dfc


# ---------------------------------------------------------------------------
# DCA (survival, IPCW)
# ---------------------------------------------------------------------------
def dca_survival(risk, T, E, t, thresholds):
    T = np.asarray(T, float); E = np.asarray(E, bool)
    n = len(T)
    km_c = KaplanMeierFitter()
    km_c.fit(T, event_observed=(1 - E).astype(int))
    G_at = lambda tt: float(np.atleast_1d(km_c.predict(np.atleast_1d(tt)))[0])
    G_T = km_c.predict(T).values
    G_T = np.clip(G_T, 1e-6, 1.0)
    Gt = max(G_at(t), 1e-6)

    event_by_t = (T <= t) & E
    nonevent_by_t = T > t

    def nb(treat):
        tp = np.sum(event_by_t[treat] / G_T[treat]) / n
        fp = np.sum(nonevent_by_t[treat] / Gt) / n
        return tp, fp

    out = {"threshold": thresholds}
    # treat-all
    tp_all, fp_all = nb(np.ones(n, bool))
    out["all"] = np.array([tp_all - fp_all * (pt / (1 - pt)) for pt in thresholds])
    out["none"] = np.zeros_like(thresholds)
    for name, r in risk.items():
        vals = []
        for pt in thresholds:
            treat = r >= pt
            tp, fp = nb(treat)
            vals.append(tp - fp * (pt / (1 - pt)))
        out[name] = np.array(vals)
    return out


def draw_dca(dca, y, t, fname):
    th = dca["threshold"]
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    colors = {"nomogram": "#c0392b", "mStage": "#1f77b4", "mTstage": "#2ca02c",
              "all": "gray", "none": "black"}
    for k in ["nomogram", "mStage", "mTstage", "all", "none"]:
        if k not in dca:
            continue
        ls = "--" if k in ("all", "none") else "-"
        ax.plot(th, dca[k], ls=ls, color=colors.get(k, None), lw=1.8,
                label={"nomogram": "Nomogram", "mStage": "mStage only",
                       "mTstage": "mTstage only", "all": "Treat all",
                       "none": "Treat none"}.get(k, k))
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlim(min(th), max(th)); ax.set_ylim(-0.05, max(0.4, np.max(dca["all"]) * 1.05))
    ax.set_xlabel(f"Threshold probability ({t//12}-year)", fontsize=11)
    ax.set_ylabel("Net benefit", fontsize=11)
    ax.set_title(f"Decision Curve Analysis — {y}, {t//12}-year", fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {fname}")


# ---------------------------------------------------------------------------
# C-index (bootstrap CI + optimism-corrected)
# ---------------------------------------------------------------------------
def cindex_stats(cph, df, feats, y, n_boot=1000, seed=42):
    T = df["T"].values.astype(float); E = df["E"].values.astype(bool)
    lp = cph.predict_log_partial_hazard(df[feats]).values.ravel()
    c_app = fast_cindex(lp, T, E)
    rng = np.random.RandomState(seed)
    cs = []
    for _ in range(n_boot):
        i = rng.randint(0, len(T), len(T))
        cs.append(fast_cindex(lp[i], T[i], E[i]))
    lo, hi = np.nanpercentile(cs, [2.5, 97.5])
    return dict(C=c_app, lo=lo, hi=hi)


def cindex_validate(df, feats, n_boot=300, seed=42):
    """Harrell C apparent + bootstrap optimism-corrected (bootstrap 내부에서 모델 재적합)."""
    T = df["T"].values.astype(float); E = df["E"].values.astype(bool)
    cph0 = CoxPHFitter(penalizer=0.1)
    cph0.fit(df[["T", "E"] + feats], duration_col="T", event_col="E")
    lp0 = cph0.predict_log_partial_hazard(df[feats]).values.ravel()
    c_app = fast_cindex(lp0, T, E)
    rng = np.random.RandomState(seed)
    n = len(df)
    opts, c_res = [], []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        c_res.append(fast_cindex(lp0[idx], T[idx], E[idx]))
        boot = df.iloc[idx]
        try:
            cphb = CoxPHFitter(penalizer=0.1)
            cphb.fit(boot[["T", "E"] + feats], duration_col="T", event_col="E")
        except Exception:
            continue
        lp_b = cphb.predict_log_partial_hazard(boot[feats]).values.ravel()
        lp_o = cphb.predict_log_partial_hazard(df[feats]).values.ravel()
        c_b = fast_cindex(lp_b, boot["T"].values, boot["E"].values.astype(bool))
        c_o = fast_cindex(lp_o, T, E)
        opts.append(c_b - c_o)
    optimism = float(np.nanmean(opts)) if opts else np.nan
    lo, hi = np.nanpercentile(c_res, [2.5, 97.5]) if c_res else (np.nan, np.nan)
    return dict(C_app=c_app, optimism=optimism, C_corr=c_app - optimism, lo=lo, hi=hi)


def calibration_slope_ipcw(pred_risk, T, E, t):
    """D'Agostino-Nam IPCW-weighted logistic calibration at t: slope (ideal 1), intercept (ideal 0)."""
    from scipy.special import logit
    from sklearn.linear_model import LogisticRegression
    T = np.asarray(T, float); E = np.asarray(E, bool); pred_risk = np.asarray(pred_risk, float)
    km_c = KaplanMeierFitter(); km_c.fit(T, event_observed=(1 - E).astype(int))
    G_T = np.clip(km_c.predict(T).values, 1e-6, 1.0)
    Gt = max(float(np.atleast_1d(km_c.predict(np.atleast_1d(t)))[0]), 1e-6)
    case = (T <= t) & E
    ctrl = T > t
    w = np.zeros(len(T)); y = case.astype(int)
    w[case] = 1.0 / G_T[case]
    w[ctrl] = 1.0 / Gt
    mask = w > 0
    x = logit(np.clip(pred_risk, 1e-4, 1 - 1e-4)).reshape(-1, 1)
    try:
        lr = LogisticRegression(penalty="l2", C=1e6, solver="lbfgs", max_iter=5000)
        lr.fit(x[mask], y[mask], sample_weight=w[mask])
        return float(lr.coef_[0][0]), float(lr.intercept_[0])
    except Exception:
        return np.nan, np.nan


def td_auc(risk, T, E, t):
    T = np.asarray(T, float); E = np.asarray(E, bool); risk = np.asarray(risk, float)
    km_c = KaplanMeierFitter(); km_c.fit(T, event_observed=(1 - E).astype(int))
    G_T = np.clip(km_c.predict(T).values, 1e-6, 1.0)
    Gt = max(float(np.atleast_1d(km_c.predict(np.atleast_1d(t)))[0]), 1e-6)
    case = (T <= t) & E
    ctrl = T > t
    if case.sum() == 0 or ctrl.sum() == 0:
        return np.nan
    wc = 1.0 / G_T[case]
    wl = np.full(ctrl.sum(), 1.0 / Gt)
    Rc = risk[case][:, None]
    Rl = risk[ctrl][None, :]
    conc = np.sum(wc[:, None] * wl[None, :] * ((Rc > Rl) + 0.5 * (Rc == Rl)))
    return float(conc / (wc.sum() * wl.sum()))


def td_roc(risk, T, E, t):
    T = np.asarray(T, float); E = np.asarray(E, bool); risk = np.asarray(risk, float)
    km_c = KaplanMeierFitter(); km_c.fit(T, event_observed=(1 - E).astype(int))
    G_T = np.clip(km_c.predict(T).values, 1e-6, 1.0)
    Gt = max(float(np.atleast_1d(km_c.predict(np.atleast_1d(t)))[0]), 1e-6)
    case = (T <= t) & E
    ctrl = T > t
    wc = 1.0 / G_T[case]
    wl = np.full(ctrl.sum(), 1.0 / Gt)
    cuts = np.unique(np.quantile(risk, np.linspace(0, 1, 101)))
    tpr, fpr = [], []
    for c in cuts:
        tp = wc[risk[case] >= c].sum()
        fp = wl[risk[ctrl] >= c].sum()
        tpr.append(tp / wc.sum())
        fpr.append(fp / wl.sum())
    order = np.argsort(fpr)
    return np.array(fpr)[order], np.array(tpr)[order], td_auc(risk, T, E, t)


def risk_group_km(cph, df, feats, endpoint, times, fname):
    lp = cph.predict_log_partial_hazard(df[feats]).values.ravel()
    q = np.percentile(lp, [33.3, 66.7])
    grp = np.where(lp <= q[0], 0, np.where(lp <= q[1], 1, 2))
    names = ["Low risk", "Intermediate", "High risk"]
    colors = ["#2ca02c", "#1f77b4", "#d62728"]
    from lifelines.statistics import multivariate_logrank_test
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    for g in [0, 1, 2]:
        m = grp == g
        if m.sum() == 0:
            continue
        km = KaplanMeierFitter()
        km.fit(df["T"].values[m], event_observed=df["E"].values[m],
               label=f"{names[g]} (n={int(m.sum())})")
        km.plot_survival_function(ax=ax, ci_show=False, color=colors[g])
    p = multivariate_logrank_test(df["T"].values, grp, df["E"].values).p_value
    ptxt = "p<0.001" if p < 0.001 else f"p={p:.3f}"
    ax.set_title(f"{endpoint} by nomogram risk group (log-rank {ptxt})",
                 fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Months"); ax.set_ylabel("Survival probability")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {fname}")
    tvals = {int(t): td_auc(lp, df["T"].values, df["E"].values, t) for t in times}
    return tvals


def draw_td_roc(cph, df, feats, endpoint, times, fname):
    lp = cph.predict_log_partial_hazard(df[feats]).values.ravel()
    T = df["T"].values; E = df["E"].values
    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    ax.plot([0, 1], [0, 1], ls="--", color="gray", lw=1)
    for t, color in zip(times, ["#1f77b4", "#c0392b"]):
        fpr, tpr, auc = td_roc(lp, T, E, t)
        ax.plot(fpr, tpr, color=color, lw=1.9, label=f"{t//12}-year (AUC={auc:.3f})")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("1 - Specificity", fontsize=11)
    ax.set_ylabel("Sensitivity", fontsize=11)
    ax.set_title(f"Time-dependent ROC — {endpoint}", fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {fname}")


def run_endpoint(d, y, times):
    df = prep_xy(d, y)
    feats, uni = select_vars(df, y)
    cph = fit_final(df, feats)
    print(f"\n===== {y} (n={len(df)}, events={int(df['E'].sum())}) =====")
    print("selected:", feats)
    print(cph.summary[["coef", "exp(coef)", "p"]].round(4).to_string())

    # self-check
    for t in times:
        err = check_baseline_mapping(cph, df, feats, t)
        print(f"  [check] baseline mapping max|Δ| @t={t}: {err:.2e}")

    ci = cindex_validate(df, feats, n_boot=300)
    print(f"  C apparent={ci['C_app']:.3f} | optimism={ci['optimism']:.3f} | "
          f"corrected={ci['C_corr']:.3f} | refit CI [{ci['lo']:.3f}, {ci['hi']:.3f}]")
    c_txt = (f"Harrell C = {ci['C_app']:.3f} (optimism-corrected {ci['C_corr']:.3f}; "
             f"95% CI {ci['lo']:.3f}–{ci['hi']:.3f})")
    draw_nomogram(cph, df, feats, times, y, f"nomogram_{y}.png", c_txt)

    risk = {}
    cal_rows = []
    lp = cph.predict_log_partial_hazard(df[feats]).values.ravel()
    for t in times:
        risk[t] = 1 - surv_at(cph, df, t, feats)
        dfc = calibration(cph, df, feats, y, t, f"calibration_{y}_{t//12}yr.png")
        slope, intercept = calibration_slope_ipcw(risk[t], df["T"].values, df["E"].values, t)
        cal_rows.append({"time": t, "calibration_slope": round(slope, 3),
                         "calibration_intercept": round(intercept, 3)})
        print(f"  calibration slope/intercept @{t}mo = {slope:.3f} / {intercept:.3f}")
        dfc.to_csv(os.path.join(OUT, f"calibration_{y}_{t//12}yr.csv"), index=False)
    pd.DataFrame(cal_rows).to_csv(os.path.join(OUT, f"calibration_slope_{y}.csv"), index=False)

    tvals = risk_group_km(cph, df, feats, y, times, f"risk_km_{y}.png")
    draw_td_roc(cph, df, feats, y, times, f"roc_{y}.png")
    print("  time-dependent AUC:", {k: round(v, 3) for k, v in tvals.items()})

    th = np.arange(0.05, 0.81, 0.025)
    dca = dca_survival({"nomogram": risk[times[0]]}, df["T"].values, df["E"].values,
                       times[0], th)
    # mStage / mTstage 단독 모델(단변량 Cox) 기반 위험도로 DCA 비교
    for f in ["mStage", "mTstage"]:
        if f in d.columns:
            m = d.set_index("연구번호").reindex(df["연구번호"]).reset_index()
            dfk = pd.DataFrame({"T": df["T"].values, "E": df["E"].values,
                                f: m[f].values.astype(float)})
            try:
                cphk = CoxPHFitter(penalizer=0.1)
                cphk.fit(dfk, duration_col="T", event_col="E")
                rk = 1 - surv_at(cphk, dfk, times[0], [f])
                dca[f] = np.array([_nb_at(rk, dfk["T"].values, dfk["E"].values,
                                          times[0], pt) for pt in th])
            except Exception as e:
                print(f"  [dca] {f} skip: {e}")
    draw_dca(dca, y, times[0], f"dca_{y}_{times[0]//12}yr.png")

    # summary
    summ = cph.summary[["coef", "exp(coef)", "p"]].copy()
    summ.to_csv(os.path.join(OUT, f"model_{y}.csv"), encoding="utf-8-sig")
    with open(os.path.join(OUT, f"summary_{y}.txt"), "w", encoding="utf-8") as f:
        f.write(f"[{y}] n={len(df)}, events={int(df['E'].sum())}\n")
        f.write(f"selected vars: {feats}\n")
        f.write(c_txt + "\n\n")
        f.write(summ.round(4).to_string() + "\n")
    return ci


def _nb_at(risk, T, E, t, pt):
    """단일 threshold net benefit (IPCW)."""
    T = np.asarray(T, float); E = np.asarray(E, bool); risk = np.asarray(risk, float)
    n = len(T)
    km_c = KaplanMeierFitter(); km_c.fit(T, event_observed=(1 - E).astype(int))
    G_T = np.clip(km_c.predict(T).values, 1e-6, 1.0)
    Gt = max(float(np.atleast_1d(km_c.predict(np.atleast_1d(t)))[0]), 1e-6)
    treat = risk >= pt
    tp = np.sum(((T <= t) & E)[treat] / G_T[treat]) / n
    fp = np.sum((T > t)[treat] / Gt) / n
    return tp - fp * (pt / (1 - pt))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--endpoints", default="OS,PFS")
    ap.add_argument("--times", default="36,60")
    args = ap.parse_args()
    times = [int(x) for x in args.times.split(",")]
    d = load()
    for y in [x.strip() for x in args.endpoints.split(",") if x.strip()]:
        run_endpoint(d, y, times)
    print("\n[출력]", OUT)


if __name__ == "__main__":
    main()
