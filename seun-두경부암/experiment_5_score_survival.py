# -*- coding: utf-8 -*-
"""
실험 5: Score ↔ 실제 생존기간 상관 분석 (Model-agnostic)
------------------------------------------------------------------
사용자 아이디어: "생존 분석을 꼭 regression 모델로 할 필요는 없다.
모델 score 값과 실제 생존 기간을 직접 상관분석할 수 있다."

분석:
  5a stage-score ↔ 생존   : mStage/ajcc8th 값을 score로, uncensored만 Spearman + 산점도
  5b RMST                 : score 4분위 → RMST(t*=36) + bootstrap 95% CI
  5c TabICL score ↔ 생존  : TabICL regression score(CV) ↔ 실제 생존기간 (24개월 이진화 없이)
  5d score-quartile KM    : score 4분위 KM + log-rank
  5e 산점도               : score vs 생존월 (censored 별도 마커)

주의: 중도절단 편향 → uncensored 상관(보조) + C-index(주) + RMST 병행 보고
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

from experiment_utils import fast_cindex

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_CSV = os.path.join(BASE_DIR, "preprocessed_data.csv")
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")
os.makedirs(OUT_DIR, exist_ok=True)

SEED = 42
RMST_T = 36.0

Y_LABELS = ["PFS", "DSS", "OS", "LRRFS"]
PAIRS = [("mStage", "ajcc8th_STAGE", "mStage vs ajcc8th"),
         ("mTstage", "T stage", "mTstage vs T stage"),
         ("mNstage", "N stage_val", "mNstage vs N stage")]


def load_prep():
    return pd.read_csv(DATA_CSV)


def _rmst_boot_one(args):
    time, event, t, seed = args
    rng = np.random.RandomState(seed)
    n = len(time)
    i = rng.randint(0, n, n)
    k2 = KaplanMeierFitter()
    k2.fit(time[i], event_observed=event[i])
    s2 = k2.survival_function_
    i2 = s2.index <= t
    if i2.sum() == 0:
        return np.nan
    x2 = np.concatenate([[0], s2.index[i2]])
    y2 = np.concatenate([[1], s2["KM_estimate"].values[i2]])
    return float(np.trapz(y2, x2))


def rmst_estimate(time, event, t=RMST_T, n_boot=100, seed=SEED):
    """제한평균생존기간 RMST(t) = KM 곡선 0~t 적분 + bootstrap CI."""
    km = KaplanMeierFitter()
    km.fit(time, event_observed=event)
    sf = km.survival_function_
    idx = sf.index <= t
    if idx.sum() == 0:
        return np.nan, np.nan, np.nan
    x = np.concatenate([[0], sf.index[idx]])
    y = np.concatenate([[1], sf["KM_estimate"].values[idx]])
    rmst = float(np.trapz(y, x))

    rng = np.random.RandomState(seed)
    boots = []
    for k in range(n_boot):
        i = rng.randint(0, len(time), len(time))
        k2 = KaplanMeierFitter()
        k2.fit(time[i], event_observed=event[i])
        s2 = k2.survival_function_
        i2 = s2.index <= t
        if i2.sum() == 0:
            continue
        x2 = np.concatenate([[0], s2.index[i2]])
        y2 = np.concatenate([[1], s2["KM_estimate"].values[i2]])
        boots.append(float(np.trapz(y2, x2)))
    boots = np.array(boots)
    if len(boots) == 0:
        return rmst, np.nan, np.nan
    return rmst, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def rmst_diff_boot(t_hi, e_hi, t_lo, e_lo, t=RMST_T, n_boot=1000, seed=SEED):
    """두 그룹 RMST 차이(고위험−저위험)의 bootstrap CI와 양측 p (H0: diff=0)."""
    def rmst_arr(time, event):
        km = KaplanMeierFitter()
        km.fit(time, event_observed=event)
        sf = km.survival_function_
        idx = sf.index <= t
        if idx.sum() == 0:
            return np.nan
        x = np.concatenate([[0], sf.index[idx]])
        y = np.concatenate([[1], sf["KM_estimate"].values[idx]])
        return float(np.trapz(y, x))

    def rmst_pair(i_hi, i_lo):
        a = rmst_arr(time_hi[i_hi], event_hi[i_hi])
        b = rmst_arr(time_lo[i_lo], event_lo[i_lo])
        return a - b

    time_hi, event_hi = np.asarray(t_hi, float), np.asarray(e_hi, bool)
    time_lo, event_lo = np.asarray(t_lo, float), np.asarray(e_lo, bool)
    obs = rmst_pair(np.arange(len(time_hi)), np.arange(len(time_lo)))
    if np.isnan(obs):
        return obs, np.nan, np.nan, np.nan
    rng = np.random.RandomState(seed)
    diffs = []
    for _ in range(n_boot):
        i_hi = rng.randint(0, len(time_hi), len(time_hi))
        i_lo = rng.randint(0, len(time_lo), len(time_lo))
        v = rmst_pair(i_hi, i_lo)
        if not np.isnan(v):
            diffs.append(v)
    if len(diffs) < 200:
        return obs, np.nan, np.nan, np.nan
    lo, hi = np.percentile(diffs, 2.5), np.percentile(diffs, 97.5)
    diffs = np.array(diffs)
    p = 2.0 * min(float((diffs <= 0).mean()), float((diffs >= 0).mean()))
    p = min(p, 1.0)
    return obs, float(lo), float(hi), float(p)


def score_survival_analysis(d, score_col, y):
    """score ↔ 생존: C-index(주) + uncensored Spearman(보조) + RMST 4분위."""
    y_time, y_event = f"{y}_time", f"{y}_event"
    dfx = d[[score_col, y_time, y_event]].dropna()
    s = dfx[score_col].values.astype(float)
    t = dfx[y_time].values.astype(float)
    e = dfx[y_event].values.astype(bool)
    if len(np.unique(s)) < 2:
        return None

    c_idx = float(fast_cindex(s, t, e)) if s.max() != s.min() else np.nan

    unc = e  # uncensored = event 발생
    rho, p_rho = (np.nan, np.nan)
    if unc.sum() >= 5 and s[unc].std() > 0:
        rho, p_rho = spearmanr(s[unc], t[unc])

    # score 4분위 RMST (양수 보정: stage는 클수록 위험, score↑=위험↑)
    try:
        q = pd.qcut(s, 4, labels=False, duplicates="drop")
    except Exception:
        q = None
    rmst_rows = []
    if q is not None:
        for qi in sorted(pd.unique(q)):
            m = q == qi
            if m.sum() < 3:
                continue
            r, lo, hi = rmst_estimate(t[m], e[m])
            rmst_rows.append({"score_col": score_col, "y": y, "quartile": int(qi),
                              "n": int(m.sum()), "RMST": r, "RMST_lo": lo, "RMST_hi": hi})

    rmst_diff = None
    if q is not None:
        qmin, qmax = int(pd.unique(q).min()), int(pd.unique(q).max())
        m_hi, m_lo = q == qmax, q == qmin
        if m_hi.sum() >= 3 and m_lo.sum() >= 3:
            diff, dlo, dhi, dp = rmst_diff_boot(t[m_hi], e[m_hi], t[m_lo], e[m_lo])
            rmst_diff = {"score_col": score_col, "y": y,
                         "Q_high": qmax, "Q_low": qmin,
                         "n_high": int(m_hi.sum()), "n_low": int(m_lo.sum()),
                         "RMST_diff": diff, "RMST_diff_lo": dlo, "RMST_diff_hi": dhi,
                         "RMST_diff_p": dp}

    return {"score_col": score_col, "y": y, "n": len(dfx),
            "C_index": c_idx,
            "unc_n": int(unc.sum()), "spearman_rho": rho, "spearman_p": p_rho,
            "rmst_q": rmst_rows, "rmst_diff": rmst_diff}


def make_scatter(d, score_col, y):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = ["AppleGothic", "sans-serif"]
    y_time, y_event = f"{y}_time", f"{y}_event"
    dfx = d[[score_col, y_time, y_event]].dropna()
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for ev, marker, color, lab in [(True, "o", "#d62728", "Event"),
                                   (False, "^", "#1f77b4", "Censored")]:
        m = dfx[y_event] == ev
        ax.scatter(dfx.loc[m, score_col], dfx.loc[m, y_time],
                   marker=marker, color=color, alpha=0.75, label=lab, s=38)
    ax.set_xlabel(f"Score ({score_col})")
    ax.set_ylabel(f"{y} survival months")
    ax.set_title(f"{score_col} vs {y} (censored 별도 마커)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fname = os.path.join(OUT_DIR,
                         f"score_survival_{y}_{score_col.replace(' ','_')}.png")
    fig.savefig(fname, dpi=130)
    plt.close(fig)
    return fname


def nonmonotonicity_diagnostics(d, y, var_cols):
    """5f-① 개별 변수 quintile별 event rate → 단조/임계/U자 분류.

    각 변수를 quintile로 나눠 y 사건율을 계산하고:
      - Spearman(변수, 사건율)로 단조 방향/강도
      - 중간 quintile이 최고/최저보다 높으면 U자(비단조) 후보
    """
    y_time, y_event = f"{y}_time", f"{y}_event"
    rows = []
    for v in var_cols:
        dfx = d[[v, y_event]].dropna()
        if dfx[v].nunique() < 3 or dfx[y_event].sum() < 5:
            continue
        try:
            q = pd.qcut(dfx[v], 5, labels=False, duplicates="drop")
        except Exception:
            continue
        rates = []
        for qi in sorted(pd.unique(q)):
            m = q == qi
            rates.append(float(dfx.loc[m, y_event].mean()))
        rho, p_rho = spearmanr(
            [float(dfx.loc[q == qi, v].median()) for qi in sorted(pd.unique(q))],
            rates)
        is_monotone = p_rho < 0.05
        peak = int(np.argmax(rates))
        u_shape = (peak != 0) and (peak != len(rates) - 1) and is_monotone is False
        rows.append({"y": y, "var": v, "n_q": len(rates),
                     "event_rates": [round(r, 3) for r in rates],
                     "spearman_rho": round(float(rho), 3),
                     "spearman_p": round(float(p_rho), 4),
                     "peak_quintile": peak,
                     "classification": "U자(비단조)" if u_shape else
                                       ("단조" if is_monotone else "비단조(무패턴)")})
    return rows


def linear_score_cindex_compare(d, y, covs):
    """5f-② 선형 Cox score vs (플래그 포함) score C-index 비교.

    선형 모델(Cox)의 partial hazard를 score로 써서 C-index 계산.
    TabICL score는 experiment_4에서 별도 비교.
    """
    from lifelines import CoxPHFitter
    y_time, y_event = f"{y}_time", f"{y}_event"
    dfx = d[covs + [y_time, y_event]].dropna()
    if len(dfx) < 20 or dfx[y_event].sum() < 10:
        return None
    cph = CoxPHFitter(penalizer=0.0)
    cph.fit(dfx, duration_col=y_time, event_col=y_event)
    score = cph.predict_partial_hazard(dfx[covs]).values.ravel()
    c_idx = float(fast_cindex(score, dfx[y_time].values,
                              dfx[y_event].values.astype(bool)))
    return {"y": y, "n": len(dfx), "linear_Cox_C": c_idx,
            "note": "TabICL score C-index는 experiment_4에서 비교"}


def main():
    d = load_prep()
    rows = []
    rmst_all = []
    rmst_diff_all = []
    for y in Y_LABELS:
        for new_c, old_c, pair_name in PAIRS:
            for tag, score_col in [("modified", new_c), ("conventional", old_c)]:
                res = score_survival_analysis(d, score_col, y)
                if res is None:
                    continue
                rows.append({"pair": pair_name, "system": tag,
                             "y": y, "C_index": res["C_index"],
                             "unc_n": res["unc_n"],
                             "spearman_rho": res["spearman_rho"],
                             "spearman_p": res["spearman_p"]})
                for r in res["rmst_q"]:
                    rmst_all.append(r)
                if res["rmst_diff"] is not None:
                    rmst_diff_all.append(res["rmst_diff"])
                make_scatter(d, score_col, y)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT_DIR, "exp1_score_survival.csv"), index=False,
               encoding="utf-8-sig")
    pd.DataFrame(rmst_all).to_csv(os.path.join(OUT_DIR, "exp1_rmst.csv"),
                                  index=False, encoding="utf-8-sig")
    if rmst_diff_all:
        pd.DataFrame(rmst_diff_all).to_csv(os.path.join(OUT_DIR, "exp1_rmst_diff.csv"),
                                           index=False, encoding="utf-8-sig")

    # 5f-① 비단조성 진단 (개별 변수 quintile event rate)
    from scipy.stats import spearmanr
    var_cols = [c for c in d.columns
                if c not in ("연구번호",) and not c.endswith(("_event", "_time"))
                and not c.startswith(("subsite_", "Tx_3tier_", "HPV/P16_"))]
    diag_rows = []
    for y in Y_LABELS:
        diag_rows += nonmonotonicity_diagnostics(d, y, var_cols)
    if diag_rows:
        diag = pd.DataFrame(diag_rows)
        diag.to_csv(os.path.join(OUT_DIR, "exp1_nonmonotonicity.csv"), index=False,
                    encoding="utf-8-sig")
        print("\n=== 5f-① 개별 변수 관계 형태 (quintile event rate) ===")
        print(diag[["y", "var", "classification", "spearman_rho",
                    "spearman_p", "event_rates"]].to_string(index=False))

    # 5f-② 선형 Cox score C-index (TabICL 비교 기준)
    covs = ["age", "성별", "tumor size (cm)", "DOI (mm)", "differentiation",
            "budding_01vs23", "PNI", "LVI", "RM", "TIL", "TSR",
            "WPOI5_2tier", "HPV/P16_1", "HPV/P16_2", "CCRT_bin"]
    lin_rows = []
    for y in Y_LABELS:
        r = linear_score_cindex_compare(d, y, covs)
        if r:
            lin_rows.append(r)
    if lin_rows:
        lin = pd.DataFrame(lin_rows)
        lin.to_csv(os.path.join(OUT_DIR, "exp1_linear_cindex.csv"), index=False,
                   encoding="utf-8-sig")
        print("\n=== 5f-② 선형 Cox score C-index (TabICL과 비교 기준) ===")
        print(lin.to_string(index=False))

    print("\n=== Score ↔ 생존 (C-index, uncensored Spearman) ===")
    print(out.to_string(index=False))
    print("\n=== RMST(t=36) by score quartile ===")
    if rmst_all:
        print(pd.DataFrame(rmst_all).to_string(index=False))
    if rmst_diff_all:
        print("\n=== RMST 차이 (최고위험 Q_high − 최저위험 Q_low, bootstrap 1,000) ===")
        print(pd.DataFrame(rmst_diff_all).to_string(index=False))
    print(f"\n저장: {OUT_DIR}/exp1_score_survival.csv, exp1_rmst.csv, "
          f"exp1_rmst_diff.csv, score_survival_*.png")


if __name__ == "__main__":
    main()