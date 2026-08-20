# -*- coding: utf-8 -*-
"""
추가 검증 도구 (실험 1-5 외 보강)
------------------------------------------------------------------
실험 3(다변량) 모델 A/B/C의 진단·교정·공선성 검증 + 추가 CV:

  T1. VIF (다중공선성)          : 수정병기/기존병기 + 공변량 공선성 진단
  T2. Schoenfeld 잔차           : Cox 비례위험 가정 검정 (lifelines)
  T3. Calibration plot          : 24개월 예측 vs 관측 (gowun 브랜치 스타일)
  T4. DCA (Decision Curve)      : 순이익 곡선 (gowun 브랜치 스타일)
  T5. LOOCV C-index             : S1/S2 보강 (n=133 → Leave-One-Out 가능)
  T6. 시간종속 Brier score      : 24개월 예측 정확도
  T7. 시드 안정성 요약          : 4개 시드(42/123/2026/777) 결과 변동 보고
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from experiment_utils import X_MULTI_BASE, fast_cindex

plt.rcParams["font.family"] = ["AppleGothic", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_CSV = os.path.join(BASE_DIR, "preprocessed_data.csv")
RES = os.path.join(BASE_DIR, "results", "exp1")
OUT = os.path.join(RES, "validation")
os.makedirs(OUT, exist_ok=True)

DPI = 300
SEED = 42
T_HORIZON = 24.0

COVS = ["age", "성별", "tumor size (cm)", "DOI (mm)", "differentiation",
        "budding_01vs23", "PNI", "LVI", "RM", "TIL", "TSR",
        "WPOI5_2tier", "HPV/P16_1", "HPV/P16_2", "CCRT_bin"]
STAGE = {"B_ajcc8th": "ajcc8th_STAGE", "C_modified": "mStage"}


def load_prep():
    d = pd.read_csv(DATA_CSV)
    return d, X_MULTI_BASE


def _save(fig, fname):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[검증] {fname}")


# ----------------------------------------------------------------------------
# T1. VIF (다중공선성)
# ----------------------------------------------------------------------------
def vif_table(d, X, stage_col):
    """VIF (numpy 구현 — statsmodels는 np.MachAr 호환 문제로 사용 불가)."""
    cols = X + [stage_col]
    dfx = d[cols].dropna()
    Xm = dfx[cols].values.astype(float)
    Xc = np.column_stack([np.ones(len(Xm)), Xm])   # 상수항 추가
    rows = []
    for j, c in enumerate(cols):
        col_idx = j + 1
        y = Xc[:, col_idx]
        Xo = np.delete(Xc, col_idx, axis=1)
        try:
            beta = np.linalg.lstsq(Xo, y, rcond=None)[0]
            resid = y - Xo @ beta
            r2 = 1 - (resid**2).sum() / ((y - y.mean())**2).sum()
            v = 1 / (1 - r2) if r2 < 1 else np.inf
        except Exception:
            v = np.nan
        rows.append({"var": c, "VIF": round(float(v), 2)})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# T2. Schoenfeld 잔차 (비례위험 가정)
# ----------------------------------------------------------------------------
def schoenfeld_test(dfx, covs, y_time, y_event):
    from lifelines import CoxPHFitter
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(dfx[["time", "event"] + covs], duration_col="time", event_col="event")
    try:
        test = cph.check_assumptions(dfx[["time", "event"] + covs],
                                     show_plots=False, p_value_print=False)
    except Exception:
        return None
    rows = []
    for v, res in test.results.items():
        rows.append({"var": v, "test_stat": round(float(res.test_statistic), 3),
                     "p_value": round(float(res.p_value), 4)})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# T3. Calibration plot (24개월)
# ----------------------------------------------------------------------------
def calibration(dfx, covs, y_time, y_event, t=T_HORIZON):
    from lifelines import CoxPHFitter
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(dfx[["time", "event"] + covs], duration_col="time", event_col="event")
    sf = cph.predict_survival_function(dfx[covs], times=[t])
    pred_risk = (1.0 - sf.values[0]).ravel()
    obs = ((dfx[y_event].values == 1) & (dfx[y_time].values <= t)).astype(float)
    dfc = pd.DataFrame({"pred": pred_risk, "obs": obs}).dropna()
    if len(dfc) < 20:
        return None, None
    dfc["grp"] = pd.qcut(dfc["pred"], 5, labels=False, duplicates="drop")
    g = dfc.groupby("grp").mean()
    return g["pred"].values, g["obs"].values


# ----------------------------------------------------------------------------
# T4. DCA (Decision Curve, 24개월)
# ----------------------------------------------------------------------------
def dca(dfx, covs, y_time, y_event, t=T_HORIZON):
    from lifelines import CoxPHFitter
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(dfx[["time", "event"] + covs], duration_col="time", event_col="event")
    sf = cph.predict_survival_function(dfx[covs], times=[t])
    risk = (1.0 - sf.values[0]).ravel()
    obs = ((dfx[y_event].values == 1) & (dfx[y_time].values <= t)).astype(float)
    n = len(obs)
    event_rate = obs.mean()
    thresholds = np.arange(0.05, 0.80, 0.01)
    nb_model, nb_all, nb_none = [], [], []
    for pt in thresholds:
        tp = ((risk >= pt) & (obs == 1)).sum()
        fp = ((risk >= pt) & (obs == 0)).sum()
        nb_model.append(tp / n - fp / n * (pt / (1 - pt)))
        nb_all.append(event_rate - (1 - event_rate) * (pt / (1 - pt)))
        nb_none.append(0.0)
    return thresholds, nb_model, nb_all, nb_none


# ----------------------------------------------------------------------------
# T5. LOOCV C-index
# ----------------------------------------------------------------------------
def _loocv_fit_one(args):
    i, dfx, covs = args
    from lifelines import CoxPHFitter
    tr = dfx.drop(index=dfx.index[i])   # dropna 후 인덱스 레이블 사용
    try:
        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(tr[["time", "event"] + covs], duration_col="time",
                event_col="event")
        return i, float(cph.predict_partial_hazard(
            dfx.iloc[[i]][covs]).values[0])
    except Exception:
        return i, np.nan


def loocv_cindex(dfx, covs, y_time, y_event):
    import multiprocessing as mp
    try:
        mp.set_start_method("fork", force=True)
    except RuntimeError:
        pass
    T = dfx[y_time].values.astype(float)
    E = dfx[y_event].values.astype(bool)
    n = len(dfx)
    tasks = [(i, dfx, covs) for i in range(n)]
    nproc = min(8, mp.cpu_count())
    with mp.Pool(processes=nproc) as pool:
        results = pool.map(_loocv_fit_one, tasks)
    risks = np.full(n, np.nan)
    for i, r in results:
        risks[i] = r
    valid = ~np.isnan(risks)
    if valid.sum() < 20:
        return np.nan
    return float(fast_cindex(risks[valid], T[valid], E[valid]))


# ----------------------------------------------------------------------------
# T6. 시간종속 Brier score (24개월, IPCW)
# ----------------------------------------------------------------------------
def brier_score(dfx, covs, y_time, y_event, t=T_HORIZON):
    from lifelines import CoxPHFitter, KaplanMeierFitter
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(dfx[["time", "event"] + covs], duration_col="time", event_col="event")
    sf = cph.predict_survival_function(dfx[covs], times=[t])
    surv = sf.values[0].ravel()
    T = dfx[y_time].values.astype(float)
    E = dfx[y_event].values.astype(bool)
    km_c = KaplanMeierFitter().fit(T, event_observed=(~E))
    g = km_c.survival_function_at_times(pd.Series(np.minimum(T, t))).values.ravel()
    g = np.clip(g, 1e-3, 1.0)
    ind = (T <= t) & E
    ind_c = T > t
    w = np.zeros(len(T))
    w[ind] = 1.0 / g[ind]
    w[ind_c] = 1.0 / g[ind_c]
    bs = np.mean(w * (ind.astype(float) - (1.0 - surv)) ** 2)
    return float(bs)


# ----------------------------------------------------------------------------
# T7. 시드 안정성 요약 (실험 2/3 결과)
# ----------------------------------------------------------------------------
def seed_stability():
    uni = os.path.join(RES, "exp1_univariate.csv")
    multi = os.path.join(RES, "exp1_multivariate.csv")
    if not os.path.exists(uni):
        return None, None
    d = pd.read_csv(uni)
    summ = d.groupby(["label", "pair"])[["dC", "NRI"]].agg(
        ["mean", "std", "min", "max"]).round(3)
    m = pd.read_csv(multi) if os.path.exists(multi) else None
    return summ, m


def main():
    d, X = load_prep()
    print(f"검증 도구 | 데이터 {d.shape[0]}행, X {len(X)}개")

    # T1. VIF (모델 B/C)
    vif_rows = []
    for tag, sc in STAGE.items():
        v = vif_table(d, X, sc)
        v.insert(0, "model", tag)
        vif_rows.append(v)
    vif = pd.concat(vif_rows, ignore_index=True)
    vif.to_csv(os.path.join(OUT, "T1_VIF.csv"), index=False, encoding="utf-8-sig")
    print("\n=== T1. VIF (10 초과 = 공선성 의심) ===")
    print(vif[vif["VIF"] > 5].to_string(index=False))

    # T2-T4, T6: y별 다변량 모델 (A 기준)
    for y in ["PFS", "DSS", "LRRFS"]:
        y_time, y_event = f"{y}_time", f"{y}_event"
        dfx = d[X + [y_time, y_event]].dropna().rename(
            columns={y_time: "time", y_event: "event"})
        # T2. Schoenfeld
        sch = schoenfeld_test(dfx, X, "time", "event")
        if sch is not None:
            sch.insert(0, "y", y)
            sch.to_csv(os.path.join(OUT, f"T2_Schoenfeld_{y}.csv"), index=False,
                       encoding="utf-8-sig")
            sig = sch[sch["p_value"] < 0.05]
            print(f"\n=== T2. Schoenfeld [{y}] 비례위험 위반: {len(sig)}개 ===")
            print(sig.to_string(index=False) if len(sig) else "(없음)")

        # T3. Calibration plot
        pred, obs = calibration(dfx, X, "time", "event")
        if pred is not None:
            fig, ax = plt.subplots(figsize=(6.5, 6))
            ax.plot(pred, obs, "o-", color="#2ca02c", linewidth=2, markersize=8)
            ax.plot([0, 1], [0, 1], "k--", alpha=0.6, label="Perfect")
            ax.set_xlabel(f"Predicted {y} risk @ {int(T_HORIZON)}mo")
            ax.set_ylabel("Observed proportion")
            ax.set_title(f"Calibration - {y} (모델 A)")
            ax.grid(alpha=0.3)
            ax.legend()
            _save(fig, f"T3_calibration_{y}.png")

        # T4. DCA plot
        thr, nb_m, nb_all, nb_none = dca(dfx, X, "time", "event")
        fig, ax = plt.subplots(figsize=(7, 5.5))
        ax.plot(thr, nb_m, color="#d62728", linewidth=2, label="Model A")
        ax.plot(thr, nb_all, color="gray", linestyle="--", label="Treat all")
        ax.plot(thr, nb_none, color="black", linestyle=":", label="Treat none")
        ax.set_xlabel("Threshold probability")
        ax.set_ylabel("Net benefit")
        ax.set_title(f"DCA - {y} @ {int(T_HORIZON)}mo (모델 A)")
        ax.legend()
        ax.grid(alpha=0.3)
        _save(fig, f"T4_DCA_{y}.png")

        # T5. LOOCV (A/B/C)
        for tag, sc in STAGE.items():
            dfx2 = d[X + [sc, y_time, y_event]].dropna().rename(
                columns={y_time: "time", y_event: "event"})
            c = loocv_cindex(dfx2, X + [sc], "time", "event")
            print(f"=== T5. LOOCV C-index [{y}] {tag}: {c:.4f} ===")

        # T6. Brier
        bs = brier_score(dfx, X, "time", "event")
        print(f"=== T6. Brier score [{y}] 모델 A @ {int(T_HORIZON)}mo: {bs:.4f} ===")

    # T7. 시드 안정성
    summ, m = seed_stability()
    if summ is not None:
        summ.to_csv(os.path.join(OUT, "T7_seed_stability.csv"), encoding="utf-8-sig")
        print("\n=== T7. 시드 안정성 (dC, NRI 변동) ===")
        print(summ.to_string())

    print(f"\n검증 결과 저장: {OUT}")


if __name__ == "__main__":
    main()