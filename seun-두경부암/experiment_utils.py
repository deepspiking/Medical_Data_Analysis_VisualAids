# -*- coding: utf-8 -*-
"""
공통 유틸: 실험 2·3에서 공유하는 함수 모음
- fast_cindex (벡터화 Harrell C)
- 전처리 CSV 로드 / X_MULTI_BASE / 희소 범주 병합
- Cox fit (penalizer + None 반환)
"""
import os
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PREP_CSV = os.path.join(BASE_DIR, "preprocessed_data.csv")

# 다변량 X (실험 3): 수정병기 구성요소 제외 + 공선성 제거
#   - CCRT_bin 제외: Tx_3tier_2와 완전 공선성 (VIF=∞, CCRT_bin == Tx_3tier_2)
#   - 희소 범주는 penalizer(0.1)로 안정화 (제거하지 않음 — 전체 변수 유지)
X_MULTI_BASE = ["age", "TIL", "differentiation", "tumor size (cm)", "DOI (mm)",
                "성별", "budding_01vs23", "PNI", "LVI", "RM", "TSR",
                "WPOI5_2tier",
                "subsite_2", "subsite_3", "subsite_4", "subsite_5", "subsite_6",
                "subsite_7", "subsite_8",
                "Tx_3tier_1", "Tx_3tier_2",
                "HPV/P16_1", "HPV/P16_2"]


def fast_cindex(S, T, E):
    """벡터화 Harrell C-index. 높은 score = 나쁜 예후."""
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


def fit_multi(dfx, covs):
    """Cox fit (ridge penalizer). 실패 시 None 반환."""
    cph = CoxPHFitter(penalizer=0.1)
    try:
        cph.fit(dfx[["time", "event"] + list(covs)], duration_col="time",
                event_col="event")
        return cph
    except Exception:
        return None


def build_multi_df(label):
    """실험 3용 데이터: 전처리 CSV + X_MULTI_BASE + stage 컬럼."""
    d = pd.read_csv(PREP_CSV)
    y_time, y_event = f"{label}_time", f"{label}_event"
    keep = X_MULTI_BASE + [y_time, y_event, "ajcc8th_STAGE", "mStage"]
    m = d[keep].rename(columns={y_time: "time", y_event: "event"})
    return m.dropna().reset_index(drop=True)
