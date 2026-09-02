# -*- coding: utf-8 -*-
"""
전수(n=133) 분석용 데이터셋 생성 — NX 4명을 'x'→0(absence=N0) 가정으로 스테이징 채움
=====================================================================================
preprocessed_data.csv(=n=129, mStage/mNstage NaN 4명) → preprocessed_data_full.csv(=n=133)

채우는 환자: 연구번호 2, 37, 81, 126 (N stage 7/x = NX, mNstage/mStage NaN)
정책(D1 가정, 의사 확정 대기): 림프절 미평가(NX)를 'x'→0 규칙에 따라 absence로 간주
  - mNstage = 0 (림프절 전이 없음)
  - mStage  = Sheet1 조합(mTstage & mN0): mT1-2→I, mT3→II, mT4→III
  - N stage_val = 0 (기존/수정 N 병기 비교를 n=133으로 통일)
  - largest node(mm)·LN tumor size(mm)·ENE·LN meta count·contra_bilateral: val=0, known=1
  - ajcc8th_STAGE는 원본 raw에 이미 존재(2,2,3,1) → 변경 없음
"""
import os
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, "preprocessed_data.csv")
OUT = os.path.join(BASE, "preprocessed_data_full.csv")
NX_IDS = [2, 37, 81, 126]

T_MN0_TO_MSTAGE = {1: 1, 2: 1, 3: 2, 4: 3}   # mTstage & mN0 → mStage (Sheet1)


def main():
    d = pd.read_csv(SRC)
    m = d[d["연구번호"].isin(NX_IDS)].copy()
    if len(m) != 4:
        raise SystemExit(f"NX 환자 수가 4가 아님: {len(m)}")
    for _, r in m.iterrows():
        idx = d["연구번호"] == r["연구번호"]
        mt = int(r["mTstage"])
        d.loc[idx, "mNstage"] = 0.0
        d.loc[idx, "mStage"] = float(T_MN0_TO_MSTAGE[mt])
        d.loc[idx, "N stage_val"] = 0.0
        d.loc[idx, "N stage_known"] = 1
        for c in ["largest node (mm)_val", "LN tumor size (mm)_val", "ENE_val",
                  "LN meta count_val", "contra_bilateral_val"]:
            d.loc[idx, c] = 0.0
        for c in ["largest node (mm)_known", "LN tumor size (mm)_known", "ENE_known",
                  "LN meta count_known", "contra_bilateral_known"]:
            d.loc[idx, c] = 1
    d.to_csv(OUT, index=False, encoding="utf-8-sig")
    print("저장:", OUT, d.shape)
    print(d[d["연구번호"].isin(NX_IDS)][
        ["연구번호", "N stage_val", "mNstage", "mTstage", "mStage", "ajcc8th_STAGE"]]
        .to_string(index=False))
    print("mStage 분포:", d["mStage"].value_counts().sort_index().to_dict())
    print("ajcc8th NaN:", int(d["ajcc8th_STAGE"].isna().sum()),
          "| mStage NaN:", int(d["mStage"].isna().sum()),
          "| mNstage NaN:", int(d["mNstage"].isna().sum()))


if __name__ == "__main__":
    main()
