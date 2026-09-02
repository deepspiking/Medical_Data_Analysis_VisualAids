# -*- coding: utf-8 -*-
"""
수정병기(mNstage/mStage) 재산출 검증 — 정세운 선생님 확정(D1·D3)용 증거
====================================================================
목적:
  'x'→0 정책이 "mNstage 파워 증가 / 전체 stage 연쇄 변화" 기대에 얼마나
  부합하는지 정량화. 원본 Sheet2의 mNstage/mStage 컬럼(의사 제공)을
  Sheet1 규칙으로 재현해 보고:
    1) 재산출값 vs 원본 컬럼 일치율 (규칙 재현 검증)
    2) NX 4명 정책별 시나리오: A=현행 제외(n=129) / B='x'→0으로 mN0 간주(n=133)

Sheet1 규칙:
  mN0: 림프절 전이 없음 (LN count=0)
  mN1: LN 1개 AND 최대크기 ≤10mm
  mN2: LN 2개 또는 크기 >10mm
  mN3: LN≥3 또는 ENE(+) 또는 양측/대측(contra_bilateral)
  mStage I : mT1-2 & mN0
  mStage II: mT3 & mN0 또는 mT1-2 & mN1
  mStage III: mT4 & mN0 또는 mT3 & mN1 또는 mT1-3 & mN2
  mStage IV: mN3 또는 (mT4 & mN1-3)

출력: results/exp1/validation/stage_recompute_check.csv (환자별 상세)
      results/exp1/stage_recompute_summary.txt (요약)
"""
import os
import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE_DIR, "정박사님께 드릴 raw data.xlsx")
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")
os.makedirs(os.path.join(OUT_DIR, "validation"), exist_ok=True)


def mN_from(lnc, size, ene, contra, x0):
    if x0:
        lnc = 0.0 if pd.isna(lnc) else lnc
        size = 0.0 if pd.isna(size) else size
        ene = 0.0 if pd.isna(ene) else ene
        contra = 0.0 if pd.isna(contra) else contra
    if pd.isna(lnc) or pd.isna(size) or pd.isna(ene) or pd.isna(contra):
        return np.nan
    lnc = float(lnc)
    if lnc == 0:
        return 0
    if lnc >= 3 or ene == 1 or contra == 1:
        return 3
    if lnc == 1 and size <= 10:
        return 1
    if lnc == 2 or size > 10:
        return 2
    return np.nan


def mStage_from(mt, mn):
    if pd.isna(mt) or pd.isna(mn):
        return np.nan
    mt, mn = int(mt), int(mn)
    if mt <= 2 and mn == 0:
        return 1
    if (mt == 3 and mn == 0) or (mt <= 2 and mn == 1):
        return 2
    if (mt == 4 and mn == 0) or (mt == 3 and mn == 1) or (mt <= 3 and mn == 2):
        return 3
    if mn == 3 or (mt == 4 and mn in (1, 2, 3)):
        return 4
    return np.nan


def num(s):
    return pd.to_numeric(s.astype(str).str.strip().replace("x", "nan"), errors="coerce")


def main():
    df = pd.read_excel(SRC, sheet_name="Sheet2")
    out = pd.DataFrame({
        "연구번호": df["연구번호"], "N stage": df["N stage"],
        "LN_count": num(df["LN meta count"]), "largest_node": num(df["largest node (mm)"]),
        "ENE": num(df["ENE"]), "contra": num(df["contra_bilateral"]),
        "mTstage_col": num(df[" mTstage"]), "mNstage_col": num(df["mNstage"]),
        "mStage_col": num(df["mStage"]),
    })

    for x0 in (False, True):
        tag = "x0" if x0 else "keepNA"
        out[f"mN_re_{tag}"] = [mN_from(r.LN_count, r.largest_node, r.ENE, r.contra, x0)
                               for r in out.itertuples()]
        out[f"mStage_re_{tag}"] = [mStage_from(mt, mn)
                                   for mt, mn in zip(out[f"mTstage_col"], out[f"mN_re_{tag}"])]

    # 1) 규칙 재현 검증: 원본 컬럼 존재(129명) 구간에서 재산출(x0 아님) 일치율
    m = out.dropna(subset=["mNstage_col"]).copy()
    m["mn_agree"] = m["mN_re_keepNA"] == m["mNstage_col"]
    m["ms_agree"] = m["mStage_re_keepNA"] == m["mStage_col"]
    print("== 규칙 재현 vs 원본 컬럼 (mNstage/mStage 존재 129명) ==")
    print(f"mNstage 일치: {int(m.mn_agree.sum())}/129 "
          f"(불일치 {int((~m.mn_agree).sum())})")
    print(f"mStage  일치: {int(m.ms_agree.sum())}/129 "
          f"(불일치 {int((~m.ms_agree).sum())})")
    bad = m[~m.mn_agree | ~m.ms_agree]
    if len(bad):
        print("\n불일치 환자 목록:")
        print(bad[["연구번호", "N stage", "LN_count", "largest_node", "ENE", "contra",
                   "mTstage_col", "mNstage_col", "mN_re_keepNA",
                   "mStage_col", "mStage_re_keepNA"]].to_string(index=False))
    else:
        print("불일치 없음 — Sheet1 규칙이 원본 컬럼을 정확히 재현.")

    # 2) 시나리오: A(현행, NX 제외 n=129) vs B(x0, NX→mN0 포함 n=133)
    print("\n== 시나리오별 mStage 분포 ==")
    print("A(현행, NX 제외):")
    print(out["mStage_col"].value_counts(dropna=False).sort_index().to_string())
    print("B('x'→0, NX 포함 재산출):")
    print(out["mStage_re_x0"].value_counts(dropna=False).sort_index().to_string())
    print("\nNX 4명 처리 (B에서 mN0/mStage 부여 여부):")
    nx = out[out["mNstage_col"].isna()]
    print(nx[["연구번호", "N stage", "LN_count", "largest_node", "ENE", "contra",
              "mTstage_col", "mN_re_x0", "mStage_re_x0"]].to_string(index=False))

    # 환자별 상세 저장
    out.to_csv(os.path.join(OUT_DIR, "validation", "stage_recompute_check.csv"),
               index=False, encoding="utf-8-sig")
    print(f"\n저장: results/exp1/validation/stage_recompute_check.csv")


if __name__ == "__main__":
    main()
