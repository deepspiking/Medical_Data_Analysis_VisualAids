# -*- coding: utf-8 -*-
"""
전처리: raw data.xlsx Sheet2 → preprocessed_data.csv (+ 컬럼 사전)

원칙
  1) 'x'/'?' 마커 = 미기재/해당없음. 행을 버리지 않고
     {col}_val (숫자 판독 가능 시 값, 아니면 NaN) + {col}_known (판독 여부 0/1) 생성
  2) 서수(ordinal) 카테고리 → 숫자 유지
  3) 명목(nominal) 카테고리 → one-hot (drop_first, 기준 레벨 제거)
  4) Y label: PFS/DSS/LRRFS + time(수술일 기준 월)
  5) CCRT: 날짜 존재 여부 → CCRT_bin

출력
  - preprocessed_data.csv     : 분석용 데이터
  - preprocessed_columns.csv  : 컬럼 사전 (역할/타입/원본 컬럼/설명)
"""
import os
import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE_DIR, "정박사님께 드릴 raw data.xlsx")
OUT_CSV = os.path.join(BASE_DIR, "preprocessed_data.csv")
OUT_COL = os.path.join(BASE_DIR, "preprocessed_columns.csv")

N_STAGE_NX = {7: np.nan}          # 7 = NX (pTNM T2NX 확인). 6 = N3b(유지)
X_COLS = ["Node ", "bone invasion", "depth of bone invasion (mm)",
          "LN meta count", "contra_bilateral", "largest node (mm)",
          "LN tumor size (mm)", "ENE", "N stage"]

ORDINAL_COLS = ["age", "TIL", "differentiation", "T stage", "ajcc8th_STAGE",
                " mTstage", "mNstage", "mStage", "tumor size (cm)", "DOI (mm)"]
BINARY_COLS = ["성별", "budding_01vs23", "PNI", "LVI", "RM", "TSR",
               "WPOI5_2tier", "PD_01vs2"]
NOMINAL_COLS = ["subsite", "Tx_3tier", "HPV/P16"]
Y_LABELS = ["PFS", "DSS", "OS", "LRRFS"]          # OS 추가 (raw: death/Death_month)
Y_SOURCE = {"PFS": ("PFS", "PFS_month"),
            "DSS": ("DSS", "DSS_month"),
            "OS": ("death", "Death_month"),
            "LRRFS": ("LRRFS", "LRRFS_month")}


def _to_numeric_or_nan(s):
    return pd.to_numeric(s.astype(str).str.replace("?", "nan", regex=False),
                         errors="coerce")


def main():
    xl = pd.ExcelFile(SRC)
    df = xl.parse("Sheet2").copy()

    # ---- Y label 정리 (PFS/DSS/OS/LRRFS) ----
    for y in Y_LABELS:
        src_ev, src_tm = Y_SOURCE[y]
        df[f"{y}_event"] = df[src_ev].astype(int)
        df[f"{y}_time"] = df[src_tm].astype(float)

    # ---- 'x' 마커 컬럼 → {col}_val + {col}_known ----
    for c in X_COLS:
        val = _to_numeric_or_nan(df[c])
        known = (~val.isna()).astype(int)
        df[f"{c.strip()}_val"] = val
        df[f"{c.strip()}_known"] = known
    # N stage: 7=NX 처리 (위에서 7은 숫자로 남으므로 여기서 NaN)
    df["N stage_val"] = df["N stage_val"].replace(N_STAGE_NX)

    # ---- CCRT: 날짜 존재 여부 ----
    ccrt = pd.to_datetime(df["1st OP 후 CCRT 날짜"], errors="coerce")
    df["CCRT_bin"] = ccrt.notna().astype(int)

    # ---- 서수: 숫자 유지 ----
    ord_out = pd.DataFrame(index=df.index)
    for c in ORDINAL_COLS:
        ord_out[c.strip()] = _to_numeric_or_nan(df[c])

    # ---- 이분형: 숫자 유지 (0/1) ----
    bin_out = pd.DataFrame(index=df.index)
    for c in BINARY_COLS:
        bin_out[c.strip()] = _to_numeric_or_nan(df[c])

    # ---- 명목: one-hot (drop_first) ----
    oh_out = pd.DataFrame(index=df.index)
    for c in NOMINAL_COLS:
        dummies = pd.get_dummies(df[c].astype(str), prefix=c.strip(), drop_first=True)
        oh_out = pd.concat([oh_out, dummies], axis=1)

    # ---- 수정병기 NaN 4명 (NX) ----
    # mNstage/mStage 는 ORDINAL_COLS에 이미 포함, NX 4명은 NaN 유지.

    out = pd.concat([
        df[["연구번호"]],
        ord_out,
        bin_out,
        oh_out,
        df[[f"{y}_{s}" for y in Y_LABELS for s in ("event", "time")]],
        df[["CCRT_bin"]],
        df[[c for c in df.columns if c.endswith("_val") or c.endswith("_known")]],
    ], axis=1)

    # 중복 컬럼 정리 (N stage_val 등은 이미 포함, ord_out과 겹치지 않도록)
    out = out.loc[:, ~out.columns.duplicated()]
    out.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")

    # ---- 컬럼 사전 ----
    doc = []
    def add(c, role, origin, desc):
        doc.append({"column": c, "role": role, "origin": origin, "description": desc})
    for c in ord_out.columns:
        add(c, "ordinal", c, "서수형 수치 (높을수록 중증/진행)")
    for c in bin_out.columns:
        add(c, "binary", c, "이분형 0/1")
    for c in oh_out.columns:
        add(c, "onehot", c.split("_")[0], "명목 one-hot (drop_first)")
    for y in Y_LABELS:
        add(f"{y}_event", "y_label", y, f"{y} 사건 여부 (0/1)")
        add(f"{y}_time", "y_label", f"{y}_month", f"{y} 시간 (수술일 기준 월)")
    add("CCRT_bin", "binary", "1st OP 후 CCRT 날짜", "CCRT 시행 여부 (날짜 존재=1)")
    for c in df.columns:
        if c.endswith("_val"):
            add(c, "x_value", c[:-4], "원본 값 ('x'/'?'는 NaN)")
        elif c.endswith("_known"):
            add(c, "x_flag", c[:-6], "값 판독 가능 여부 (1=있음, 0='x'/미기재)")
    pd.DataFrame(doc).to_csv(OUT_COL, index=False, encoding="utf-8-sig")

    print(f"저장 완료: {OUT_CSV}")
    print(f"컬럼 사전: {OUT_COL}")
    print(f"shape: {out.shape}")
    n_nx = df["mStage"].isna().sum()
    print(f"mStage/mNstage NaN (NX): {n_nx}명 — mStage 비교에서 제외, 모델에는 _known으로 유지")


if __name__ == "__main__":
    main()