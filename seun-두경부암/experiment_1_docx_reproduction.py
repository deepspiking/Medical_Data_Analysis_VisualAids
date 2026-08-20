# -*- coding: utf-8 -*-
"""
실험 1: oralSCC_table.docx (정세운 선생님 기존 분석) 재현 + Bootstrap 검증
------------------------------------------------------------------
docx Table 1: PFS/DSS/OS/LRRFS 각각 univariate + multivariate Cox
  변수        : tumor size(cm), DOI(mm), LN meta count, LN tumor size(mm),
                ENE, PD-HC(PD_01vs2), PNI, LVI, WPOI5, budding
  multivariate 구성 (docx에 명시된 최종 모델):
    PFS   : tumor size + LN tumor size + PD-HC
    DSS   : tumor size + LN tumor size + PD-HC + LVI
    OS    : PD-HC
    LRRFS : LN tumor size
검증: 재현 HR/CI/p vs docx 값 방향·유의성 일치 + bootstrap 1,000회 안정성

docx 참조값 (Table 1에서 추출):
  PFS  uni: size 1.296(1.123-1.495) / DOI 1.054(1.029-1.080) / #LN 1.061(1.028-1.095)
            / dep 1.047(1.011-1.084) / ENE 2.332(1.328-4.094) / PD-HC 3.867(2.284-6.546)
            / PNI 2.112(1.139-3.915) / LVI 2.361(1.311-4.250)
  PFS  multi: size 1.430(1.063-1.924) / dep 1.083(1.024-1.145) / PD-HC 2.602(1.186-5.711)
  DSS  uni: size 1.577(1.342-1.854) / DOI 1.081(1.050-1.112) / #LN 1.096(1.051-1.144)
            / dep 1.052(1.008-1.099) / ENE 3.633(1.681-7.853) / PD-HC 8.114(3.804-17.309)
            / PNI 3.067(1.347-6.980) / LVI 3.481(1.571-7.711)
  DSS  multi: size 1.992(1.348-2.944) / dep 1.094(1.007-1.189) / PD-HC 13.115(2.517-68.348)
              / LVI 6.158(1.784-21.253)
  OS   multi: PD-HC 3.372(1.631-6.970)
  LRRFS uni: size 1.098(0.922-1.306) / DOI 1.039(1.009-1.071) / dep 1.057(1.016-1.100)
            / ENE 2.260(1.168-4.372) / PD-HC 3.281(1.759-6.118)
  LRRFS multi: dep 1.073(1.014-1.135)
"""
import os
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(BASE_DIR, "정박사님께 드릴 raw data.xlsx")
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")
os.makedirs(OUT_DIR, exist_ok=True)

N_BOOT = 1_000
SEED = 42

# docx Table 1 변수명 → raw data 컬럼
VAR_MAP = {
    "tumor_size": "tumor size (cm)",
    "DOI": "DOI (mm)",
    "LN_count": "LN meta count",
    "LN_deposit": "LN tumor size (mm)",
    "ENE": "ENE",
    "PD_HC": "PD_01vs2",
    "PNI": "PNI",
    "LVI": "LVI",
    "WPOI5": "WPOI5_2tier",
    "budding": "budding_01vs23",
}
UNI_VARS = list(VAR_MAP.keys())

# docx multivariate 최종 모델 구성 (y label → 변수)
MULTI_MODELS = {
    "PFS": ["tumor_size", "LN_deposit", "PD_HC"],
    "DSS": ["tumor_size", "LN_deposit", "PD_HC", "LVI"],
    "OS": ["PD_HC"],
    "LRRFS": ["LN_deposit"],
}

# docx 참조값 (변수 → (HR, lo, hi)) — univariate
DOCX_UNI = {
    "PFS": {"tumor_size": (1.296, 1.123, 1.495), "DOI": (1.054, 1.029, 1.080),
            "LN_count": (1.061, 1.028, 1.095), "LN_deposit": (1.047, 1.011, 1.084),
            "ENE": (2.332, 1.328, 4.094), "PD_HC": (3.867, 2.284, 6.546),
            "PNI": (2.112, 1.139, 3.915), "LVI": (2.361, 1.311, 4.250)},
    "DSS": {"tumor_size": (1.577, 1.342, 1.854), "DOI": (1.081, 1.050, 1.112),
            "LN_count": (1.096, 1.051, 1.144), "LN_deposit": (1.052, 1.008, 1.099),
            "ENE": (3.633, 1.681, 7.853), "PD_HC": (8.114, 3.804, 17.309),
            "PNI": (3.067, 1.347, 6.980), "LVI": (3.481, 1.571, 7.711)},
    "LRRFS": {"tumor_size": (1.098, 0.922, 1.306), "DOI": (1.039, 1.009, 1.071),
              "LN_deposit": (1.057, 1.016, 1.100), "ENE": (2.260, 1.168, 4.372),
              "PD_HC": (3.281, 1.759, 6.118)},
}
DOCX_MULTI = {
    "PFS": {"tumor_size": (1.430, 1.063, 1.924), "LN_deposit": (1.083, 1.024, 1.145),
            "PD_HC": (2.602, 1.186, 5.711)},
    "DSS": {"tumor_size": (1.992, 1.348, 2.944), "LN_deposit": (1.094, 1.007, 1.189),
            "PD_HC": (13.115, 2.517, 68.348), "LVI": (6.158, 1.784, 21.253)},
    "OS": {"PD_HC": (3.372, 1.631, 6.970)},
    "LRRFS": {"LN_deposit": (1.073, 1.014, 1.135)},
}

Y_OUTCOMES = {  # outcome → (time_col, event_col)
    "PFS": ("PFS_month", "PFS"),
    "DSS": ("DSS_month", "DSS"),
    "OS": ("Death_month", "death"),
    "LRRFS": ("LRRFS_month", "LRRFS"),
}


def load_outcome(outcome):
    xl = pd.ExcelFile(DATA_PATH)
    df = xl.parse("Sheet2")
    time_col, event_col = Y_OUTCOMES[outcome]
    keep = ["연구번호"] + list(VAR_MAP.values()) + [time_col, event_col]
    d = df[keep].copy()
    d = d.rename(columns={time_col: "time", event_col: "event"})
    d["event"] = d["event"].astype(int)
    d["time"] = pd.to_numeric(d["time"], errors="coerce").astype(float)
    for v in VAR_MAP.values():
        d[v] = pd.to_numeric(d[v].astype(str).str.replace("x", "nan", regex=False),
                             errors="coerce")
    d = d.dropna(subset=["time"]).reset_index(drop=True)
    return d


def fit_cox(d, vars_list):
    dfx = d[["time", "event"] + [VAR_MAP[v] for v in vars_list]].dropna()
    cph = CoxPHFitter(penalizer=0.0)
    cph.fit(dfx, duration_col="time", event_col="event")
    out = {}
    for v in vars_list:
        col = VAR_MAP[v]
        ci = cph.confidence_intervals_.loc[col]
        out[v] = (float(np.exp(cph.params_[col])),
                  float(np.exp(ci["95% lower-bound"])),
                  float(np.exp(ci["95% upper-bound"])),
                  float(cph.summary.loc[col, "p"]))
    return out


def _bootstrap_one(args):
    d, v, n_boot, seed = args
    rng = np.random.RandomState(seed + (abs(hash(v)) % 1000))  # 고정 오프셋
    dfx = d[["time", "event", VAR_MAP[v]]].dropna()
    n = len(dfx)
    col = VAR_MAP[v]
    hr_boot = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.randint(0, n, n)
        boot = dfx.iloc[idx]
        try:
            cph = CoxPHFitter(penalizer=0.0)
            cph.fit(boot, duration_col="time", event_col="event")
            hr_boot[b] = float(np.exp(cph.params_[col]))
        except Exception:
            hr_boot[b] = np.nan
    ok = hr_boot[~np.isnan(hr_boot)]
    if len(ok) == 0:
        return None
    return {"var": v, "boot_HR_lo": float(np.percentile(ok, 2.5)),
            "boot_HR_hi": float(np.percentile(ok, 97.5)),
            "boot_n": len(ok)}


def bootstrap_hr(d, vars_list, n_boot=N_BOOT, seed=SEED):
    import multiprocessing as mp
    tasks = [(d, v, n_boot, seed) for v in vars_list]
    nproc = min(8, mp.cpu_count())
    with mp.Pool(processes=nproc) as pool:
        results = pool.map(_bootstrap_one, tasks)
    rows = [r for r in results if r is not None]
    print(f"  [bootstrap] {len(rows)} 변수 완료 ({n_boot:,}회, {nproc}코어)", flush=True)
    return pd.DataFrame(rows)


def _init_parallel():
    import multiprocessing as mp
    try:
        mp.set_start_method("fork", force=True)
    except RuntimeError:
        pass
    print(f"[실험 1] 병렬 모드: fork, {min(8, mp.cpu_count())} 프로세스", flush=True)


def main():
    import time
    _init_parallel()
    t0 = time.time()
    out_rows = []
    boot_rows = []
    for outcome, (time_col, event_col) in Y_OUTCOMES.items():
        print(f"[실험 1] {outcome} 시작...", flush=True)
        d = load_outcome(outcome)
        uni = fit_cox(d, UNI_VARS)
        for v, (hr, lo, hi, p) in uni.items():
            ref = DOCX_UNI.get(outcome, {}).get(v)
            out_rows.append({"outcome": outcome, "type": "univariate", "var": v,
                             "HR": hr, "lo": lo, "hi": hi, "p": p,
                             "docx_HR": ref[0] if ref else np.nan,
                             "docx_lo": ref[1] if ref else np.nan,
                             "docx_hi": ref[2] if ref else np.nan})
        multi_vars = MULTI_MODELS[outcome]
        multi = fit_cox(d, multi_vars)
        for v, (hr, lo, hi, p) in multi.items():
            ref = DOCX_MULTI.get(outcome, {}).get(v)
            out_rows.append({"outcome": outcome, "type": "multivariate", "var": v,
                             "HR": hr, "lo": lo, "hi": hi, "p": p,
                             "docx_HR": ref[0] if ref else np.nan,
                             "docx_lo": ref[1] if ref else np.nan,
                             "docx_hi": ref[2] if ref else np.nan})
        br = bootstrap_hr(d, multi_vars)
        br.insert(0, "outcome", outcome)
        br.insert(1, "type", "multivariate")
        boot_rows.append(br)
        print(f"[실험 1] {outcome} 완료 ({time.time()-t0:.0f}s)", flush=True)

    res = pd.DataFrame(out_rows)
    res["direction_match"] = np.sign(res["HR"] - 1) == np.sign(res["docx_HR"] - 1)
    res["sig_match"] = (res["p"] < 0.05) == \
        ((res["docx_lo"] > 1) | (res["docx_hi"] < 1))
    # 비교 가능한 행만: docx 참조값 존재 + 재현 HR 존재
    cmp = res[res["docx_HR"].notna() & res["HR"].notna()].copy()
    cmp["direction_match"] = np.sign(cmp["HR"] - 1) == np.sign(cmp["docx_HR"] - 1)
    cmp["sig_match"] = (cmp["p"] < 0.05) == \
        ((cmp["docx_lo"] > 1) | (cmp["docx_hi"] < 1))
    res.to_csv(os.path.join(OUT_DIR, "exp1_docx_reproduction.csv"), index=False,
               encoding="utf-8-sig")

    boot = pd.concat(boot_rows, ignore_index=True)
    boot.to_csv(os.path.join(OUT_DIR, "exp1_docx_bootstrap.csv"), index=False,
                encoding="utf-8-sig")

    print("=== 재현 vs docx (univariate, 참조값 있는 경우만) ===")
    u = cmp[cmp["type"] == "univariate"]
    for _, r in u.iterrows():
        dm = "O" if r["direction_match"] else "X"
        sm = "O" if r["sig_match"] else "X"
        print(f"{r['outcome']:6s} {r['var']:12s} 재현 {r['HR']:6.3f} ({r['lo']:.2f}-{r['hi']:.2f}) "
              f"| docx {r['docx_HR']:6.3f} | 방향{dm} 유의{sm}")
    print("\n=== multivariate (docx 최종 모델, 참조값 있는 경우만) ===")
    m = cmp[cmp["type"] == "multivariate"]
    for _, r in m.iterrows():
        dm = "O" if r["direction_match"] else "X"
        sm = "O" if r["sig_match"] else "X"
        print(f"{r['outcome']:6s} {r['var']:12s} 재현 {r['HR']:7.3f} ({r['lo']:.2f}-{r['hi']:.2f}) "
              f"| docx {r['docx_HR']:7.3f} | 방향{dm} 유의{sm}")
    print(f"\n방향 일치율 (참조값 있는 {len(cmp)}개): "
          f"{(cmp['direction_match'].sum()/len(cmp))*100:.0f}% "
          f"({cmp['direction_match'].sum()}/{len(cmp)})")
    print(f"유의성 일치율 (참조값 있는 {len(cmp)}개): "
          f"{(cmp['sig_match'].sum()/len(cmp))*100:.0f}% "
          f"({cmp['sig_match'].sum()}/{len(cmp)})")
    print(f"\n저장: {OUT_DIR}/exp1_docx_reproduction.csv, exp1_docx_bootstrap.csv")


if __name__ == "__main__":
    main()