# -*- coding: utf-8 -*-
"""
새 N grouping 제안·검증 (정세운 선생님 피드백 #2 대응)
=====================================================================
입력: preprocessed_data_full.csv (n=133, #52 mNstage=3 정정 반영)
사용 변수: LN meta count, LN tumor size(=전이된 tumor deposit), ENE, contra_bilateral

목표:
  1) 위 4개 변수로 AJCC N과 다른 ≥4군 grouping 제안
  2) 기존 AJCC N 및 선생님 mN과의 C-index 비교 (bootstrap paired ΔC, 95% CI)
  3) 단조성(사건율 저→고) · log-rank p · 군별 구성 보고

출력:
  results/exp1/ngroup_proposal_summary.csv  : 후보 스킴별 지표
  results/exp1/ngroup_proposal_groups.csv   : 제안 스킴 군별 상세
  results/exp1/ngroup_proposal_report.txt   : 요약 리포트
"""
import os
import numpy as np
import pandas as pd
from lifelines.statistics import multivariate_logrank_test

from experiment_utils import fast_cindex

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")
DATA = os.path.join(BASE_DIR, "preprocessed_data_full.csv")
Y_LABELS = ["PFS", "DSS", "OS", "LRRFS"]
N_BOOT = 2000
BOOT_SEED = 42

# AJCC N raw code → ordinal (0=N0,1=N1,2=N2a,3=N2b,4=N2c,6=N3b)
AJCC_MAP = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 6: 5}


def load():
    d = pd.read_csv(DATA)
    D = pd.DataFrame({
        "id": d["연구번호"],
        "ajccN": d["N stage_val"].map(AJCC_MAP),
        "mN": d["mNstage"],
        "cnt": d["LN meta count_val"].fillna(0),
        "dep": d["LN tumor size (mm)_val"].fillna(0),
        "node": d["largest node (mm)_val"].fillna(0),
        "ene": d["ENE_val"].fillna(0),
        "contra": d["contra_bilateral_val"].fillna(0),
    })
    for y in Y_LABELS:
        D[f"{y}_t"] = d[f"{y}_time"]
        D[f"{y}_e"] = d[f"{y}_event"].astype(int)
    return D


def scheme_ndep_4(D, dep_cut=10):
    """제안: N0 / (cnt1-2 & dep<=cut & ENE- & 단측) / (cnt1-2 & dep>cut & ENE- & 단측) / (cnt>=3 or ENE+ or 양측)."""
    g = np.zeros(len(D), dtype=int)
    npos = D.cnt >= 1
    high = (D.ene == 1) | (D.contra == 1) | (D.cnt >= 3)
    g[npos & ~high & (D.dep <= dep_cut)] = 1
    g[npos & ~high & (D.dep > dep_cut)] = 2
    g[npos & high] = 3
    return g


def scheme_ndep_5(D, dep_cut=10, cnt_cut=3):
    """5군 변형: 고위험을 (cnt>=cut & ENE- & 단측) vs (ENE+ or 양측)로 분리."""
    g = np.zeros(len(D), dtype=int)
    npos = D.cnt >= 1
    veryhigh = (D.ene == 1) | (D.contra == 1)
    g[npos & ~veryhigh & (D.cnt < cnt_cut) & (D.dep <= dep_cut)] = 1
    g[npos & ~veryhigh & (D.cnt < cnt_cut) & (D.dep > dep_cut)] = 2
    g[npos & ~veryhigh & (D.cnt >= cnt_cut)] = 3
    g[npos & veryhigh] = 4
    return g


def scheme_count4(D):
    return np.select([D.cnt == 0, D.cnt == 1, D.cnt == 2], [0, 1, 2], 3)


def scheme_count_tier(D):
    return np.select([D.cnt == 0, D.cnt == 1, D.cnt <= 3], [0, 1, 2], 3)


SCHEMES = {
    "AJCC N (기존, 6군)": lambda D: D.ajccN.values.astype(float),
    "mN (선생님 기존, 4군)": lambda D: D.mN.values.astype(float),
    "N-dep 4군 (제안)": lambda D: scheme_ndep_4(D, 10),
    "N-dep 5군 (변형)": lambda D: scheme_ndep_5(D, 10, 3),
    "count-only 4군 (0/1/2/3+)": scheme_count4,
    "count-tier 4군 (0/1/2-3/4+)": scheme_count_tier,
}


def cindex_by_group(g, T, E):
    return fast_cindex(np.asarray(g, dtype=float), T, E)


def bootstrap_dc(g_new, g_ref, T, E, n_boot=N_BOOT, seed=BOOT_SEED):
    n = len(T)
    rng = np.random.RandomState(seed)
    out = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.randint(0, n, n)
        out[b] = cindex_by_group(g_new[idx], T[idx], E[idx]) - \
            cindex_by_group(g_ref[idx], T[idx], E[idx])
    return np.nanpercentile(out, [2.5, 97.5])


def monotone_flags(D, g):
    flags = {}
    rates = {}
    for y in Y_LABELS:
        r = pd.Series(D[f"{y}_e"].values).groupby(g).mean().sort_index()
        rates[y] = r.values
        flags[y] = bool(np.all(np.diff(r.values) >= -1e-9))
    return flags, rates


def main():
    D = load()
    ref = D.ajccN.values.astype(float)

    rows = []
    group_tables = {}
    for name, fn in SCHEMES.items():
        g = np.asarray(fn(D), dtype=float)
        row = {"scheme": name, "k": int(pd.Series(g).nunique())}
        flags, rates = monotone_flags(D, g)
        for y in Y_LABELS:
            T, E = D[f"{y}_t"].values, D[f"{y}_e"].values.astype(bool)
            row[f"C_{y}"] = round(cindex_by_group(g, T, E), 3)
            try:
                row[f"p_lr_{y}"] = round(multivariate_logrank_test(T, g, E).p_value, 5)
            except Exception:
                row[f"p_lr_{y}"] = np.nan
            row[f"mono_{y}"] = flags[y]
            if name != "AJCC N (기존, 6군)":
                lo, hi = bootstrap_dc(g, ref, T, E)
                row[f"dC_{y}"] = round(cindex_by_group(g, T, E) - cindex_by_group(ref, T, E), 3)
                row[f"dC_{y}_lo"] = round(lo, 3)
                row[f"dC_{y}_hi"] = round(hi, 3)
                row[f"dC_{y}_excl0"] = bool(lo > 0 or hi < 0)
        row["C_mean"] = round(np.mean([row[f"C_{y}"] for y in Y_LABELS]), 4)
        row["mono_all"] = all(flags[y] for y in Y_LABELS)
        rows.append(row)

        agg = D.assign(_g=g).groupby("_g").agg(
            n=("id", "size"),
            PFS_ev=("PFS_e", "sum"), PFS_rate=("PFS_e", "mean"),
            DSS_ev=("DSS_e", "sum"), DSS_rate=("DSS_e", "mean"),
            OS_rate=("OS_e", "mean"), LRRFS_rate=("LRRFS_e", "mean"),
            dep_median=("dep", "median"), cnt_median=("cnt", "median"),
            ENE=("ene", "sum"), contra=("contra", "sum"))
        group_tables[name] = agg.round(3)

    res = pd.DataFrame(rows).sort_values("C_mean", ascending=False)
    res.to_csv(os.path.join(OUT_DIR, "ngroup_proposal_summary.csv"),
               index=False, encoding="utf-8-sig")
    with open(os.path.join(OUT_DIR, "ngroup_proposal_groups.csv"), "w",
              encoding="utf-8-sig") as f:
        for name, agg in group_tables.items():
            f.write(f"# {name}\n")
            agg.to_csv(f)
            f.write("\n")

    lines = []
    lines.append("=== 새 N grouping 후보 검증 (n=133, #52 정정 반영) ===")
    cols = ["scheme", "k", "C_PFS", "C_DSS", "C_OS", "C_LRRFS", "C_mean", "mono_all"]
    lines.append(res[cols].to_string(index=False))
    lines.append("")
    for name in SCHEMES:
        lines.append(f"--- {name} 군별 ---")
        lines.append(group_tables[name].to_string())
        lines.append("")
    txt = "\n".join(lines)
    with open(os.path.join(OUT_DIR, "ngroup_proposal_report.txt"), "w",
              encoding="utf-8") as f:
        f.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
