# -*- coding: utf-8 -*-
"""
실험 4 보조 단계: SHAP(4-7)·RSF(4-8) 단독 실행
------------------------------------------------------------------
실험 4의 핵심(CV 안정성 + score↔생존)은 이미 완료·저장됨.
이 스크립트는 보조 단계만 독립 실행 — 부모 segfault와 무관하게
별도 프로세스에서 수행하고 결과만 저장한다.

  SHAP: TabICL fit + shap.Explainer → 변수 중요도 vs Cox HR 순위
  RSF : Random Survival Forest OOB C-index
"""
import os
import sys
import json
import time
import tempfile
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_CSV = os.path.join(BASE_DIR, "preprocessed_data.csv")
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")
os.makedirs(OUT_DIR, exist_ok=True)

SEED = 42
T_HORIZON = 24.0
Y_LABELS = ["PFS", "DSS", "LRRFS"]
STAGE_VARS = ["ajcc8th_STAGE", "mStage"]

COVS = ["age", "성별", "tumor size (cm)", "DOI (mm)", "differentiation",
        "budding_01vs23", "PNI", "LVI", "RM", "TIL", "TSR",
        "WPOI5_2tier", "HPV/P16_1", "HPV/P16_2", "CCRT_bin"]


def load_prep():
    d = pd.read_csv(DATA_CSV)
    flag_cols = [c for c in d.columns if c.endswith("_known")]
    oh_cols = [c for c in d.columns if c.startswith(("subsite_", "Tx_3tier_", "HPV/P16_"))]
    return d, COVS + flag_cols + oh_cols


def shap_rank_subprocess(d, X, y_time, y_event):
    """SHAP 변수 중요도 vs Cox HR 순위 — subprocess 격리 실행."""
    dfx = d[X + [y_time, y_event]].dropna().reset_index(drop=True)
    dfx["y24"] = ((dfx[y_event] == 1) & (dfx[y_time] <= T_HORIZON)).astype(int)
    Xm = dfx[X].values.astype(float)
    ym = dfx["y24"].values.astype(int)

    tmpdir = tempfile.mkdtemp(prefix="tabicl_shap2_")
    x_path = os.path.join(tmpdir, "X.npy")
    y_path = os.path.join(tmpdir, "y.npy")
    out_path = os.path.join(tmpdir, "shap.json")
    np.save(x_path, Xm)
    np.save(y_path, ym)

    script = f'''
import os, sys, json, numpy as np
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, {BASE_DIR!r})
from tabicl import TabICLClassifier
import shap
Xm = np.load({x_path!r})
ym = np.load({y_path!r})
clf = TabICLClassifier(n_estimators=1, random_state={SEED})
clf.fit(Xm, ym)
def pred_fn(x):
    return clf.predict(x)
explainer = shap.Explainer(pred_fn, Xm, feature_names={X!r})
shap_vals = explainer(Xm, max_evals=1000)
shap_map = {{v: float(np.abs(shap_vals.values).mean(axis=0)[i]) for i, v in enumerate({X!r})}}
json.dump(shap_map, open({out_path!r}, "w"))
os._exit(0)
'''
    env = dict(os.environ)
    env["KMP_DUPLICATE_LIB_OK"] = "TRUE"
    env["OMP_NUM_THREADS"] = "1"
    try:
        p = subprocess.Popen([sys.executable, "-c", script],
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                             start_new_session=True, env=env)
        try:
            p.communicate(timeout=1800)
        except subprocess.TimeoutExpired:
            p.kill()
            p.communicate()
        if not os.path.exists(out_path):
            return None
        shap_map = json.load(open(out_path))
    except Exception:
        return None
    finally:
        for f in (x_path, y_path, out_path):
            try:
                os.remove(f)
            except OSError:
                pass

    # 단변량 Cox HR (부모에서 — lifelines만 사용, segfault 위험 없음)
    from lifelines import CoxPHFitter
    hr_map = {}
    for v in X:
        sub = dfx[[v, y_time, y_event]].dropna()
        if sub[y_event].sum() < 5 or sub[v].nunique() < 2:
            continue
        try:
            cph = CoxPHFitter(penalizer=0.0)
            cph.fit(sub, duration_col=y_time, event_col=y_event)
            hr_map[v] = float(np.exp(cph.params_[v]))
        except Exception:
            continue
    common = [v for v in X if v in shap_map and v in hr_map]
    if len(common) < 5:
        return None
    sh = np.array([shap_map[v] for v in common])
    hr = np.array([np.abs(np.log(hr_map[v])) for v in common])
    corr = float(np.corrcoef(sh, hr)[0, 1])
    rho = float(spearmanr(sh, hr)[0])
    return {"n": len(common), "pearson": corr, "spearman": rho}


def rsf_cindex(d, X, y_time, y_event):
    """RSF OOB C-index (sksurv — 비교적 안정)."""
    from sksurv.ensemble import RandomSurvivalForest
    dfx = d[X + [y_time, y_event]].dropna()
    Xm = dfx[X].values.astype(float)
    struct = np.array([(bool(e), t) for t, e in
                       zip(dfx[y_event].values, dfx[y_time].values.astype(float))],
                      dtype=[("event", bool), ("time", float)])
    rsf = RandomSurvivalForest(n_estimators=500, random_state=SEED, n_jobs=-1)
    rsf.fit(Xm, struct)
    return float(rsf.score(Xm, struct))


def main():
    t0 = time.time()
    d, X = load_prep()
    print(f"실험 4 보조 | X {len(X)}개", flush=True)

    # SHAP (4-7)
    shap_rows = []
    for y in Y_LABELS:
        y_time, y_event = f"{y}_time", f"{y}_event"
        try:
            r = shap_rank_subprocess(d, X, y_time, y_event)
        except Exception:
            r = None
        if r is not None:
            shap_rows.append({"y": y, "n": r["n"], "pearson": r["pearson"],
                              "spearman": r["spearman"]})
            print(f"  SHAP vs Cox|HR| [{y}]: ρ={r['spearman']:+.3f} (n={r['n']})",
                  flush=True)
        else:
            print(f"  SHAP [{y}]: 실패 (segfault/시간초과)", flush=True)
    if shap_rows:
        pd.DataFrame(shap_rows).to_csv(os.path.join(OUT_DIR, "exp1_tabicl_shap.csv"),
                                       index=False, encoding="utf-8-sig")

    # RSF (4-8)
    rsf_rows = []
    for y in Y_LABELS:
        y_time, y_event = f"{y}_time", f"{y}_event"
        try:
            c = rsf_cindex(d, X, y_time, y_event)
        except Exception:
            c = np.nan
        rsf_rows.append({"y": y, "model": "RSF_full", "OOB_C": c})
        for s in STAGE_VARS:
            try:
                c2 = rsf_cindex(d, X + [s], y_time, y_event)
            except Exception:
                c2 = np.nan
            rsf_rows.append({"y": y, "model": f"RSF_{s}", "OOB_C": c2})
        print(f"  RSF [{y}] 완료", flush=True)
    pd.DataFrame(rsf_rows).to_csv(os.path.join(OUT_DIR, "exp1_rsf_cindex.csv"),
                                  index=False, encoding="utf-8-sig")

    print(f"\n[완료] exp1_tabicl_shap.csv, exp1_rsf_cindex.csv (총 {time.time()-t0:.0f}s)",
          flush=True)


if __name__ == "__main__":
    main()