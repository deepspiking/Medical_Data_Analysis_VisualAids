# -*- coding: utf-8 -*-
"""
실험 4d: 단순 모델 비교 (Logistic / Lasso / Decision Tree) — TabICL 대체
------------------------------------------------------------------
TabICL이 n=133에서 불안정(순위 ρ≈0.5) → 더 단순한 모델들과 비교.

각 모델: 시드 4개 × 5-fold CV (실험 4와 동일 구조)
  - y target: 24개월 사건 여부 (T_HORIZON)
  - fold별: train bootstrap 확장(100명) → fit → test fold score
  - 시드 간 score 순위(Spearman ρ) = rank 안정성
  - 주요 변수·중요도:
      Logistic : |표준화 계수|
      Lasso    : |계수| (0 = 선택 안 됨)
      DecisionTree: feature_importances_

출력:
  exp1_simple_models.csv          : 모델×y×시드 순위 안정성
  exp1_simple_importance.csv      : 모델×y×시드×변수 중요도 (fold 평균)
"""
import os
import time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_CSV = os.path.join(BASE_DIR, "preprocessed_data.csv")
OUT_DIR = os.path.join(BASE_DIR, "results", "exp1")

SEEDS = [42, 123, 2026, 777]
N_FOLDS = 5
T_HORIZON = 24.0
N_EXPAND = 100
Y_LABELS = ["PFS", "DSS", "LRRFS"]

COVS = ["age", "성별", "tumor size (cm)", "DOI (mm)", "differentiation",
        "budding_01vs23", "PNI", "LVI", "RM", "TIL", "TSR",
        "WPOI5_2tier", "HPV/P16_1", "HPV/P16_2", "CCRT_bin"]


def load_prep():
    d = pd.read_csv(DATA_CSV)
    flag = [c for c in d.columns if c.endswith("_known")]
    oh = [c for c in d.columns if c.startswith(("subsite_", "Tx_3tier_", "HPV/P16_"))]
    X = list(dict.fromkeys(COVS + flag + oh))
    return d, X


def fold_split(n, y24, seed, n_folds=N_FOLDS):
    rng = np.random.RandomState(seed)
    pos = np.nonzero(y24 == 1)[0]
    neg = np.nonzero(y24 == 0)[0]
    rng.shuffle(pos)
    rng.shuffle(neg)
    pos_folds = np.array_split(pos, n_folds)
    neg_folds = np.array_split(neg, n_folds)
    folds = []
    for f in range(n_folds):
        test_idx = np.concatenate([pos_folds[f], neg_folds[f]])
        train_idx = np.concatenate([
            np.concatenate([pos_folds[j] for j in range(n_folds) if j != f]),
            np.concatenate([neg_folds[j] for j in range(n_folds) if j != f])])
        folds.append((train_idx, test_idx))
    return folds


def fit_model(name, X_train, y_train, X_test, seed):
    """모델 fit → (test score, 변수 중요도 dict). sklearn은 빠르고 안정적."""
    from sklearn.linear_model import LogisticRegression, Lasso
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.preprocessing import StandardScaler

    sc = StandardScaler()
    Xtr_s = sc.fit_transform(X_train)
    Xte_s = sc.transform(X_test)

    if name == "Logistic":
        m = LogisticRegression(max_iter=2000, C=1.0, random_state=seed)
        m.fit(Xtr_s, y_train)
        score = m.predict_proba(Xte_s)[:, 1]
        imp = {i: abs(c) for i, c in enumerate(m.coef_[0])}
    elif name == "Lasso":
        m = Lasso(alpha=0.01, max_iter=5000, random_state=seed)
        m.fit(Xtr_s, y_train)
        score = m.predict(Xte_s)
        imp = {i: abs(c) for i, c in enumerate(m.coef_)}
    elif name == "DecisionTree":
        m = DecisionTreeClassifier(max_depth=4, min_samples_leaf=5,
                                   random_state=seed)
        m.fit(Xtr_s, y_train)
        score = m.predict_proba(Xte_s)[:, 1]
        imp = {i: v for i, v in enumerate(m.feature_importances_)}
    else:
        raise ValueError(name)
    return score, imp


def main():
    t0 = time.time()
    d, X = load_prep()
    print(f"단순 모델 비교 | X {len(X)}개, 시드 {len(SEEDS)} × fold {N_FOLDS} × y {len(Y_LABELS)}",
          flush=True)

    MODELS = ["Logistic", "Lasso", "DecisionTree"]
    ckpt_csv = os.path.join(OUT_DIR, "exp1_simple_models.csv")
    imp_csv = os.path.join(OUT_DIR, "exp1_simple_importance.csv")

    done = set()
    if os.path.exists(ckpt_csv):
        prev = pd.read_csv(ckpt_csv)
        done = set(zip(prev["model"], prev["y"], prev["seed"]))
        print(f"체크포인트: {len(done)}개 시드 완료", flush=True)

    rows = []
    imp_rows = []
    if os.path.exists(ckpt_csv):
        rows = pd.read_csv(ckpt_csv).to_dict("records")
    if os.path.exists(imp_csv):
        imp_rows = pd.read_csv(imp_csv).to_dict("records")

    for y in Y_LABELS:
        y_time, y_event = f"{y}_time", f"{y}_event"
        dfx = d[X + [y_time, y_event]].dropna().reset_index(drop=True)
        dfx["y24"] = ((dfx[y_event] == 1) & (dfx[y_time] <= T_HORIZON)).astype(int)
        X_all = dfx[X].values.astype(float)
        y_all = dfx["y24"].values.astype(int)

        for model in MODELS:
            seed_scores = {}
            for seed in SEEDS:
                key = (model, y, seed)
                if key in done:
                    continue
                folds = fold_split(len(dfx), y_all, seed)
                all_scores = np.full(len(dfx), np.nan)
                imp_acc = {}
                for f, (train_idx, test_idx) in enumerate(folds):
                    rng_try = np.random.RandomState(seed * 1000 + f * 10)
                    boot_idx = rng_try.choice(train_idx, size=N_EXPAND, replace=True)
                    score, imp = fit_model(model, X_all[boot_idx], y_all[boot_idx],
                                           X_all[test_idx], seed * 100 + f * 10)
                    all_scores[test_idx] = score
                    for v_idx, v in enumerate(X):
                        imp_acc[v] = imp_acc.get(v, 0.0) + imp.get(v_idx, 0.0)
                if np.isnan(all_scores).sum() > len(dfx) * 0.2:
                    print(f"  [{y}] {model} seed={seed} score 누락", flush=True)
                    continue
                seed_scores[seed] = all_scores
                # 주요 변수 (fold 평균 중요도)
                n_done = sum(1 for _ in folds)
                imp_mean = {v: val / n_done for v, val in imp_acc.items()}
                for v in X:
                    imp_rows.append({"model": model, "y": y, "seed": seed,
                                     "var": v, "importance": imp_mean.get(v, 0.0)})
                print(f"  [{y}] {model} seed={seed} 완료 ({time.time()-t0:.0f}s)",
                      flush=True)

            if len(seed_scores) == len(SEEDS):
                mats = np.column_stack([seed_scores[s] for s in SEEDS])
                valid = ~np.isnan(mats).any(axis=1)
                M = mats[valid]
                rho_pairs = []
                for a in range(len(SEEDS)):
                    for b in range(a + 1, len(SEEDS)):
                        rho_pairs.append(float(spearmanr(M[:, a], M[:, b])[0]))
                rows.append({"model": model, "y": y,
                             "n": int(valid.sum()),
                             "seed_rho_mean": round(float(np.mean(rho_pairs)), 3),
                             "seed_rho_min": round(float(np.min(rho_pairs)), 3),
                             "seed_rho_max": round(float(np.max(rho_pairs)), 3)})
                pd.DataFrame(rows).to_csv(ckpt_csv, index=False, encoding="utf-8-sig")
                pd.DataFrame(imp_rows).to_csv(imp_csv, index=False,
                                              encoding="utf-8-sig")

    res = pd.DataFrame(rows)
    res.to_csv(ckpt_csv, index=False, encoding="utf-8-sig")
    pd.DataFrame(imp_rows).to_csv(imp_csv, index=False, encoding="utf-8-sig")

    # 요약 출력
    print("\n=== 모델별 score 순위 안정성 (시드 간 Spearman ρ) ===")
    print(res.to_string(index=False))
    print("\n=== 모델별 주요 변수 (중요도 상위 5) ===")
    for model in MODELS:
        for y in Y_LABELS:
            sub = pd.DataFrame(imp_rows)
            sub = sub[(sub["model"] == model) & (sub["y"] == y)]
            if sub.empty:
                continue
            vm = sub.groupby("var")["importance"].mean().sort_values(ascending=False)
            print(f"  {model} | {y}: " + ", ".join(
                f"{v}({val:.3f})" for v, val in vm.head(5).items()))

    print(f"\n[완료] {ckpt_csv}, {imp_csv} (총 {time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()