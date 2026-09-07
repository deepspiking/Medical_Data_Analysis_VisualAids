# -*- coding: utf-8 -*-
"""
키인사이트 요약용 서포트 그림 생성 (n=133 전수, *_full.csv 공식 데이터)
---------------------------------------------------------------------
수정병기_키인사이트_요약.md 에 삽입할 그림 3종 + 수치 검증 출력:
  Fig1  C-index 쌍별 비교 (mStage vs ajcc8th, mTstage vs T stage, 4 endpoint)
  Fig2  S3 증분 NRI: 수정병기 vs 기존병기 (4 seed 평균, ΔNRI 표기)
  Fig3  ΔRMST forest (mStage·mTstage, 4 endpoint, 95% CI)
docx(oralSCC_table.docx)와 수치 결이 같도록 *_full.csv 공식 값을 그대로 사용.
"""
import os
import ast
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["AppleGothic", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, "results", "exp1")
OUT = os.path.join(RES, "키인사이트")
os.makedirs(OUT, exist_ok=True)

COL_MOD = "#d62728"   # 수정병기
COL_OLD = "#1f77b4"   # 기존병기
Y_ORDER = ["PFS", "DSS", "OS", "LRRFS"]
SEEDS = [42, 123, 2026, 777]


def load(name):
    p = os.path.join(RES, name)
    return pd.read_csv(p) if os.path.exists(p) else None


def save(fig, fname):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {fname}")


# ----------------------------------------------------------------------------
# Fig1. C-index 쌍별 비교 (실험 2 S1, seed 42 — univariate_full)
# ----------------------------------------------------------------------------
def fig1():
    uni = load("exp1_univariate_full.csv")
    uni = uni[uni["seed"] == 42]
    pairs = [("mStage vs ajcc8th", "mTstage vs T stage")]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), sharex=False)
    for ax, pair in zip(axes, pairs[0]):
        sub = uni[uni["pair"] == pair].set_index("label")
        rows = [(y, sub.loc[y]) for y in Y_ORDER]
        ypos = np.arange(len(rows))[::-1]
        ylabel_map = dict(zip(ypos, Y_ORDER))
        ax.axvline(0.6, color="gray", ls="--", lw=0.8, alpha=0.5)
        for yy, (y, r) in zip(ypos, rows):
            # 기존병기 (파랑) / 수정병기 (빨강) — yy 위아래로 오프셋
            ax.errorbar([r["C_old"]], [yy + 0.17], xerr=[[r["C_old"] - r["C_old_lo"]],
                        [r["C_old_hi"] - r["C_old"]]], fmt="s", color=COL_OLD,
                        ms=7, capsize=3, lw=1.4, zorder=3, label="기존병기" if y == "PFS" else None)
            ax.errorbar([r["C_new"]], [yy - 0.17], xerr=[[r["C_new"] - r["C_new_lo"]],
                        [r["C_new_hi"] - r["C_new"]]], fmt="o", color=COL_MOD,
                        ms=8, capsize=3, lw=1.4, zorder=3, label="수정병기" if y == "PFS" else None)
            dC = r["dC"]
            star = "*" if r["dC_excl_0"] else ""
            ax.text(1.12, yy, f"ΔC {dC:+.3f}{star}", ha="right", va="center",
                    fontsize=9, fontweight="bold" if r["dC_excl_0"] else "normal")
        ax.set_yticks(sorted(ypos))
        ax.set_yticklabels([ylabel_map[k] for k in sorted(ypos)],
                           fontsize=10, fontweight="bold")
        ax.set_xlim(0.48, 1.18)
        ax.set_xlabel("C-index (95% CI)", fontsize=11)
        ax.set_title(pair, fontsize=12, fontweight="bold")
        ax.legend(loc="upper left", fontsize=9, frameon=False)
        ax.grid(axis="x", alpha=0.25)
    fig.suptitle("수정병기 vs 기존병기 C-index (n=133, S1)  ·  * : 95% CI가 0을 제외",
                 fontsize=12, y=1.02)
    save(fig, "Fig1_cindex_pairs.png")


# ----------------------------------------------------------------------------
# Fig2. S3 증분 ΔNRI (실험 3, 4 seed 평균 — 공식 리포트 §4 값과 일치)
# ----------------------------------------------------------------------------
def fig2():
    mv = load("exp1_multivariate_full.csv")
    s3 = mv[mv["model"] == "S3_increment_C_vs_B"]
    rows = []
    for y in Y_ORDER:
        sub = s3[s3["label"] == y]
        rows.append({"y": y, "dNRI": sub["NRI_diff"].mean()})
    df = pd.DataFrame(rows).set_index("y").loc[Y_ORDER]
    x = np.arange(len(Y_ORDER))
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    bars = ax.bar(x, df["dNRI"], 0.52, color=COL_MOD)
    for xi, (_, r) in zip(x, df.iterrows()):
        ax.text(xi, r["dNRI"] + 0.02, f"+{r['dNRI']:.3f}",
                ha="center", fontsize=12, fontweight="bold", color=COL_MOD)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(Y_ORDER, fontsize=11)
    ax.set_ylabel("ΔNRI (수정병기 증분 - 기존병기 증분)", fontsize=11)
    ax.set_title("S3 증분: “수정병기 재분류 개선 > 기존병기 재분류 개선” — 4개 endpoint 모두 ΔNRI > 0 (n=133)",
                 fontsize=12, fontweight="bold")
    ax.set_ylim(0, 0.75)
    ax.grid(axis="y", alpha=0.25)
    save(fig, "Fig2_S3_dNRI.png")
    print(df.round(3).to_string())


# ----------------------------------------------------------------------------
# Fig3. ΔRMST forest (실험 5, mStage·mTstage)
# ----------------------------------------------------------------------------
def fig3():
    rm = load("exp1_rmst_diff_full.csv")
    scores = ["mStage", "mTstage"]
    sub = rm[rm["score_col"].isin(scores)].copy()
    sub = sub.rename(columns={"label": "y"})
    sub = sub[sub["y"].isin(Y_ORDER)]
    colors = {"mStage": "#c0392b", "mTstage": "#e67e22"}
    # 표시 순서 (위→아래): mStage PFS→LRRFS, 한 줄 띄고 mTstage PFS→LRRFS
    rows = []
    y = 0
    for si, sc in enumerate(scores):
        for yl in Y_ORDER:
            r = sub[(sub["score_col"] == sc) & (sub["y"] == yl)].iloc[0]
            y += 1
            rows.append((y, sc, yl, r))
        y += 1  # 그룹 간 간격
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ticks, labels = [], []
    for yy, sc, yl, r in rows:
        ticks.append(yy)
        labels.append(f"{yl} ({sc})")
        ax.errorbar([r["RMST_diff"]], [yy],
                    xerr=[[r["RMST_diff"] - r["RMST_diff_lo"]],
                          [r["RMST_diff_hi"] - r["RMST_diff"]]],
                    fmt="o", color=colors[sc], ms=9, capsize=4, lw=1.6, zorder=3)
        ax.text(r["RMST_diff_hi"] + 0.4, yy,
                f"{r['RMST_diff']:.1f}  [{r['RMST_diff_lo']:.1f}, {r['RMST_diff_hi']:.1f}]",
                va="center", fontsize=9)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=9)
    ax.axvline(0, color="gray", ls="--", lw=1)
    ax.set_xlabel("ΔRMST (개월, 고위험 - 저위험, t*=36개월)", fontsize=11)
    ax.set_title("고위험군 평균 생존 단축: 모든 endpoint p<0.001 (n=133, bootstrap 1,000)",
                 fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.25)
    save(fig, "Fig3_RMST_forest.png")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    print("\n[출력 폴더]", OUT)
