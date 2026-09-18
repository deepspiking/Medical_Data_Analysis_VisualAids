# Nomogram 해석 가이드 — 개념과 읽는 법 (2026-09-18)

> 정세운 선생님 질문(개념 중복 / 변수 중요도 해석)에 대한 정리.
> 관련 코드: `nomogram_analysis.py` · 그림: `results/exp1/nomogram/`

---

## 1. Nomogram 읽는 순서 (환자 1명 계산법)

1. 환자의 각 변수 값을 **각 축에서 찾아** → 바로 위 **Points 자(ruler)** 에서 점수를 읽음
2. 모든 변수 점수를 **합산** → **Total Points** 축에서 위치 확인
3. Total Points에서 **세로로 내려가** 3년/5년 생존확률 축을 읽음
4. (보조) **Risk group 컬러 밴드**(Low 초록 / Medium 주황 / High 빨강)로 위험군 즉시 확인

예) OS nomogram에서 mTstage=4, mNstage=3이면 → 총점이 High-risk(빨강) 구간에 들어가 3·5년 생존확률이 낮게 읽힘.

---

## 2. 변수 중요도 = "축 길이(점수 범위)" — HR이 아님

- nomogram에서 한 변수의 영향력은 **그 변수 축의 길이(0~최대 점수)** 로 표현됩니다.
- **축 길이 = |회귀계수 β| × (그 변수의 관측 값 범위)**.
- 따라서 **연속형 변수는 범위가 넓어 축이 길어집니다** — 단위당 HR이 작아도.

**primary(주요변수 5개) nomogram의 축 길이 (각 endpoint 최대=100으로 정규화)**

| endpoint | 1위 | 2위 | 3위 | 4위 | 5위 |
|---|---|---|---|---|---|
| DSS | Size 100 | Deposit 73 | DOI 58 | PD 50 | Bilateral 25 |
| OS | Size 100 | DOI 85 | Deposit 73 | PD 55 | Bilateral 45 |
| PFS | DOI 100 | Deposit 89 | PD 65 | Size 61 | Bilateral 32 |
| LRRFS | Deposit 100 | DOI 98 | PD 55 | Size 47 | Bilateral 9 |

> **주의(중요)**: 축 길이는 `|β| × 값의 "범위"`라서 **이상치/왜도·단위에 민감**합니다.
> 예: `Size`·`DOI`는 연속형이라 범위가 넓어 축이 길게 나옵니다.
> 반면 **표준화(|β|×SD) 기준으로는 `PD`가 4개 endpoint 모두 1위**입니다
> (`정세운선생님_구성요소_중요도_20260918.md`).
> → **"nomogram 축 길이"와 "통계적 중요도(|β|×SD)"가 다를 수 있음**을 함께 해석해야 합니다.

---

## 3. "stage와 그 구성요소를 같이 넣는 것" — 개념 정리

**선생님 지적이 맞습니다: 중복(double-counting)이라 개념적으로 어색하며, 주의가 필요합니다.**

| 항목 | 설명 |
|---|---|
| 왜 문제인가 | mTstage는 size·DOI·bone invasion·PD의 **함수** → 같은 정보가 두 번 들어감 |
| 증상 | size/DOI가 긴 축을 차지해 **"가장 중요한 변수"를 잘못 읽게** 만듦; 계수 해석 왜곡; 공선성 |
| 예시 nomogram은? | T stage + tumor size를 같이 넣는 경우가 많음 = **관행**이지 "옳아서"가 아님 |
| 언제 허용? | 순수 예측 성능 목적 + 교차검증에서 이득이 확인될 때 |
| 언제 피하나? | **병기 자체를 보여주는 nomogram** → composite만 (우리 선택) |

**최종 결정(정세운 선생님, 2026-09-18)**: 핵심은 "**PD component 개념으로 stage가 바뀌었다**"이고,
**T stage의 기능은 앞선 그래프·수치에서 이미 검증**되었습니다. 따라서 T stage와 구성요소가 겹치는
상황에서는 **stage를 빼고 개별 독립변수(구성요소)를 넣은 nomogram을 primary**로 사용합니다
(개별 변수 중요도 강조). **복합 병기(mTstage·mNstage) nomogram은 supplement**(`nomogram_*_stage.png`).
- 두 버전 성능은 동등 수준(구성요소 버전 apparent C가 근소 우세: DSS 0.878 vs 0.859 등).

---

## 4. 자주 하는 오해

| 오해 | 실제 |
|---|---|
| "HR이 큰 변수가 nomogram에서 가장 중요" | ❌ **축 길이**가 중요도. HR은 단위당 위험이고, 축 길이는 HR×값 범위 |
| "stage와 component를 같이 넣으면 틀린 분석" | ❌ 통계적으로 금지는 아님. 단 **중복·해석 주의** |
| "축 길이가 길면 임상적으로 꼭 더 중요" | ⚠️ 영향력 크기일 뿐, 인과/중요도와 동일하지 않음 |
| "mNstage가 nomogram에서 유의하니 기존 N보다 우월" | ❌ 유의(예후인자)와 **기존 N 대비 우월(ΔC≈0)** 은 다른 얘기 |

---

## 5. 4개 endpoint nomogram 성능 요약 (주요변수 모델, n=133)

| endpoint | 모델 | C-index (corrected) | Calibration slope (3·5년) | Time-AUC (3·5년) | Risk-group log-rank |
|---|---|---|---|---|---|
| **DSS** | 주요변수 5개 | **0.869 (0.856)** [0.806–0.929] | 1.31 / 1.38 | 0.869 / 0.870 | p<0.001 |
| **OS** | 주요변수 5개 | 0.805 (0.795) [0.726–0.877] | 0.96 / 1.02 | 0.793 / 0.788 | p<0.001 |
| **PFS** | 주요변수 5개 | 0.726 (0.714) [0.653–0.792] | 1.08 / 0.92 | 0.737 / 0.727 | p<0.001 |
| **LRRFS** | 주요변수 5개 | 0.718 (0.694) [0.627–0.790] | 1.61 / 1.24 | 0.767 / 0.737 | p<0.001 |

> 주요변수 = PD_01vs2 · LN deposit · DOI · tumor size · bilateral.
> Companion(복합 병기 mTstage·mNstage, 본문 동반): DSS 0.859 / OS 0.803 / PFS 0.718 / LRRFS 0.697.
> **DSS가 가장 강함, LRRFS가 가장 약함** — endpoint 난이도 차이.

---

## 6. "그럼 stage를 빼고 구성요소만 넣으면?" — 개별 인자 중요도

병기 nomogram에서는 중복 때문에 구성요소를 빼지만, **stage를 전부 빼고 구성요소만 넣은 모델**로
"어떤 인자가 가장 중요한가"를 별도로 확인할 수 있습니다(`정세운선생님_구성요소_중요도_20260918.md`).

**결과(구성요소만, 표준화 |β|×SD): PD_01vs2가 4개 endpoint 중 3개에서 1위, LRRFS에서 2위.**
임상인자까지 넣어도 **PD가 4/4에서 1위**. 그 뒤로 LN deposit · DOI · tumor size.

![구성요소 중요도](results/exp1/component_importance/component_importance.png)

→ **PD_01vs2가 가장 강력한 구성요소** = 수정 T staging(PD-high → mT4)의 근거를 정면 지지.
→ 단, 구성요소 간 상관(size↔DOI, count↔deposit)으로 개별 p값은 불안정 → **순위는 모델 내 상대 기여도로 해석**.
