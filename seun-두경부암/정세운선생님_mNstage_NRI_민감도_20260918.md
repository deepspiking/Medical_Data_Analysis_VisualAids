# mNstage LRRFS NRI 신호 — 민감도·다중비교 분석 (2026-09-18)

> 질문: mNstage의 LRRFS NRI(0.497)가 **다중비교를 견디는 진짜 신호인가, 우연인가?**
> 코드: `nri_sensitivity.py` · 결과: `results/exp1/nri_sensitivity/`

## 결론 (먼저)

> **취약한 경계선 신호(borderline)로, "mNstage가 기존 N보다 우월하다"는 근거로는 부족합니다.**
> - Holm 다중비교 보정 후 LRRFS는 **간신히 유의(p_holm=0.047)**, permutation test는 **p=0.054(비유의)**.
> - 유의성은 **24~36개월에서만** 나타나고(12·48개월은 비유의), 개선은 **non-event 재분류에서만** 발생.
> - C-index는 개선되지 않음(ΔC CI 0 포함).
> - 반면 **mTstage는 4개 endpoint 전부 Holm 보정 후에도 유의** → 우월성의 실제 근거는 T축.

---

## 1. 12개 비교 NRI + Holm-Bonferroni (t=24, fixed-score bootstrap B=3,000)

| pair | y | NRI | 95% CI | p | p_Holm | Holm 유의 |
|---|---|---|---|---|---|---|
| mStage vs ajcc8th | PFS | 0.848 | [0.508, 1.152] | <0.001 | <0.001 | ✅ |
| mStage vs ajcc8th | LRRFS | 1.007 | [0.647, 1.334] | <0.001 | <0.001 | ✅ |
| mStage vs ajcc8th | DSS / OS | 0.250 / 0.135 | 0 포함 | 0.28 / 0.54 | ns | ❌ |
| **mTstage vs T** | PFS / DSS / OS / LRRFS | 0.717 / 1.061 / 0.883 / 0.631 | 모두 0 제외 | <0.001 | **<0.001** | ✅ **4/4** |
| mNstage vs N | PFS | 0.592 | [0.243, 0.925] | 0.001 | 0.004 | ✅* |
| mNstage vs N | OS | 0.503 | [0.094, 0.900] | 0.015 | 0.061 | ❌ |
| **mNstage vs N** | **LRRFS** | **0.497** | **[0.126, 0.867]** | **0.009** | **0.047** | **✅ (경계)** |
| mNstage vs N | DSS | 0.442 | [−0.007, 0.887] | 0.056 | 0.168 | ❌ |

\* 주의: mNstage PFS의 유의는 **fixed-score bootstrap**에서만. 모델 재적합을 반영한 **refit bootstrap에서는 CI가 0을 포함**(0.592 [−0.307, 0.957], 공식 실험2) → **robust하지 않음**.

## 2. mNstage: refit bootstrap(모델 불확실성 반영) 기준

| y | NRI | 95% CI (refit, 실험2 B=1,000) | 0 제외 |
|---|---|---|---|
| PFS | 0.592 | [−0.307, 0.957] | ❌ |
| DSS | 0.442 | [−0.566, 0.925] | ❌ |
| OS | 0.503 | [−0.772, 0.894] | ❌ |
| **LRRFS** | **0.497** | **[0.070, 0.919]** | **✅ (유일)** |

→ **mNstage에서 refit bootstrap을 견디는 것은 LRRFS 하나뿐**.

## 3. 시간대별 (mNstage LRRFS)

| t (개월) | NRI | 95% CI | p |
|---|---|---|---|
| 12 | 0.339 | [−0.044, 0.722] | 0.076 |
| **24** | **0.497** | [0.126, 0.867] | **0.009** |
| **36** | **0.429** | [0.037, 0.796] | **0.031** |
| 48 | 0.358 | [−0.072, 0.765] | 0.107 |

→ **24·36개월에서만 유의**, 12·48개월은 비유의(일관성 부족).

## 4. Permutation test

- mNstage 위험점수를 무작위 재배정(5,000회) → **permutation p = 0.054** (비유의, 경계).

## 5. 왜 취약한가 — 분해(decomposition)

LRRFS mNstage (t=24): **NRI = NRI_event − NRI_nonevent = (−0.128) − (−0.625)**
- event 재분류는 오히려 **약간 악화**(−0.128), 개선은 **non-event(무재발자)를 저위험으로 재분류**한 데서 발생.
- 임상적으로는 "사건을 더 잘 잡는다"기보다 "사건 없는 환자를 덜 치료"하는 쪽 → 근거 약함.

## 6. 종합 판정

| 주장 | 근거 강도 |
|---|---|
| **mTstage > 기존 T** | **강함** (4/4 endpoint, Holm 유의, refit 견딤) |
| mStage > ajcc8th | 중간 (PFS·LRRFS Holm 유의) |
| **mNstage > 기존 N** | **약함/취약** — LRRFS NRI만 경계선, C-index는 개선 없음 |

**권고**: 논문에서는 **"mNstage는 AJCC N과 구동등(ΔC≈0)"**으로 본문 기술하고,
LRRFS NRI는 **exploratory/secondary**로만 언급(다중비교·경계성·non-event driven 명시).
수정병기의 **주된 기여는 mTstage(PD cellularity)**로 정리하는 것이 리뷰어 방어에 안전합니다.
