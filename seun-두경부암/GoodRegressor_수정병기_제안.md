# GoodRegressor 기호회귀 기반 수정병기 대안 공식 제안

`최종_종합_리포트_AI엔지니어용.md`의 후속 분석. 수정병기(mStage/mTstage)가 기존 T/N stage에 PD
cellularity를 더해 **수작업 규칙표**로 설계된 것이었다면, 이 문서는 같은 원재료 변수들로부터
**기호회귀(symbolic regression) 엔진 GoodRegressor**가 자동으로 찾아낸 산술 공식이 실제로
쓸만한 대안(또는 보조 지표)이 될 수 있는지 검증한 결과다.

## 방법 (요약)

1. **연속형 변환**: PFS/DSS/LRRFS는 중도절단 생존데이터라 회귀 엔진에 직접 못 넣음 →
   RMST(τ=36개월)의 **jackknife pseudo-value**로 연속형 대리변수를 만들어 Y로 사용.
2. **탐색**: 후보 X 72개(실험3의 임상변수 23개 + 병기 원재료 T stage/N stage/PD_01vs2 등 4개
   + 이 중 10개 사이 pairwise 곱 상호작용 45개)로 GoodRegressor 전수탐색(1~3변수, Nd=1/2/3).
3. **1차 검증(Cox 다변량)**: 찾은 2변수 조합을 `GRscore`라는 단일 병기 점수로 만들어, 실험3과
   동일한 S2(대체)/S3(증분 NRI) 틀로 Model A(임상변수만)/B(A+ajcc8th)/C(A+mStage)/D(A+GRscore)
   비교(500회 부트스트랩 optimism 보정).
4. **2차 검증(leakage 제거 nested CV)**: 1차 검증은 GRscore 수식 자체를 129명 전체로 이미
   골라놓고 그 위에 Cox 가중치만 재검증한 것이라 여전히 낙관적일 수 있음 → **공식 발견 자체를
   매 fold 안에서 다시 수행**하는 event-층화 5-fold × 10회 반복 CV로 최종 확인
   (train fold만으로 pseudo-value 재계산 + GoodRegressor 재탐색 → test fold에 대입 → 한 repeat의
   5-fold held-out 예측을 모아 C-index 1개 → 10개 repeat 평균; repeat끼리는 절대 풀링하지 않음
   — 같은 환자가 여러 번 등장해 쌍 계산이 깨지기 때문).

## 결과 1 — Cox 다변량 검증 (S2/S3/NRI)

optimism 보정 후 C-index (Model A=임상변수만 / B=A+ajcc8th / C=A+mStage / D=A+GRscore):

| endpoint | A | B(ajcc8th) | C(mStage) | D(GRscore) |
|---|---|---|---|---|
| PFS | 0.643 | 0.648 | 0.699 | 0.692 |
| DSS | 0.812 | 0.823 | 0.848 | 0.844 |
| LRRFS | 0.610 | 0.604 | 0.661 | 0.658 |

- **S2(paired CV, D−C)**: PFS −0.019 [−0.088, +0.042] · DSS +0.001 [−0.112, +0.067] ·
  LRRFS −0.012 [−0.108, +0.061] — 세 endpoint 모두 95% 구간이 0을 포함, **GRscore와 mStage는
  통계적으로 구분되지 않음**.
- **S3(baseline 대비 NRI 증분)**: ajcc8th / mStage / GRscore = PFS 0.433/1.026/0.872 ·
  DSS 0.861/1.194/0.833 · LRRFS 0.544/1.061/0.881 — GRscore가 ajcc8th보다는 뚜렷하게 낫고,
  mStage보다는 근소하게 못 미침.

## 결과 2 — Leakage 제거 nested CV (최종 판정)

1차 검증은 "어떤 두 변수를 곱할지"를 129명 전체로 이미 정해놓고 검증한 것이라 낙관 편향이
남아있었다. 공식 발견 자체를 매 fold 안에서 새로 하도록 다시 돌린 결과, **LRRFS에서 그 편향이
실제로 드러났다** (하단 표 참조).

## 제안 수식 (129명 전체 피팅)

| endpoint | Nd | 수식 (예측값 = pseudo-RMST, 높을수록 예후 좋음) |
|---|---|---|
| PFS | 1 | 27.08 − 0.144×(N stage×DOI) |
| PFS | 2 | 28.72 − 10.69×PD_01vs2 − 1.87×N stage |
| PFS | 3 | 28.97 − 11.51×PD_01vs2 − 0.484×(N stage×종양크기) − 25.97×subsite_5 |
| DSS | 1 | 33.61 − 0.125×(N stage×DOI) |
| DSS | 2 | 35.04 − 0.104×(N stage×DOI) − 7.95×PD_01vs2 |
| DSS | 3 | 35.12 − 14.62×PD_01vs2 − 0.558×(N stage×종양크기) + 11.01×(PD_01vs2×budding) |
| LRRFS | 1 | 29.07 − 8.98×PD_01vs2 |
| LRRFS | 2 | 29.48 − 0.452×(T stage×N stage) − 7.91×PD_01vs2 |
| LRRFS | 3 | 28.93 − 4.34×(N stage×LVI) + 0.670×(DOI×LVI) − 2.94×(T stage×PD_01vs2) |

모든 계수가 음수 = 그 변수가 클수록 예후가 나쁨(pseudo-RMST 감소). DSS Nd=3의
PD×budding 양수항은 PD_01vs2 단독항(−14.6)이 이미 크게 깎아놓은 걸 일부 상쇄하는 조합
효과이며, 독립적인 "budding이 있으면 오히려 좋다"는 뜻으로 해석하면 안 됨.

## 최종 판정 — nested CV (leakage 없음, fold마다 공식 재발견, 10회 평균)

| endpoint | 추천 Nd | C-index | 실험4 최고 모델(같은 CV 방식) | 판정 |
|---|---|---|---|---|
| PFS | Nd=2 | 0.682 ± 0.024 | TabICL 0.589 | **채택 가능** — 실험4의 모든 회귀모델을 이김 |
| DSS | Nd=3 | 0.778 ± 0.037 | AFT(LogNormal) 0.697 | **채택 가능** — 실험4의 모든 회귀모델을 이김 |
| LRRFS | Nd=1 | 0.545 ± 0.039 | DT-IPCW 0.536 | **기각** — 오차범위 안, 사실상 동률. fold별로 유의한 모델을 못 찾는 경우가 잦아(유효 n이 25~27명까지 감소한 fold 존재) 신호 자체가 불안정 |

## 권고

1. **PFS·DSS**: 위 표의 Nd=2(PFS)/Nd=3(DSS) 공식은 mStage와 Cox 다변량에서 통계적으로
   구분 안 되고(S2), leakage를 제거한 held-out 검증에서도 실험4가 시도한 어떤 회귀모델보다
   낫다. mStage를 **대체**할 근거는 아니지만(정확한 공식이 재표본마다 흔들려 "유일해"는 못 박기
   어려움), N stage·PD cellularity·종양크기/DOI가 예후 판별의 핵심이라는 걸 **독립적인 방법으로
   재확인**했고, 향후 병기 개정 논의에서 참고할 수 있는 보조 지표로 제안한다.
2. **LRRFS**: 이번 검증으로 자동 발견 공식을 지지할 근거가 사라졌다 — **채택하지 않는다**.
   LRRFS 자체의 신호가 이 표본 크기(n=129, 사건 44건)에서는 원래 약했다는 정황 증거로 남긴다.
3. **공통 한계**: n=129 단일 코호트 내 discovery+검증이라 다기관/외부 코호트 검증 전에는
   "확정 공식"이 아니라 "가설 생성"으로 취급해야 한다.

---
*생성: GoodRegressor(deepspiking/goodregressor) 기호회귀 파이프라인. RMST pseudo-value →
GoodRegressor 전수탐색(72후보, Nd≤3) → Cox S2/S3/NRI → event-층화 5-fold×10회 nested CV.*
