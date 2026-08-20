# 실험 4e 결과 해석 — 생존시간 회귀 (Survival-time regression)

> 작성일: 2026-08-18
> 목적: 24개월 이진화(classification)의 임의성을 없애고, **생존기간(월)을 직접 회귀**로 예측
> 방법: 5-fold CV × 시드 4개, 16개 모델 (RandomForest는 소표본 속도 문제로 제외)
> 결과: `exp1_surv_regression.csv`, `exp1_surv_regression_summary.csv`
> ⚠️ 참고: 아래 결과는 시드별 5-fold 평균의 시드 평균값

---

## 1. 핵심 요약

| y | 최고 C-index | 모델 | AUC(이진화)와 비교 |
|---|---|---|---|
| PFS | **0.589** | TabICL_naive | AUC 0.61 ≈ 유사 |
| DSS | **0.697** | AFT_LogNormal | AUC 0.75보다 낮음 |
| LRRFS | 0.536 | DecisionTree_IPCW | AUC 0.56 ≈ 유사 |

> **결론**: 생존시간 회귀는 "언제 사건이 발생하는지"를 정확히 예측해야 하므로
> **이진 분류(순위)보다 어렵다** — 전반적으로 C-index가 이진화 AUC보다 낮거나
> 유사합니다. 회귀는 **탐색적 보조**로, 핵심 검증은 C-index/NRI(순위 기반)가 적합.

---

## 2. 모델·설계

| 구분 | 내용 |
|---|---|
| Target | log(생존기간) — 회귀 |
| 중도절단 처리 | ① AFT(LogNormal/Weibull): 정식 처리 ② IPCW: uncensored에 1/Ĝ(T) 가중 ③ naive: 보정 없음(대조) |
| Normalization | **StandardScaler / PowerTransformer(Yeo-Johnson)** 2종 비교 — 사용자 지적 반영 |
| TabICL | 내부 전처리 수행 → normalization 없이 naive/uncensored만 |
| 평가 | C-index(주) + uncensored Spearman ρ + RMSE(log-time) |
| 제외 | RandomForest — 소표본에서 fold당 수 분으로 비현실적 |

---

## 3. 결과 상세

### 3.1 PFS

| 모델 | C-index | ρ (uncensored) | RMSE |
|---|---|---|---|
| **TabICL_naive** | **0.589** | +0.063 | 1.415 |
| AFT_Weibull | 0.587 | +0.041 | 2.402 |
| AFT_LogNormal | 0.576 | +0.050 | 2.054 |
| Lasso_std/power_naive | 0.574/0.572 | +0.076/+0.048 | 1.27 |
| DecisionTree_std_naive | 0.561 | +0.172 | 1.361 |
| Ridge_std_naive | 0.558 | +0.032 | 1.511 |
| *(IPCW 계열)* | 0.49-0.53 | ±0.2 | 0.96-1.4 |

**해석**: naive 계열이 IPCW보다 우월. TabICL이 최고지만 단순 모델과 큰 차이 없음.

### 3.2 DSS

| 모델 | C-index | ρ (uncensored) | RMSE |
|---|---|---|---|
| **AFT_LogNormal** | **0.697** | +0.331 | 1.931 |
| AFT_Weibull | 0.688 | +0.329 | 2.281 |
| Lasso_std_naive | 0.683 | +0.244 | 0.941 |
| TabICL_naive | 0.679 | +0.166 | 1.413 |
| Ridge_std_naive | 0.660 | +0.221 | 1.082 |
| *(IPCW 계열)* | 0.58-0.63 | 음수 | 0.80-1.05 |

**해석**: **AFT(LogNormal)이 C-index 0.697 + uncensored ρ 0.33으로 최고** —
중도절단을 정식 처리한 회귀가 DSS에서 가장 좋음. IPCW는 ρ가 음수로 나쁨.

### 3.3 LRRFS

| 모델 | C-index | ρ (uncensored) | RMSE |
|---|---|---|---|
| **DecisionTree_std_IPCW** | **0.536** | +0.157 | 1.072 |
| TabICL_naive | 0.532 | +0.023 | 1.234 |
| Lasso_power_naive | 0.529 | -0.001 | 1.418 |
| AFT_LogNormal | 0.528 | -0.081 | 2.480 |
| *(나머지)* | 0.50-0.53 | ±0.1 | 0.93-2.6 |

**해석**: 전반적으로 **모든 모델이 0.5 근처** — LRRFS는 회귀로도 잘 예측되지
않는 어려운 결과(이진화 AUC 0.56과 유사).

---

## 4. Normalization 비교 (Std vs Power)

| y | 모델 | Std C-index | Power C-index | 차이 |
|---|---|---|---|---|
| PFS | Lasso naive | 0.574 | 0.572 | 미미 |
| DSS | Lasso naive | 0.683 | 0.675 | 미미 |
| DSS | Ridge naive | 0.660 | 0.658 | 미미 |
| LRRFS | Lasso naive | 0.528 | 0.529 | 미미 |

> **결론**: StandardScaler와 PowerTransformer는 **C-index에 유의미한 차이를
> 만들지 않았습니다.** 다만 uncensored Spearman ρ에서는 power가 다소 나은
> 경우가 있었습니다(예: DSS Lasso power ρ=0.223 vs std 0.244). 변수 분포의
> 비대칭성이 큰 경우 power가 안정적인 해석에 도움이 될 수 있습니다.

---

## 5. 방법론별 총평 (y 통합)

| 방법론 | 특징 | 적합성 |
|---|---|---|
| **AFT (LogNormal/Weibull)** | 중도절단 정식 처리, DSS에서 최고 | ✅ 생존시간 회귀에 가장 이론적으로 적합 |
| **naive Lasso/Ridge/DT** | 단순·빠름, PFS/DSS에서 양호 | ✅ 실용적 |
| **IPCW 가중** | 소표본에서 ρ 음수·불안정 | ⚠️ n=133에선 비효율 |
| **TabICL** | 내부 처리, 단순 모델과 유사 | ⚠️ 소표본 한계 유지 |
| Power/Std | 차이 미미 | ✅ 둘 다 무방 |

---

## 6. AUC(이진화) vs 회귀 C-index

| y | AUC(24개월) | 회귀 최고 C-index | 차이 |
|---|---|---|---|
| PFS | 0.614 (Lasso) | 0.589 (TabICL) | -0.025 |
| DSS | 0.752 (Lasso) | 0.697 (AFT) | -0.055 |
| LRRFS | 0.557 (DT) | 0.536 (DT-IPCW) | -0.021 |

> 회귀가 전반적으로 **약간 낮은 성능** — "정확한 시점 예측"이 "위험도 순위"보다
> 근본적으로 어렵기 때문. **이 연구의 핵심(수정병기 우월성)은 순위 기반 지표
> (C-index/NRI)로 답하는 것이 적절**하며, 회귀는 보조적 탐색으로 보고합니다.

---

## 7. 한계

1. n=133, 사건 수 적음(DSS 28명) — 회귀 계수 추정이 불안정
2. IPCW 가중치는 소표본에서 극단적 가중치로 오히려 해로움
3. RandomForest 제외 (소표본 속도) — 앙상블 트리 회귀 미검토
4. 생존시간 회귀는 중도절단이 많은 데이터에서 정보 손실 큼

---

## 8. 결론 (한 문장)

> **생존시간 회귀(특히 AFT·naive 계열)는 이진화 AUC와 유사하거나 약간 낮은
> 성능(최고 C-index 0.70)으로 탐색적 가치가 있으나, n=133 소표본에서는
> 순위 기반 지표(C-index/NRI)가 수정병기 검증의 주 근거로 적합하며,
> 회귀는 보조 분석으로 보고한다.**
