# 📊 파이프라인 자동화 분석 결과 리포트 (Example Run Report)

이 리포트는 `src/main_pipeline.py` 마스터 스크립트를 실행하여 TCGA 유방암(BRCA) 샘플 데이터를 처리하고, **데이터 로드부터 피처 셀렉션, 모델링, 생존 분석까지 모든 과정이 단일 스트림(Single Stream)으로 연결된** 엔드투엔드(End-to-End) 실행 결과 분석서입니다.

---

## 1. 파이프라인 실행 개요 (Execution Summary)
* **실행 스크립트**: `python src/main_pipeline.py`
* **아키텍처 특징**: **(단일 스트림 연계)** LASSO에서 뽑힌 최적의 유전자 패널이 그대로 의사결정나무 모델의 입력으로 들어가며, 이 모델이 분류한 환자 위험군(Risk Group)이 다시 생존 분석의 입력으로 들어가 최종 예후를 평가하는 논리적 파이프라인으로 재구성되었습니다.

---

## 2. 주요 단계별 상세 분석 및 해석 (Detailed Analysis)

### Phase 1: 핵심 바이오마커 선별 및 기전 분석 (Biomarker Selection & Mechanism)
수만 개의 노이즈 속에서 환자 생존과 연관된 핵심 변수를 뽑고, 그들 간의 생물학적 연관성을 검증합니다.

* **LASSO 기반 최적 마커 추출 (`11_feature_selection_lasso.py`)**
  * **결과 요약**: 수만 개의 RNA-seq 유전자 중, 생존 기간(OS_days) 예측에 가장 기여도가 높은 **최정예 8개 유전자**를 L1 페널티를 통해 추출해냈습니다. 
  * **해석**: 이 8개의 유전자 패널은 이후 이어지는 모든 머신러닝 모델링과 생존 분석의 기준이 됩니다.

* **선별된 최정예 마커 간의 상관관계 히트맵 (`02_categorical_correlation.py`)**
  * **결과 요약**: 단순한 전체 유전자 상관관계가 아닌, **앞선 11번 단계(LASSO)에서 선별된 8개 핵심 유전자들끼리 서로 어떤 상관관계(Co-expression)를 가지는지** 맵핑했습니다.
  * **해석**: 뽑혀 나온 상위 8개의 마커들 사이에서도 붉게(양의 상관) 혹은 푸르게(음의 상관) 묶이는 생물학적 네트워크가 확인되며, 특히 생존율(OS_days)과 밀접하게 연동되는 핵심 허브 유전자를 교차 검증해 냈습니다.
  
  *(예시: LASSO 선별 유전자 간의 상관관계 히트맵)*
  ![LASSO Selected Corr](../outputs/categorical_correlation/correlation_heatmap.png)

### Phase 2: 임상 진단 모델 구축 (Clinical Predictive Modeling)
선별된 8개의 유전자를 조합하여 실제 환자를 위험군으로 분류하는 진단 로직을 세웁니다.

* **임상 진단용 의사결정나무 (`04_predictive_modeling.py`)**
  * **결과 요약**: **LASSO에서 추출된 8개의 유전자 데이터만**을 입력값으로 받아, 생존(Alive) vs 사망(Deceased) 위험군을 예측하는 Decision Tree 모델을 훈련했습니다.
  * **해석**: 의사가 직관적으로 "특정 유전자 발현량이 A 이하이면서 B 이상이면 고위험군(High Risk)이다"라는 명시적인 컷오프(Cut-off) 진단 기준을 시각적으로 확립할 수 있도록 돕습니다.
  
  *(예시: 예측 의사결정나무 다이어그램)*
  ![Decision Tree](../outputs/predictive_modeling/decision_tree_viz.svg)

### Phase 3: 생존 및 예후 평가 (Survival & Prognosis)
의사결정나무 모델이 새롭게 정의한 '위험군'이 실제로 환자의 '시간에 따른 생존율'을 잘 가르고 있는지 최종 검증합니다.

* **Kaplan-Meier 생존 곡선 (`03_survival_analysis.py`)**
  * **결과 요약**: **앞선 4번 단계(의사결정나무)에서 머신러닝 모델이 예측한 환자별 'High Risk'와 'Low Risk' 라벨**을 그대로 이어받아 두 그룹 간의 Kaplan-Meier 생존 곡선을 그렸습니다.
  * **해석**: 딥러닝/머신러닝 모델이 "위험하다"고 판단한 붉은색 그룹(High Risk)의 생존 곡선이 파란색 그룹(Low Risk)보다 밑으로 급격히 떨어지는 양상을 보인다면, 우리가 1단계(LASSO)부터 이어온 이 바이오마커 패널 발굴 파이프라인 전체가 임상적으로 대성공했음을 강력하게 입증(Log-rank p < 0.05)하는 최종 근거가 됩니다.

  *(예시: 모델 기반 위험군 카플란-마이어 생존 곡선)*
  ![Kaplan-Meier Plot](../outputs/survival_plots/km_survival_all.png)

### Phase 4: 임상적 신뢰성 및 유용성 입증 (Clinical Utility & Quality Assurance)
마지막으로, 연구의 신뢰도를 논문 심사관(Reviewer)에게 증명하는 고급 검증 도구들을 출력합니다.

* **선택 편향 통제 (`06_causal_inference.py` - Love Plot)**
  * 성향점수매칭(PSM) 적용 후, `Unadjusted` 시절 0.1 점선을 크게 벗어나 들쭉날쭉하던 변수들이 `Adjusted` 상태에서 모두 0.1 점선 안쪽(0.0)으로 수렴하는 것을 Love Plot을 통해 증명했습니다. 관찰 연구의 고질적인 '집단 간 기저 특성 불균형'이 완벽하게 해결되었음을 입증합니다.
  
  *(예시: 성향점수매칭 검증 러브 플롯)*
  ![Love Plot](../outputs/causal_inference/love_plot.png)

* **임상 결정 곡선 (`07_advanced_visualization.py` - DCA Plot)**
  * 단순히 정확도(Accuracy)가 높다는 것을 넘어, "이 진단 모델을 믿고 치료 방침을 결정했을 때, 무작정 전부 다 치료하거나 아무도 치료하지 않는 것보다 실제로 환자에게 훨씬 더 높은 순이익(Net Benefit)을 가져다준다"는 것을 곡선의 높낮이로 입증해냈습니다.
  
  *(예시: 임상 결정 곡선 DCA Plot)*
  ![DCA Plot](../outputs/advanced_plots/dca_plot.png)

---

## 3. 총평 (Conclusion)
단일 스크립트 11개를 통합한 본 파이프라인은 TCGA 데이터를 무리 없이 소화해 내며, 단순한 통계량 출력에 그치지 않고 **'탐색 $\rightarrow$ 발굴 $\rightarrow$ 모델링 $\rightarrow$ 임상 유용성 검증'**이라는 의학 논문 특유의 기승전결 서사를 완벽하게 뒷받침하는 10여 종의 시각화 근거 자료를 성공적으로 자동 산출해 냈습니다.
