# 🧬 End-to-End 의료 오믹스 데이터 분석 파이프라인 아키텍처

본 프로젝트는 개별적으로 동작하는 11개의 분석 스크립트들을 유기적으로 연결하여, **원시 데이터(Raw Data) 입력부터 최종 임상 논문용 시각화 출력까지 전 과정을 자동화**한 통합 분석 파이프라인(Master Pipeline)을 구축했습니다.

`src/main_pipeline.py`를 단 한 번 실행하는 것만으로, 아래의 **4단계(Phase) 핵심 연구 파이프라인**이 순차적으로 가동됩니다.

---

## 🌊 파이프라인 워크플로우 (Pipeline Workflow)

### Phase 1: 데이터 탐색 및 기초 통계 (Exploratory Data Analysis & Basic Stats)
임상 정보와 다중 오믹스(Multi-omics) 데이터를 병합한 후, 데이터의 전반적인 분포와 환자 그룹 간의 기초적인 생물학적 차이를 탐색합니다.
* **`09_dimensionality_reduction.py` (PCA/t-SNE)**: 
  * 고차원 유전자 발현 데이터를 2차원으로 축소하여, 고면역군/저면역군 등의 환자 그룹이 생물학적으로 뚜렷하게 군집(Clustering)을 이루는지 거시적으로 파악합니다.
* **`01_comparative_analysis.py` (비교 분석)**: 
  * 환자의 병기나 등급에 따라 특정 바이오마커의 발현량이 통계적으로 유의미한 차이(P-value)를 보이는지 탐색합니다.
* **`02_categorical_correlation.py` (상관 및 범주 분석)**: 
  * 환자의 임상 그룹 간 분포 비율을 확인하고, 여러 바이오마커와 생존 레이블(`OS_event`) 간의 상관관계 네트워크를 히트맵으로 시각화하여 다중공선성을 사전 점검합니다.

### Phase 2: 핵심 바이오마커 발굴 및 기전 분석 (Biomarker Discovery & Mechanism)
수만 개의 노이즈 데이터 속에서 환자의 예후를 결정짓는 핵심 유전자(Signature)를 찾아내고, 그 생물학적 원인을 규명합니다.
* **`08_hierarchical_clustering.py` (계층적 군집화)**: 
  * 발굴된 유전자들의 패턴을 바탕으로 환자들 간의 종양 이질성(Tumor Heterogeneity)을 히트맵으로 시각화합니다.
* **`11_feature_selection_lasso.py` (LASSO 피처 셀렉션)**: 
  * 기계학습 기반 L1 정규화를 통해 생존 예측에 가장 핵심적인 소수 정예의 '다중 유전자 마커 패널(Multi-gene signature)'을 추출합니다.
* **`10_pathway_enrichment.py` (패스웨이 분석)**: 
  * 추출된 핵심 유전자들이 mTOR, MAPK 등 어떤 특정 생물학적 기전(Pathway)의 고장과 연관되어 있는지 통계적으로 입증합니다.

### Phase 3: 임상 예측 모델링 및 생존 분석 (Clinical Prediction & Survival)
발굴된 바이오마커 패널을 바탕으로 실제 환자의 예후(생존/사망)를 예측하는 임상 통계 모델을 구축합니다.
* **`03_survival_analysis.py` (생존 분석)**: 
  * 발굴된 마커를 기준으로 환자군을 나누고, 시간에 따른 생존 확률의 차이(Kaplan-Meier)와 사망 위험비(Cox HR)를 도출합니다.
* **`04_predictive_modeling.py` (진단 모델링)**: 
  * 바이오마커 수치를 결합하여 예측 알고리즘을 만들고, 의사가 직관적으로 진단에 활용할 수 있도록 명시적인 컷오프(Cut-off)를 제공하는 의사결정나무(Decision Tree)를 생성합니다.

### Phase 4: 임상적 유용성 검증 및 고급 시각화 (Clinical Utility & Validation)
구축된 예측 모델이 실제 병원과 진료실에서 얼마나 믿을만하고 이득이 되는지(Clinical Utility)를 다각도로 평가합니다.
* **`06_causal_inference.py` (인과 추론 및 PSM)**: 
  * 후향적 데이터의 한계인 선택 편향을 성향점수매칭(PSM)으로 통제하여, 모델의 결과가 인과적 신뢰성을 갖도록 증명합니다 (Love Plot).
* **`05_concordance_analysis.py` (일치도 평가)**: 
  * 새로운 마커 측정법이 기존 표준 검사법과 구조적인 오차 없이 일치하는지(Bland-Altman) 검증합니다.
* **`07_advanced_visualization.py` (고급 시각화)**: 
  * 최종적으로 의사가 진료실에서 환자의 생존율을 직접 수기 계산할 수 있는 **예측 노모그램(Nomogram)**, 모델의 정확도를 평가하는 **보정 곡선(Calibration)**, 그리고 치료 결정 시 실질적 이득을 보여주는 **임상 결정 곡선(DCA Plot)**을 출력하여 논문의 마침표를 찍습니다.

---

## 🛠 파이프라인 실행 방법

이 모든 4단계 과정은 다음 명령어 단 한 줄로 자동 실행되며, 모든 결과물은 `outputs/` 폴더에 일목요연하게 저장됩니다.

```bash
python src/main_pipeline.py
```
