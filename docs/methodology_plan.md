# 의료 데이터 분석 및 시각화 구현 계획서

## 개요
이 문서는 제시된 임상 및 중개연구 데이터 분석의 필수 통계 방법론과 핵심 시각화 도구를 Python으로 구현하기 위한 계획서입니다. 분석은 프로젝트 내의 `BRCA` (Breast Invasive Carcinoma) 폴더에 있는 샘플 데이터를 활용합니다.

## 구현 파일 및 포함 내용 (`src/` 폴더 내 작성)

### 1. `01_comparative_analysis.py` (그룹 간 비교 분석)
* **목표**: 정규/비정규 분포를 따르는 두 개 또는 세 개 이상의 집단 간 차이 비교 및 유의성 시각화
* **구현 방법론**:
  * Student's t-test, ANOVA, Mann-Whitney U test, Kruskal-Wallis test
  * **Box Plot (P-value annotation 포함)**: 집단 간 발현량 분포와 통계적 유의성을 한눈에 보여주는 시각화
* **사용 데이터**: `BRCA_phenotype.txt` (임상 그룹 변수), `BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` (유전자 발현량 등 연속형 변수)

### 2. `02_categorical_correlation.py` (범주형 데이터 및 상관관계 분석)
* **목표**: 임상 인자 간의 연관성 및 연속형 변수 간의 상관관계 검증 및 시각화
* **구현 방법론**:
  * Chi-Square test, Fisher's Exact test
  * Pearson 상관분석, Spearman 상관분석
  * **Scatter Plot (Regression line & p-value 포함)**: 두 변수 간의 선형 회귀 추세 및 상관계수 표시
* **사용 데이터**: `BRCA_phenotype.txt` (병기, 성별 등 범주형 변수), `BRCA_proteomics...txt` (단백질 발현 등 연속형 변수)

### 3. `03_survival_analysis.py` (생존 및 예후 분석 + Kaplan-Meier Plot)
* **목표**: 환자의 누적 생존율 추적 및 예후 인자 식별
* **구현 방법론**:
  * Kaplan-Meier 생존 곡선 (Kaplan-Meier Plot)
  * Log-rank test
  * Cox Proportional Hazards Model (Cox 회귀분석)
* **사용 데이터**: `BRCA_survival.txt` (생존 시간 및 상태), `BRCA_phenotype.txt` (비교를 위한 환자 변수)

### 4. `04_predictive_modeling.py` (진단 및 예측 모델링 + Decision Tree Diagram)
* **목표**: 특정 임상 결과 발생 확률 예측 및 분류
* **구현 방법론**:
  * Logistic Regression (로지스틱 회귀)
  * Decision Tree (의사결정나무 모델 및 시각화 Diagram)
  * 다중 모달 융합 딥러닝 템플릿 코드 (구조적 예시 제공)
* **사용 데이터**: `BRCA_phenotype.txt`, 유전자/단백질 발현 데이터 융합

### 5. `05_concordance_analysis.py` (검사법 간 일치도 + Bland-Altman Plot)
* **목표**: 두 진단/판독 간의 일치도 측정 및 바이어스 확인
* **구현 방법론**:
  * Cohen's Kappa 계수
  * Bland-Altman Plot
* **사용 데이터**: (유사한 두 가지 연속 측정치 또는 범주형 측정치 생성하여 비교)

### 6. `06_causal_inference.py` (혼란 인자 통제 + Love Plot)
* **목표**: 관찰 연구의 선택 편향 보정
* **구현 방법론**:
  * Propensity Score Matching (PSM) 및 Love Plot 시각화
  * Inverse Probability of Treatment Weighting (IPTW)
* **사용 데이터**: `BRCA_phenotype.txt` (치료군 vs 대조군 모사)

### 7. `07_advanced_visualization.py` (고급 임상 시각화 도구)
* **목표**: 논문 품질의 필수 시각화 도구 구현
* **구현 방법론**:
  * 예측 노모그램 (Nomogram 템플릿)
  * 보정 곡선 (Calibration Curve)
  * 임상 결정 곡선 (Decision Curve Analysis, DCA Plot)
  * 볼케이노 플롯 (Volcano Plot - DEG 분석 시각화)
  * 구획별 버블 플롯 (Bubble Plot)
  * 설명 가능한 AI (Grad-CAM 템플릿 - 가상의 이미지/텐서 활용)
* **사용 데이터**: `BRCA_RNAseq...`, `BRCA_survival.txt`, 모형 예측 결과 등

### 8. `08_hierarchical_clustering.py` (계층적 군집화)
* **목표**: 채택된 다수의 유전자/바이오마커에 대해 집단(Case) 간 이질성을 파악하는 히트맵 시각화
* **구현 방법론**:
  * Hierarchical Clustering (Clustermap)
  * 환자 그룹별(예: 고병기 vs 저병기) 라벨링을 통한 표현형-유전자 발현 연관성 분석
* **사용 데이터**: `BRCA_phenotype.txt`, `BRCA_RNAseq_gene...` 또는 Proteomics 데이터

## 개발 및 실행 환경
* Python 3.8+
* 필요 라이브러리: `pandas`, `numpy`, `scipy`, `statsmodels`, `scikit-learn`, `lifelines` (생존 분석 용), `matplotlib`, `seaborn`, `psmatch` 등.
* 데이터 로드 시 각 데이터셋의 공통 환자 ID를 기준으로 병합(Merge)하여 사용합니다.
* **결과물 저장**: 모든 스크립트 실행을 통해 생성된 시각화 결과물(Plot 등)은 `outputs/` 폴더 내에 각 분석 카테고리별 하위 폴더(예: `survival_plots/`, `concordance_plots/`, `predictive_modeling/`, `causal_inference/`, `advanced_plots/`)로 구분되어 저장됩니다.