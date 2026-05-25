# 📊 단일 스트림 바이오마커 발굴 및 이중 검증 리포트 (Single-Stream Biomarker Pipeline)

이 리포트는 `src/comprehensive_biomarker_pipeline.py`를 구동하여, **데이터 로드부터 이중 피처 셀렉션(LASSO + Decision Tree), 머신러닝 예측(DT & LR), 그리고 생존 및 임상 유용성 검증까지 모든 과정이 하나로 꿰어진(Single Stream)** 최종 분석 결과서입니다. 

사용자님의 강력한 피드백을 수용하여, **"단순히 묶어둔 파이프라인"을 넘어 "의사결정나무가 리프(Leaf) 노드 분기에 실제 사용한 최정예 변수들만 역추적하여 최종 바이오마커로 확정하고, 이를 기반으로 로지스틱 회귀 계수 테이블 및 두 모델(DT, LR)의 스코어를 동시 비교 검증하는 구조"**로 전면 업그레이드 되었습니다.

---

## 📌 분석 디자인 (Target & Feature Definition)
* **타겟 변수 (y)**: `ESTIMATE_ImmuneScore`의 중앙값을 기준으로 나눈 **고면역군(1) vs 저면역군(0)** 이진 분류.
* **독립 변수 (X)**: 수만 개의 RNA-seq 유전자 발현량 전체.

---

## 1. 2단계 유효 변수 선택 (Dual Feature Selection)
수만 개의 유전자(X) 중에서 타겟(y)을 가장 잘 분류할 수 있는 핵심 유전자를 선별합니다. 단순 LASSO만 돌리지 않고, 모델의 작동 원리를 활용하여 마커를 압축합니다.

**1차 선택 (LASSO 기반)**
상관관계(p<0.01)로 1차 필터링된 유전자들을 L1-정규화(LASSO)에 통과시켜 후보군을 추려냅니다.
![LASSO](../outputs/comprehensive_analysis/1_lasso_coefficients.png)

**2차 선택 (Decision Tree 분기 역추적)**
LASSO 후보군을 다시 의사결정나무(Decision Tree)에 학습시킵니다. 트리가 불순도(Gini)를 낮추기 위해 **실제 가지치기(Split)에 사용한 핵심 유전자들(최종 4개)**만을 완벽한 '최종 바이오마커(Final Biomarkers)'로 선포합니다. 이들끼리의 생물학적 동반 발현(Correlation) 네트워크는 아래와 같습니다.
![Correlation](../outputs/comprehensive_analysis/1_correlation_heatmap.png)

---

## 2. 진단 모델 구축 (Diagnostic Modeling)
최종 선별된 4개의 바이오마커만을 입력값으로 받아 두 가지 머신러닝 모델을 구축합니다.

**1) 의사결정나무 (Decision Tree)**: 
의사가 직관적으로 "이 유전자의 발현량이 얼마 이상이면 고면역군이다"라는 명시적인 진단 기준(Cut-off)을 볼 수 있게 해줍니다.
![Decision Tree](../outputs/comprehensive_analysis/2_decision_tree.svg)

**2) 로지스틱 회귀 (Logistic Regression)**: 
이 4개의 최종 유전자만으로 LR 모델을 재학습시켰으며, 각 유전자가 고면역군 진단에 미치는 가중치(Coefficient) 및 오즈비(Odds Ratio)를 표로 산출했습니다.
![LR Coefficients](../outputs/comprehensive_analysis/2_lr_coefficients_table.png)

---

## 3. 해당 바이오마커에 대한 다각도 검증 (Validation)
우리가 압축해낸 이 4개의 유전자와 두 가지 진단 모델(DT, LR)이 실제로 훌륭한지 8가지 엄격한 검증을 거칩니다.

**1) 개별 바이오마커 그룹 별 t-test**
* 최종 선택된 유전자 각각이 저면역군(0)과 고면역군(1)에서 실제로 유의미한 발현 차이를 가지는지 개별 박스플롯(t-test p-value)으로 1차 검증합니다.
![Biomarker t-test](../outputs/comprehensive_analysis/3d_biomarkers_ttest.png)

**2) 계층적 군집화 (Hierarchical Clustering)**
* 수만 개가 아닌 **오직 이 4개의 유전자 발현 패턴**만으로 환자들을 바둑판 클러스터링 했을 때, 타겟 그룹(고/저면역군)이 뚜렷하게 나뉘는지 종양 이질성을 시각화합니다.
![Clustering](../outputs/comprehensive_analysis/3b_clustering.png)

**3) 바이오마커만으로의 PCA 및 t-SNE**
* 오직 4개 유전자의 다차원 공간을 2차원으로 축소했을 때, 붉은색(고면역)과 푸른색(저면역) 점들이 섬처럼 확연히 분리되는지 확인합니다.
![PCA t-SNE](../outputs/comprehensive_analysis/3f_pca_tsne.png)

**4) 두 모델(DT vs LR) 예측 점수에 대한 그룹 별 t-test**
* 결정나무(Leaf Probability)와 로지스틱 회귀(Predicted Probability)가 뱉어낸 환자별 위험도 점수가, 실제 그룹 간에 확연한 차이를 보이는지(두 모델의 성능 대결) 증명합니다.
![Model Score t-test](../outputs/comprehensive_analysis/3c_model_scores_ttest.png)

**5) 볼케이노 플롯 (Volcano Plot) 위상 확인**
* 전체 수만 개 유전자의 화산 플롯 상에서, 우리가 방금 뽑은 4개의 최종 바이오마커(빨간 점, 라벨링)가 실제로 가장 양쪽 상단 끝(최우수 유의성)에 위치하고 있음을 시각적 쐐기로 박습니다.
![Volcano Plot](../outputs/comprehensive_analysis/3h_volcano_plot.png)

**6) 바이오마커 기반 그룹의 생존 분석 (Survival Analysis)**
* 로지스틱 회귀(LR)와 의사결정나무(DT)가 예측한 위험군 라벨이 실제 환자의 **생존 기간(OS_days)**마저도 유의미하게 가르는지 각각 2개의 카플란-마이어(Kaplan-Meier) 곡선을 나란히 배치하여 모델 간의 생존 예측 능력을 전격 비교합니다.
![KM Plot](../outputs/comprehensive_analysis/3e_survival_km.png)

**7) 임상 결정 곡선(DCA) 및 보정 곡선 (Calibration)**
* **DCA Plot (우측)**: DT(녹색선)와 LR(파란선) 중 어떤 모델을 믿고 진단했을 때 의사에게 더 높은 순이익(Net Benefit)을 가져다주는지 두 모델을 전격 비교합니다.
* **Calibration (좌측)**: 각 모델의 예측 확률이 뻥튀기되지 않고 45도 점선에 가깝게 실제와 잘 들어맞는지 평가합니다.
![DCA Calibration](../outputs/comprehensive_analysis/3g_dca_calibration.png)
