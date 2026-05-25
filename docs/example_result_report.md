# 📊 단일 스트림 바이오마커 발굴 및 이중 검증 리포트 (Single-Stream Biomarker Pipeline)

이 리포트는 `src/comprehensive_biomarker_pipeline.py`를 구동하여, **데이터 로드부터 이중 피처 셀렉션(LASSO + Decision Tree), 머신러닝 예측(DT & LR), 그리고 생존 및 임상 유용성 검증까지 모든 과정이 하나로 꿰어진(Single Stream)** 최종 분석 결과서입니다. 

---

## 📌 분석 디자인 (Target & Feature Definition)
* **타겟 변수 (y)**: `ESTIMATE_ImmuneScore`의 중앙값을 기준으로 나눈 **고면역군(1) vs 저면역군(0)** 이진 분류.
* **독립 변수 (X)**: 수만 개의 RNA-seq 유전자 발현량 전체.

---

## 1. 유효한 변수 선택 (1차 Feature Selection)
수만 개의 유전자(X) 중에서 타겟(y)을 가장 잘 분류할 수 있는 핵심 유전자를 1차적으로 선별합니다. 단순 LASSO만 돌리지 않고, 타겟(y)과의 상관관계(Pearson Correlation, p<0.01)로 1차 필터링을 거친 후 **LASSO L1-정규화**를 통해 10개의 최정예 후보군을 추려냅니다.

* **LASSO Coefficient Plot**
  수만 개의 노이즈를 뚫고 1차적으로 선별된 10개의 바이오마커 후보군입니다.
  ![LASSO](../outputs/comprehensive_analysis/1_lasso_coefficients.png)

* **선택된 10개 변수 간의 상관관계 히트맵 (Correlation Heatmap)**
  이 10개의 후보 유전자들이 타겟(`Target`)과 어떻게 연관되며, 본인들끼리 어떻게 동반 발현(Co-expression)하는지 다중공선성과 네트워크를 보여줍니다.
  ![Correlation](../outputs/comprehensive_analysis/1_correlation_heatmap.png)

---

## 2. 진단 모델 구축 및 최종 변수 확정 (Diagnostic Modeling)
앞서 1번 단계에서 선별된 10개의 유전자를 입력값으로 받아 두 가지 머신러닝 진단 모델을 구축합니다. 이 과정에서 의사결정나무의 특성을 활용해 최종 바이오마커를 한 번 더 압축(가지치기)합니다.

**1) 의사결정나무 (Decision Tree) 및 최종 바이오마커 확정**
* 10개의 유전자를 모두 넣고 트리를 학습시켰으나, 트리가 불순도(Gini)를 낮추기 위해 **실제 가지치기(Split)에 사용한 핵심 유전자는 단 4개뿐**이었습니다. 
* 따라서 우리는 이 4개의 유전자(`ENSG00000172716.16` 등)만을 **'최종 바이오마커(Final Biomarkers)'**로 확정 짓습니다.
![Decision Tree](../outputs/comprehensive_analysis/2_decision_tree.svg)

**2) 로지스틱 회귀 (Logistic Regression)** 
* 확정된 4개의 최종 유전자만으로 LR 모델을 재학습시켰으며, 각 유전자가 고면역군 진단에 미치는 가중치(Coefficient) 및 오즈비(Odds Ratio)를 표로 산출했습니다.
![LR Coefficients](../outputs/comprehensive_analysis/2_lr_coefficients_table.png)

---

## 3. 최종 바이오마커에 대한 다각도 검증 (Validation)
최종 압축해낸 4개의 유전자와 두 가지 진단 모델(DT, LR)이 실제로 훌륭한지 엄격한 검증을 거칩니다.

**1) 계층적 군집화 (Hierarchical Clustering)**
* 수만 개가 아닌 **오직 최종 4개의 유전자 발현 패턴**만으로 환자들을 바둑판 클러스터링 했을 때, 타겟 그룹(고/저면역군)이 뚜렷하게 나뉘는지 종양 이질성을 시각화합니다.
![Clustering](../outputs/comprehensive_analysis/3b_clustering.png)

**2) 바이오마커만으로의 PCA 및 t-SNE**
* 오직 4개 유전자의 다차원 공간을 2차원으로 축소했을 때, 붉은색(고면역)과 푸른색(저면역) 점들이 섬처럼 확연히 분리되는지 확인합니다.
![PCA t-SNE](../outputs/comprehensive_analysis/3f_pca_tsne.png)

**3) 두 모델(DT vs LR) 예측 점수에 대한 그룹 별 t-test**
* 결정나무(Leaf Probability)와 로지스틱 회귀(Predicted Probability)가 뱉어낸 환자별 위험도 점수가, 실제 그룹 간에 확연한 차이를 보이는지(두 모델의 성능 대결) 증명합니다.
![Model Score t-test](../outputs/comprehensive_analysis/3c_model_scores_ttest.png)

**4) 개별 바이오마커 그룹 별 t-test**
* 최종 선택된 유전자 각각이 저면역군(0)과 고면역군(1)에서 실제로 유의미한 발현 차이를 가지는지 개별 박스플롯(t-test p-value)으로 교차 검증합니다.
![Biomarker t-test](../outputs/comprehensive_analysis/3d_biomarkers_ttest.png)

**5) 바이오마커 기반 그룹의 생존 분석 (Survival Analysis)**
* 두 모델이 예측한 확률로 나눈 위험군이, 단순히 면역도를 넘어 환자의 **실제 생존 기간(OS_days)**마저도 유의미하게 가르는지 2개의 카플란-마이어 곡선을 나란히 비교하여 증명합니다.
![KM Plot](../outputs/comprehensive_analysis/3e_survival_km.png)

**6) 임상 결정 곡선(DCA) 및 보정 곡선 (Calibration)**
* **DCA Plot (우측)**: DT(녹색선)와 LR(파란선) 중 어떤 모델을 믿고 진단했을 때 의사에게 더 높은 순이익(Net Benefit)을 가져다주는지 전격 비교합니다.
* **Calibration (좌측)**: 각 모델의 예측 확률이 뻥튀기되지 않고 45도 점선에 가깝게 실제와 잘 들어맞는지 평가합니다.
![DCA Calibration](../outputs/comprehensive_analysis/3g_dca_calibration.png)

**7) 볼케이노 플롯 (Volcano Plot) 위상 확인**
* 전체 수만 개 유전자의 화산 플롯 상에서, 우리가 뽑은 4개의 최종 바이오마커(빨간 점, 라벨링)가 실제로 가장 양쪽 상단 끝(최우수 유의성)에 위치하고 있음을 시각적 쐐기로 박습니다.
![Volcano Plot](../outputs/comprehensive_analysis/3h_volcano_plot.png)
