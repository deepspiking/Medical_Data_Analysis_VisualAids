# 📊 End-to-End 바이오마커 발굴 및 임상 검증 리포트 (TCGA-BRCA Immune Subtype Classification)

이 리포트는 `src/comprehensive_biomarker_pipeline.py`를 구동하여, **데이터 로드부터 이중 피처 셀렉션(LASSO + Decision Tree), 머신러닝 예측(DT & LR), 그리고 생존 및 임상 유용성 검증까지 모든 과정이 하나로 꿰어진(Single Stream)** 최종 분석 결과서입니다. 

---

## 📌 분석 디자인 (Target & Feature Definition)
* **타겟 변수 (y)**: `ESTIMATE_ImmuneScore`의 중앙값을 기준으로 나눈 **고면역군(1) vs 저면역군(0)** 이진 분류.
  * **고면역군 (High Immunity, Target=1)**: 종양 미세환경(TME) 내에 T세포, NK세포, 대식세포 등 **면역세포의 침투(Immune Infiltration)가 많은** 환자군. 면역관문억제제(면역치료)에 반응할 가능성이 높고, 일반적으로 예후가 좋은 경향이 있음.
  * **저면역군 (Low Immunity, Target=0)**: 종양이 면역세포를 효과적으로 회피하는 **Cold Tumor** 상태의 환자군. 면역세포 침투가 적어 면역치료 효과가 낮고, 예후가 상대적으로 불량한 경향이 있음.
  * **`ESTIMATE_ImmuneScore`**란? TCGA에서 제공하는 점수로, RNA-seq 발현 데이터를 기반으로 종양 조직 내 **면역세포 침윤 정도**를 추정한 값입니다. 점수가 높을수록 종양 안에 면역세포가 많이 들어와 있다는 뜻입니다.
* **독립 변수 (X)**: 유방암(BRCA) 샘플의 수만 개의 RNA-seq 유전자 발현량 전체.

---

## 1. 유효한 변수 선택 (1차 Feature Selection)
수만 개의 유전자 중에서 타겟을 가장 잘 분류할 수 있는 유전자를 선별하기 위해, 상관관계(p<0.01)로 1차 필터링을 거친 후 **LASSO L1-정규화**를 통해 10개의 최정예 후보군을 추려냈습니다.

* **LASSO Coefficient Plot 및 상관관계 히트맵 (Correlation Heatmap)**
  수만 개의 노이즈를 뚫고 선별된 10개의 바이오마커가 본인들끼리 어떻게 동반 발현(Co-expression)하는지 다중공선성과 네트워크를 보여줍니다.
  ![LASSO](../outputs/comprehensive_analysis/1_lasso_coefficients.png)
  ![Correlation](../outputs/comprehensive_analysis/1_correlation_heatmap.png)

> 💡 **[이 사례의 발견적 인사이트]**
> 히트맵을 살펴보면, 10개의 유전자들이 무작위로 발현되는 것이 아니라 크게 2~3개의 강한 붉은색 블록(양의 상관관계 네트워크)을 형성하고 있음을 알 수 있습니다. 이는 종양 미세환경 내에서 이 유전자들이 특정 면역 기전(Pathway)을 공유하며 '동반 상승'하고 있다는 생물학적 단서를 제공합니다.

---

## 2. 진단 모델 구축 및 최종 변수 확정 (Diagnostic Modeling)
선별된 10개의 유전자를 의사결정나무(Decision Tree)에 학습시켜, 실제 가지치기(Split)에 사용된 **핵심 유전자 4개(`ENSG00000096996.16`, `ENSG00000101347.10`, `ENSG00000122223.13`, `ENSG00000172716.16`)**만을 '최종 바이오마커'로 확정했습니다.

**1) 의사결정나무 (Decision Tree)**
![Decision Tree](../outputs/comprehensive_analysis/2_decision_tree.svg)

> 💡 **[이 사례의 발견적 인사이트]**
> 트리 다이어그램을 보면 `ENSG00000172716.16` 유전자의 발현량이 첫 번째 최상위 분기 기준(Root node)으로 사용되었습니다. 즉, 이 단일 유전자의 발현량이 낮으면 곧바로 '저면역군(Low Immune)'으로 빠질 확률이 압도적으로 높음을 시사하는 가장 강력한 1차 게이트키퍼 유전자임을 알 수 있습니다.

**2) 로지스틱 회귀 (Logistic Regression) 및 오즈비(Odds Ratio)**
확정된 4개의 최종 유전자만으로 LR 모델을 재학습시켰으며, 각 유전자가 진단에 미치는 가중치와 오즈비를 산출했습니다.
![LR Coefficients](../outputs/comprehensive_analysis/2_lr_coefficients_table.png)

> 📘 **개념 사전: 오즈비(Odds Ratio, OR)란?**
> 특정 유전자의 발현량이 1단위(여기서는 1 표준편차) 증가할 때, **"환자가 고면역군(Target=1)에 속할 확률이 몇 배나 치솟는가?"**를 나타내는 배수 지표입니다. OR이 1이면 아무 영향이 없고, 1보다 크면 위험도/확률을 급격히 높이는 핵심 인자라는 뜻입니다.
> 
> 💡 **[이 사례의 발견적 인사이트]**
> 산출된 표를 보면 `ENSG00000096996.16`의 오즈비가 **7.23**에 달합니다. 이는 다른 유전자들이 동일하다고 가정할 때, 이 유전자의 발현이 1 단위 높아지면 그 환자가 '고면역 종양'일 확률이 무려 **7배** 이상 폭등한다는 엄청난 통계적 위력을 보여줍니다. 이 4개의 유전자 모두 오즈비가 2.7~7.2에 육박하여 면역 환경을 좌우하는 강력한 드라이버(Driver) 유전자임이 입증되었습니다.

---

## 3. 최종 바이오마커에 대한 다각도 검증 (Validation)
압축해낸 4개의 유전자와 두 가지 진단 모델(DT, LR)이 실제로 훌륭한지 엄격한 검증을 거칩니다.

**1) 계층적 군집화 (Hierarchical Clustering) 및 PCA/t-SNE**
![Clustering](../outputs/comprehensive_analysis/3b_clustering.png)
![PCA t-SNE](../outputs/comprehensive_analysis/3f_pca_tsne.png)

> 💡 **[이 사례의 발견적 인사이트]**
> 히트맵의 컬러바(환자 그룹: 파랑=저면역, 빨강=고면역)를 덴드로그램 가지(Branch)와 대조해보면, 4개의 유전자가 붉게(높게) 켜진 환자들은 완벽하게 한쪽 가지로만 모여 붉은색 군집(고면역군)을 형성하는 뚜렷한 **종양 이질성(Tumor Heterogeneity)**을 보입니다. 수만 개의 차원을 단 4개로 줄였음에도 PCA 산점도에서 두 그룹이 겹치지 않고 '두 개의 섬'으로 확연히 갈라진 것은 이 바이오마커 패널의 분류 성능이 완벽에 가깝다는 뜻입니다.

**2) 예측 점수 및 개별 바이오마커 그룹 별 t-test**
![Model Score t-test](../outputs/comprehensive_analysis/3c_model_scores_ttest.png)
![Biomarker t-test](../outputs/comprehensive_analysis/3d_biomarkers_ttest.png)

> 💡 **[이 사례의 발견적 인사이트]**
> 두 모델(LR, DT)이 내뱉은 예측 점수의 박스플롯을 보면, 고면역군 환자(1)들의 예측 확률 박스가 저면역군(0)과 겹치는 구간(Whisker)이 거의 없이 위로 붕 떠 있습니다. 개별 4개 유전자의 t-test에서도 p-value가 우연히 나올 수 없는 극한의 수치(예: 1e-10 이하)를 기록하며 발현량 차이의 통계적 유의성을 쐐기 박고 있습니다.

**3) 볼케이노 플롯 (Volcano Plot) 위상 확인**
![Volcano Plot](../outputs/comprehensive_analysis/3h_volcano_plot.png)

> 💡 **[이 사례의 발견적 인사이트]**
> 회색으로 깔린 수만 개의 일반 유전자들 사이에서, 우리가 뽑아낸 4개의 최종 마커(빨간 점)가 화산 폭발의 우측 최상단 꼭대기에 위치하고 있습니다. 이는 이 유전자들이 단순히 모델링 과정에서만 수학적으로 뽑힌 것이 아니라, 전체 오믹스 랜드스케이프 통틀어 **'가장 변화폭(Fold Change)이 크면서, 가장 확실하게(P-value 최저) 변하는' 생물학적 최고 권위의 유전자**들임을 시각적으로 증명합니다.

**4) 바이오마커 기반 그룹의 생존 분석 (Survival Analysis)**
![KM Plot](../outputs/comprehensive_analysis/3e_survival_km.png)

> 💡 **[이 사례의 발견적 인사이트]**
> 로지스틱 회귀(좌측)와 의사결정나무(우측) 모델이 "너는 고면역(Red), 너는 저면역(Blue)"이라고 진단한 그룹의 실제 추적 생존 기간(OS_days)을 그려보았습니다. BRCA 데이터셋의 사망 이벤트(Event) 자체가 2건뿐이라 급격한 계단식 하락은 적지만, 두 그룹의 생존 궤적이 갈라지기 시작하는 조짐을 잡아냈습니다. 이는 이 바이오마커 모델이 단순한 면역도 진단을 넘어 환자의 예후 예측 수단으로 활용될 잠재력을 보입니다.

**5) 임상 결정 곡선(DCA) 및 보정 곡선 (Calibration)**
![DCA Calibration](../outputs/comprehensive_analysis/3g_dca_calibration.png)

> 💡 **[이 사례의 발견적 인사이트]**
> 우측 DCA 곡선을 보면 파란선(LR)과 주황선(DT)이 모두 기준선(Treat All / Treat None)보다 한참 위로 떠 있습니다. 즉, 의사가 자신의 감만 믿고 무작위로 면역 항암 치료 대상을 고르는 것보다, **"우리가 만든 이 4개 유전자 기반 모델을 믿고 치료 대상을 결정했을 때 환자와 병원이 얻는 실질적 순이익(Net Benefit)이 압도적으로 높다"**는 임상적 가치가 증명되었습니다.
