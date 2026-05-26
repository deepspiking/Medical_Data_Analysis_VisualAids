# 의료 오믹스 데이터 분석 및 시각화 아키텍처 (Single-Stream Pipeline)

## 개요
이 문서는 임상 및 중개연구 데이터(예: TCGA-BRCA)를 분석하여 바이오마커를 발굴하고 검증하는 **단일 스트림(Single Stream) 통합 파이프라인**의 설계도입니다. 기존에 파편화되어 있던 11개의 통계/머신러닝 기법들을 유기적으로 결합하여, 하나의 스크립트 실행으로 논문 한 편의 서사가 완성되도록 구축되었습니다.

## 핵심 구동 스크립트
* **`src/comprehensive_biomarker_pipeline.py`**: 모든 분석 및 시각화를 관장하는 메인 엔진입니다. (또는 `src/main_pipeline.py`로 오케스트레이션 수행)
* **결과물 저장소**: 모든 시각화 결과는 `outputs/comprehensive_analysis/` 폴더에 순차적인 번호표(1_, 2_, 3_...)를 달고 생성됩니다.

---

## 파이프라인 워크플로우 (3단계 검증 아키텍처)

### Phase 1: 1차 유효 변수 선택 (Feature Selection)
수만 개의 노이즈 데이터 속에서 타겟(예: 고면역군/저면역군)과 연관된 1차 후보군을 압축합니다.
* **적용 기법**: Pearson Correlation Filter $\rightarrow$ LASSO (L1-Regularization)
* **출력물**: 
  * `1_lasso_coefficients.png`: 페널티를 견디고 살아남은 후보 유전자들의 가중치 시각화
  * `1_correlation_heatmap.png`: 선별된 후보군 간의 다중공선성 및 동반 발현(Co-expression) 네트워크 히트맵

### Phase 2: 진단 모델 구축 및 최종 마커 확정 (Diagnostic Modeling)
1차 후보군을 입력받아 두 가지 머신러닝 모델을 훈련하며, 모델의 가지치기 특성을 활용해 최종 바이오마커를 확정합니다.
* **적용 기법**: Decision Tree, Logistic Regression
* **출력물**:
  * `2_decision_tree.svg`: 모델이 실제 분기(Split)에 사용한 핵심 유전자(최종 바이오마커)의 컷오프(Cut-off) 구조도
  * `2_lr_coefficients_table.png`: 최종 확정된 유전자들의 오즈비(Odds Ratio) 및 회귀 계수 테이블

### Phase 3: 최종 바이오마커 다각도 검증 (Validation Suite)
최종 선별된 유전자들과 구축된 모델이 임상적으로 완벽한지 8가지 엄격한 통계 기법으로 교차 검증합니다.
1. **계층적 군집화 (`3b_clustering.png`)**: 유전자 발현 패턴에 따른 종양 이질성(Tumor Heterogeneity) 분리 확인.
2. **차원 축소 (`3f_pca_tsne.png`)**: PCA 및 t-SNE를 통해 다차원 공간에서 타겟 그룹 간의 군집 분리도 시각화.
3. **모델 스코어 검증 (`3c_model_scores_ttest.png`)**: 모델(LR, DT)이 예측한 확률 점수가 실제 타겟 그룹 간에 유의미한 차이(t-test)를 보이는지 박스플롯 증명.
4. **개별 마커 검증 (`3d_biomarkers_ttest.png`)**: 각 유전자 발현량의 그룹 간 t-test.
5. **생존 분석 (`3e_survival_km.png`)**: 예측 모델 기반 분류가 실제 환자의 생존 기간(OS_days)을 가르는지 Kaplan-Meier 및 Log-rank test 적용.
6. **임상 결정 곡선 및 보정 (`3g_dca_calibration.png`)**: DCA 곡선을 통한 임상적 순이익(Net Benefit) 입증 및 Calibration 정확도 평가.
7. **위상 확인 (`3h_volcano_plot.png`)**: 전체 유전자 볼케이노 플롯에서 선별된 최종 마커의 위치(최우수 유의성)를 하이라이팅.

## 개발 및 실행 환경
* Python 3.8+
* 필수 라이브러리: `pandas`, `numpy`, `scipy`, `scikit-learn`, `lifelines`, `matplotlib`, `seaborn`, `requests`
