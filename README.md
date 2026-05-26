# Medical Data Analysis & Visual Aids

본 프로젝트는 다중 오믹스(Multi-omics) 및 임상 병리 데이터를 활용하여 종양 연구와 의료 통계 분석에 필요한 핵심 방법론과 고도화된 시각화 파이프라인을 파이썬(Python)으로 구현한 코드베이스입니다.

---

## 📌 주요 특징
* **논문 출판 수준의 시각화**: Box Plot, Kaplan-Meier 생존 곡선, 의사결정나무, 노모그램, 볼케이노 플롯, LASSO Path 등 임상 논문에 필수적인 시각화 도구가 포함되어 있습니다.
* **초보자 친화적 가이드**: 각 그래프가 갖는 의료 통계학적 의미와 그래프 읽는 법이 문서(`docs/results_explanation.md`)로 상세히 정리되어 있습니다.
* **실제 BRCA 데이터 기반**: TCGA 유방암(BRCA) 샘플 데이터를 통해 코드가 바로 구동되는 것을 확인할 수 있습니다.

---

## 🚀 사용 방법 (Usage)

### 1. 환경 설정 및 의존성 설치
본 스크립트들은 Python 3.8 이상 환경에서 작동을 권장합니다. 아래 명령어를 통해 필요한 분석 패키지들을 설치하세요.
```bash
pip install -r requirements.txt
```

### 2. 마스터 파이프라인 실행 (Recommended)
개별 스크립트를 하나씩 실행할 필요 없이, `main_pipeline.py`를 실행하면 1번부터 11번까지의 전체 분석 파이프라인(데이터 로드 -> 통계 분석 -> 피처 셀렉션 -> 머신러닝 예측 -> 최종 논문용 시각화)이 자동 오케스트레이션되어 순차적으로 구동됩니다.
```bash
python src/main_pipeline.py
```

### 3. 개별 스크립트 실행
특정 분석만 필요할 경우 개별 실행도 가능합니다.
```bash
# 예시: 생존 분석 및 Kaplan-Meier Plot 생성
python src/03_survival_analysis.py

# 예시: 차원 축소(PCA/t-SNE) 클러스터링
python src/09_dimensionality_reduction.py
```

### 4. 결과 확인
* 스크립트가 실행되면 결과 이미지(PNG, SVG 등)는 모두 `outputs/` 디렉터리 하위의 카테고리별 폴더에 저장됩니다.
* 전체 파이프라인의 유기적인 흐름과 아키텍처는 `docs/pipeline_architecture.md`에서 확인하세요.
* 각 분석 그래프를 해석하는 방법은 `docs/results_explanation.md` 문서를 참고해 주십시오.

---
