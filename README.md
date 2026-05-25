# Medical Data Analysis & Visual Aids (gowun 브랜치)

본 프로젝트(gowun 브랜치)는 다중 오믹스(Multi-omics) 및 임상 병리 데이터를 활용하여 종양 연구와 의료 통계 분석에 필요한 핵심 방법론과 고도화된 시각화 파이프라인을 파이썬(Python)으로 구현한 코드베이스입니다.

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

### 2. 스크립트 실행
`src/` 디렉터리에 번호가 매겨진 11개의 파이썬 분석 스크립트가 있습니다. 터미널에서 실행하고자 하는 스크립트를 구동하세요.
```bash
# 예시: 생존 분석 및 Kaplan-Meier Plot 생성
python src/03_survival_analysis.py

# 예시: 고급 시각화 도구 (노모그램, 화산 플롯 등) 생성
python src/07_advanced_visualization.py
```

### 3. 결과 확인
* 스크립트가 실행되면 결과 이미지(PNG, SVG 등)는 모두 `outputs/` 디렉터리 하위의 카테고리별 폴더에 저장됩니다.
* 각 분석의 목적, 의료 통계적 의미 및 출력된 **그래프를 해석하는 방법**은 `docs/results_explanation.md` 문서를 참고해 주십시오.
* 전체 방법론의 기획 의도는 `docs/methodology_plan.md`에서 확인하실 수 있습니다.

---

## 📖 참고 문헌 및 자료 (References)
이 브랜치의 코드 파이프라인과 시각화 방법론은 다음 자료들을 심층 분석 및 참고하여 개발되었습니다.

* **의료 데이터 분석 필수 방법론 및 시각화 도구 요약 (Google Docs)**
  * 링크: [Methodology Summary Docs](https://docs.google.com/document/d/10X8jf6L8_PoPvcxdMODvtSFAUmv7hRngIzUBJWXtLbU/edit?usp=sharing)
* **참고 연구 논문 아카이브 (Google Drive)**
  * 링크: [Reference Papers Drive](https://drive.google.com/drive/folders/1HLuMjDPXKpNBRTDru3snUya64gqgskQR?usp=sharing)
  * *분석 대상 논문 예시: Cancer Research and Treatment, Experimental & Molecular Medicine, Modern Pathology 등 다수의 중개연구 논문*
