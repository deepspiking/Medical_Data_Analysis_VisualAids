# SPSS Bootstrapping 방법론 조사 결과

> 조사 목적: SPSS Bootstrapping 사용 시 Python 직접 분석 결과와 동일하게 나오는지 검증
> 생성일: 2026-08-17T12:59:39.838636
> 참여 에이전트: Statistician, SPSS_Expert, Reproducer, Researcher

---

## Q1: SPSS Bootstrapping의 알고리즘: simple random vs stratified, resampling 전략

### 🔹 Statistician

"Statistician"으로서, SPSS Statistics의 BOOTSTRAP 프로시저의 알고리즘에 대한 심층 분석을 통해 방법론적 정확성을 검증하고, 제시된 질문들에 대해 전문적인 답변을 드리겠습니다. 이전 에이전트들의 의견에 동의하며, 이를 보완하고 구체화하는 방향으로 진행하겠습니다.

---

## SPSS Bootstrapping 알고리즘 상세 분석

SPSS Statistics의 `BOOTSTRAP` 프로시저는 비모수적 부트스트랩(Non-parametric Bootstrap) 방식을 기반으로 하며, 통계적 추론의 견고성(robustness)을 높이는 데 기여합니다. 각 질문에 대해 상세히 분석하겠습니다.

### (1) 기본 Resampling 방법: Simple Random Sampling vs. Stratified Sampling

SPSS `BOOTSTRAP` 프로시저의 **기본(default) resampling 방법은 복원 추출을 통한 단순 무작위 표본 추출(Simple Random Sampling with replacement)**입니다. 이는 원본 데이터셋의 $N$개 관측치 중에서 $N$개의 관측치를 무작위로 복원 추출하여 하나의 부트스트랩 표본(bootstrap sample)을 생성하는 방식입니다.

그러나 SPSS는 **계층화된 부트스트랩(Stratified Bootstrap) 옵션도 제공**합니다. 이는 특정 변수(예: 성별, 처리군/대조군 등)에 따라 데이터가 여러 계층(strata)으로 나뉘어 있을 때 유용합니다. 계층화된 부트스트랩을 사용하면 각 계층 내에서 독립적으로 복원 추출을 수행하여, 원본 데이터셋의 계층별 비율을 부트스트랩 표본에서도 유지할 수 있습니다. 이는 특히 특정 계층의 표본 크기가 작거나, 계층 간의 분포가 크게 다를 때 추정치의 안정성과 정확성을 높이는 데 중요합니다.

*   **Simple Random Sampling (단순 무작위 표본 추출):** 모든 관측치가 동일한 확률로 선택됩니다.
*   **Stratified Sampling (계층화 표본 추출):** 사용자가 지정한 계층 변수(stratification variable)에 따라 데이터를 나눈 후, 각 계층 내에서 독립적으로 복원 추출을 수행합니다. SPSS에서는 `BOOTSTRAP` 대화 상자의 "Stratification" 옵션을 통해 계층 변수를 지정할 수 있습니다.

### (2) Cox 회귀분석 및 로지스틱 회귀분석에서의 적용 (원데이터 재추출 vs. 잔차 재추출)

SPSS의 `BOOTSTRAP` 프로시저는 Cox 회귀분석(Cox Regression)과 로지스틱 회귀분석(Logistic Regression)을 포함한 대부분의 모델에서 **원데이터 재추출(resampling original data/cases) 방식**을 사용합니다.

*   **원데이터 재추출 (Resampling Original Data/Cases):**
    *   이 방식은 원본 데이터셋의 각 관측치(즉, 종속변수와 모든 독립변수를 포함하는 하나의 행)를 하나의 단위로 보고, 이를 복원 추출하여 부트스트랩 표본을 생성합니다.
    *   **장점:** 모델의 가정(예: 잔차의 분포)에 대한 의존성이 낮고, 공변량(covariates) 간의 관계를 그대로 유지하며, 복잡한 모델이나 비선형 모델에도 일반적으로 적용 가능합니다.
    *   **Cox 회귀분석:** 생존 시간과 사건 발생 여부(event status), 그리고 공변량을 포함하는 전체 관측치를 재추출합니다. 이는 Censoring 패턴과 공변량 구조를 보존하는 데 필수적입니다.
    *   **로지스틱 회귀분석:** 이항 또는 다항 종속변수와 공변량을 포함하는 전체 관측치를 재추출합니다. 이는 이항 분포를 따르는 종속변수의 특성과 공변량 간의 관계를 유지합니다.

*   **잔차 재추출 (Resampling Residuals):**
    *   이 방식은 주로 선형 회귀 모델(Linear Regression)에서 잔차(residuals)가 독립적이고 동일한 분포를 따른다는 가정이 비교적 잘 충족될 때 사용될 수 있습니다. 모델을 적합한 후 잔차를 계산하고, 이 잔차를 재추출하여 예측값에 더하는 방식입니다.
    *   **Cox 회귀분석 및 로지스틱 회귀분석에서는 일반적으로 잔차 재추출을 사용하지 않습니다.**
        *   **로지스틱 회귀분석:** 잔차가 정규 분포를 따르지 않고 이항 분포의 특성을 가지므로, 잔차 재추출은 부적절합니다. 또한, 잔차의 분산이 예측 확률에 따라 달라지는 이분산성(heteroscedasticity)을 가집니다.
        *   **Cox 회귀분석:** Cox 모델은 부분 우도(partial likelihood)를 기반으로 하며, "잔차"의 개념이 선형 회귀와는 다릅니다 (예: Martingale residuals, Schoenfeld residuals 등). 이러한 잔차들은 특정 목적을 위해 사용되지만, 직접 재추출하여 모델을 재구성하는 방식은 통계적으로 복잡하고 일반적인 부트스트랩 방식이 아닙니다. 생존 데이터의 특성상 Censoring이 존재하므로, 잔차 재추출은 Censoring 패턴을 적절히 반영하기 어렵습니다.

따라서 SPSS `BOOTSTRAP` 프로시저는 Cox 및 로지스틱 회귀분석에서 **원데이터 재추출**을 통해 모델 계수의 분포를 추정합니다.

### (3) 1000회 Resampling의 구체적 의미

"1000회 resampling"은 부트스트랩 표본(bootstrap samples)을 생성하는 횟수, 즉 **부트스트랩 반복 횟수($B$)**를 의미합니다.

*   **의미:** SPSS는 원본 데이터셋에서 1000개의 독립적인 부트스트랩 표본을 생성합니다. 각 부트스트랩 표본에 대해 사용자가 지정한 통계량(예: 회귀 계수, 오즈비, 위험비 등)을 계산합니다. 결과적으로 1000개의 통계량 추정치 분포를 얻게 됩니다.
*   **목적:** 이 1000개의 추정치 분포를 사용하여 원래 통계량의 표준 오차(standard error), 편향(bias), 신뢰 구간(confidence interval) 등을 추정합니다.
*   **횟수의 중요성:**
    *   **표준 오차 추정:** 일반적으로 200~500회 정도의 반복으로도 표준 오차 추정치는 비교적 안정적입니다.
    *   **신뢰 구간 추정:** 특히 BCa(Bias-Corrected and Accelerated)와 같은 고급 신뢰 구간 추정 방식의 경우, 1000회 이상의 반복(예: 2000회, 5000회)이 더 정확하고 안정적인 결과를 제공하는 데 권장됩니다. 1000회는 많은 통계 소프트웨어에서 기본값으로 제공되는 일반적인 횟수이며, 대부분의 경우 합리적인 균형을 제공합니다. 반복 횟수가 적으면 부트스트랩 분포가 원본 통계량의 실제 분포를 잘 반영하지 못하여 추정치의 변동성이 커질 수 있습니다.

### (4) Bootstrap Confidence Interval 계산 방식

SPSS `BOOTSTRAP` 프로시저는 여러 가지 부트스트랩 신뢰 구간(Confidence Interval, CI) 계산 방식을 제공하며, 각각의 원리는 다음과 같습니다.

1.  **Percentile Interval (백분위수 신뢰 구간):**
    *   **원리:** 가장 직관적이고 구현하기 쉬운 방법입니다. $B$개의 부트스트랩 통계량 추정치들을 오름차순으로 정렬한 후, 원하는 신뢰 수준(예: 95%)에 해당하는 백분위수를 직접 취하여 신뢰 구간을 형성합니다. 예를 들어, 95% 신뢰 구간의 경우, 정렬된 부트스트랩 분포의 2.5번째 백분위수와 97.5번째 백분위수를 사용합니다.
    *   **장점:** 계산이 간단합니다.
    *   **단점:** 부트스트랩 분포가 대칭적이고 편향되지 않았다는 가정에 의존합니다. 만약 분포가 비대칭적이거나 편향되어 있다면, 이 신뢰 구간은 실제 모수(parameter)를 포함할 확률이 95%가 아닐 수 있습니다.

2.  **Bias-Corrected (BC) Interval (편향 보정 신뢰 구간):**
    *   **원리:** 백분위수 신뢰 구간의 단점 중 하나인 편향(bias)을 보정합니다. 이는 부트스트랩 분포의 중앙값이 원본 표본에서 계산된 통계량과 얼마나 다른지를 측정하는 편향 보정 계수($z_0$)를 계산하여, 백분위수 신뢰 구간의 상한과 하한을 조정합니다.
    *   **장점:** 부트스트랩 분포에 편향이 있을 때 백분위수 신뢰 구간보다 더 정확합니다.
    *   **단점:** 분포의 비대칭성(skewness) 또는 분산의 비정상성(non-constant variance)은 보정하지 못합니다.

3.  **Bias-Corrected and Accelerated (BCa) Interval (편향 보정 및 가속화 신뢰 구간):**
    *   **원리:** 현재까지 개발된 부트스트랩 신뢰 구간 중 가장 정교하고 일반적으로 권장되는 방법입니다. BC 신뢰 구간의 편향 보정($z_0$)뿐만 아니라, 부트스트랩 분포의 비대칭성(skewness) 또는 분산의 비정상성(non-constant variance)을 보정하는 가속화 계수($a$)를 추가로 고려합니다. 이 두 계수를 사용하여 백분위수 신뢰 구간의 상한과 하한을 동적으로 조정합니다.
    *   **장점:** 부트스트랩 분포가 편향되거나 비대칭적일 때 가장 정확한 신뢰 구간을 제공합니다.
    *   **단점:** 계산이 복잡하며, 안정적인 추정을 위해 더 많은 부트스트랩 반복 횟수(예: 1000회 이상)가 필요합니다.

SPSS는 이 세 가지 방법을 모두 제공하며, `BOOTSTRAP` 대화 상자의 "Confidence Intervals" 섹션에서 선택할 수 있습니다. 통계 방법론자로서, 특별한 이유가 없다면 **BCa 신뢰 구간을 우선적으로 고려**할 것을 권장합니다.

### (5) Random Seed 설정 방법

**Random Seed (난수 시드) 설정은 부트스트랩 분석의 재현성(reproducibility)을 보장하는 데 필수적입니다.** 부트스트랩은 무작위 표본 추출을 기반으로 하므로, 난수 생성기의 초기값(seed)을 고정하지 않으면 매번 다른 부트스트랩 표본이 생성되어 결과가 달라질 수 있습니다.

SPSS에서 Random Seed를 설정하는 방법은 다음과 같습니다:

1.  **`SET` 명령어 사용:** 구문(Syntax) 창에서 `SET SEED = [임의의 정수].` 명령어를 실행합니다. 예를 들어, `SET SEED = 2023.`과 같이 설정할 수 있습니다. 이 명령어는 SPSS 세션 동안 난수 생성기의 시드를 고정합니다.
2.  **`BOOTSTRAP` 프로시저 대화 상자:** `BOOTSTRAP` 대화 상자 내에 "Random Seed" 옵션이 있습니다. 여기에 특정 정수 값을 입력하여 해당 부트스트랩 분석에만 적용되는 시드를 설정할 수 있습니다. 이 방법은 특정 분석에 대한 재현성을 보장하는 데 유용합니다.

**권장 사항:** 연구 보고서나 논문 작성 시에는 사용한 Random Seed 값을 명시하여 다른 연구자들이 동일한 결과를 재현할 수 있도록 하는 것이 좋습니다.

---

### 방법론적 검증 및 Python 구현과의 수학적 동일성 판단

SPSS `BOOTSTRAP` 프로시저가 사용하는 알고리즘은 Efron과 Tibshirani의 고전적인 부트스트랩 방법론(Efron & Tibshirani, 1993)에 기반을 두고 있으며, 이는 통계학계에서 널리 인정받는 표준적인 접근 방식입니다.



### 🔹 SPSS_Expert

안녕하십니까, "SPSS 전문가" SPSS_Expert입니다.

"Statistician" 에이전트님의 심층적이고 정확한 분석에 전적으로 동의하며, SPSS Statistics의 `BOOTSTRAP` 프로시저의 알고리즘에 대한 탁월한 설명을 해주셨습니다. 저는 SPSS_Expert로서, Statistician님의 의견을 보완하고, 실제 SPSS 사용 환경에서의 구체적인 적용 방식과 추가적인 전문적 통찰을 더하여 조사 질문에 대한 답변을 드리겠습니다.

SPSS Statistics의 `BOOTSTRAP` 프로시저는 통계적 추론의 견고성(robustness)을 크게 향상시키는 강력한 도구이며, 그 알고리즘은 통계학계에서 널리 인정받는 표준적인 비모수적 부트스트랩(Non-parametric Bootstrap) 방법론을 따릅니다.

---

## SPSS Bootstrapping 알고리즘 상세 분석 (SPSS_Expert 관점)

### (1) 기본 Resampling 방법: Simple Random Sampling vs. Stratified Sampling

Statistician님의 의견에 전적으로 동의합니다.

SPSS `BOOTSTRAP` 프로시저의 **기본(default) resampling 방법은 복원 추출을 통한 단순 무작위 표본 추출(Simple Random Sampling with replacement)**입니다. 이는 원본 데이터셋의 각 관측치(case)가 동일한 확률로 선택되어 부트스트랩 표본을 구성함을 의미합니다.

그러나 SPSS는 **계층화된 부트스트랩(Stratified Bootstrap) 옵션을 명확하게 제공**하여 사용자가 특정 연구 설계에 맞춰 유연하게 적용할 수 있도록 합니다.

*   **SPSS 구현 방식:**
    *   `Analyze` > `Regression` (또는 다른 분석 프로시저) > `Bootstrap...` 대화 상자로 진입합니다.
    *   `Bootstrap` 대화 상자에서 `Perform bootstrapping`을 체크한 후, `Stratification` 탭을 클릭합니다.
    *   여기서 하나 이상의 계층 변수(stratification variable)를 지정할 수 있습니다. 예를 들어, `Gender`나 `Treatment_Group`과 같은 범주형 변수를 계층 변수로 지정하면, SPSS는 각 계층 내에서 독립적으로 복원 추출을 수행하여 원본 데이터의 계층별 비율을 부트스트랩 표본에서도 유지합니다.
*   **활용 시점:** 특정 집단(계층)의 크기가 작거나, 계층 간의 분포가 크게 달라 각 계층의 대표성을 보장해야 할 때 계층화 부트스트랩은 필수적입니다. 이는 특히 임상 시험이나 사회 과학 연구에서 중요한 고려 사항입니다.

### (2) Cox 회귀분석 및 로지스틱 회귀분석에서의 적용 (원데이터 재추출 vs. 잔차 재추출)

Statistician님의 의견에 전적으로 동의하며, SPSS의 `BOOTSTRAP` 프로시저는 Cox 회귀분석(Cox Regression)과 로지스틱 회귀분석(Logistic Regression)을 포함한 모든 모델에서 **원데이터 재추출(resampling original data/cases)** 방식을 사용합니다.

*   **원데이터 재추출 (Resampling Original Data/Cases):**
    *   SPSS는 각 관측치(즉, 종속변수와 모든 독립변수를 포함하는 하나의 완전한 데이터 행)를 하나의 단위로 간주하고, 이를 복원 추출하여 부트스트랩 표본을 생성합니다.
    *   **로지스틱 회귀분석:** 이항 또는 다항 종속변수와 공변량(covariates)의 관계를 분석합니다. 잔차(residuals)는 이항 분포의 특성을 가지며, 이분산성(heteroscedasticity)을 보이므로 선형 회귀에서와 같은 잔차 재추출은 부적절합니다. 원데이터 재추출은 종속변수의 이항적 특성과 공변량 간의 복잡한 관계를 그대로 보존하면서 모델 계수의 분포를 추정합니다.
    *   **Cox 회귀분석:** 생존 시간(survival time), 사건 발생 여부(event status), 그리고 공변량을 포함하는 전체 관측치를 재추출합니다. 생존 데이터는 Censoring(중도절단)이라는 특수한 형태의 데이터가 포함되어 있어, 잔차 재추출 방식으로는 Censoring 패턴과 공변량의 구조를 적절히 반영하기 매우 어렵습니다. Cox 모델의 "잔차" 개념(예: Martingale residuals, Schoenfeld residuals)은 모델 진단에 사용되지만, 이를 재추출하여 모델을 재구성하는 방식은 일반적인 부트스트랩 접근법이 아닙니다.
*   **SPSS의 일관된 접근:** SPSS `BOOTSTRAP` 프로시저는 일반화 선형 모델(Generalized Linear Models, GLM) 계열(로지스틱 회귀 포함), 선형 회귀, Cox 회귀 등 다양한 모델에서 일관되게 원데이터 재추출 방식을 채택합니다. 이는 모델의 가정에 대한 의존성을 최소화하고, 공변량 간의 복잡한 상관관계 및 데이터 구조를 보존하는 가장 견고한 방법입니다. **SPSS는 Cox 회귀나 로지스틱 회귀에서 잔차 재추출 옵션을 제공하지 않습니다.**

### (3) 1000회 Resampling의 구체적 의미

Statistician님의 의견에 전적으로 동의합니다.

"1000회 resampling"은 **부트스트랩 반복 횟수($B$)**를 의미하며, SPSS는 원본 데이터셋에서 1000개의 독립적인 부트스트랩 표본을 생성합니다. 각 표본에 대해 사용자가 요청한 통계량(예: 회귀 계수, 오즈비, 위험비 등)을 계산하고, 이 1000개의 통계량 추정치들로 구성된 분포를 생성합니다.

*   **SPSS 구현 방식:**
    *   `Bootstrap` 대화 상자에서 `Number of samples` 필드에 원하는 반복 횟수를 입력합니다. 기본값은 1000입니다.
*   **횟수의 중요성 (SPSS_Expert 관점):**
    *   **표준 오차(Standard Error) 추정:** 1000회는 표준 오차 추정에는 충분한 횟수로 간주됩니다.
    *   **신뢰 구간(Confidence Interval) 추정:** 특히 BCa(Bias-Corrected and Accelerated)와 같은 고급 신뢰 구간은 부트스트랩 분포의 꼬리 부분에 민감하므로, 1000회보다 더 많은 반복(예: 2000회, 5000회)이 더 안정적이고 정확한 결과를 제공할 수 있습니다. 계산 시간이 허락한다면, 특히 BCa를 사용할 때는 2000회 이상을 고려하는 것이 좋습니다. SPSS는 최대 10,000회까지 설정할 수 있습니다.

### (4) Bootstrap Confidence Interval 계산 방식

Statistician님의 의견에 전적으로 동의합니다. SPSS `BOOTSTRAP` 프로시저는 다음과 같은 세 가지 주요 신뢰 구간 계산 방식을 제공합니다.

1.  **Percentile Interval (백분위수 신뢰 구간):**
    *   **원리:** 부트스트랩 통계량 분포의 직접적인 백분위수를 사용합니다. 95% 신뢰 구간의 경우, 정렬된 1000개의 부트스트랩 통계량 중 2.5번째 백분위수와 97.5번째 백분위수를 경계로 합니다.
    *   **SPSS 선택:** `Bootstrap` 대화 상자의 `Confidence Intervals` 탭에서 `Percentile`을 선택합니다.
    *   **주의사항:** 부트스트랩 분포가 대칭적이고 편향되지 않았을 때 가장 적합합니다. 실제 데이터에서는 이러한 가정이 위배되는 경우가 많아, 다른 방법보다 정확도가 떨어질 수 있습니다.

2.  **Bias-Corrected (BC) Interval (편향 보정 신뢰 구간):**
    *   **원리:** 부트스트랩 분포의 편향(bias)을 보정하기 위해 편향 보정 계수($z_0$)를 계산하여 백분위수 신뢰 구간을 조정합니다.
    *   **SPSS 선택:** `Bootstrap` 대화 상자의 `Confidence Intervals` 탭에서 `Bias-corrected`를 선택합니다.
    *   **장점:** 편향이 있는 경우 Percentile Interval보다 개선된 결과를 제공합니다.

3.  **Bias-Corrected and Accelerated (BCa) Interval (편향 보정 및 가속화 신뢰 구간):**
    *   **원리:** 편향 보정($z_0$)뿐만 아니라, 부트스트랩 분포의 비대칭성(skewness)과 분산의 비정상성(non-constant variance)을 보정하는 가속화 계수($a$)를 추가로 고려하여 신뢰 구간을 동적으로 조정합니다.
    *   **SPSS 선택:** `Bootstrap` 대화 상자의 `Confidence Intervals` 탭에서 `BCa`를 선택합니다.
    *   **SPSS_Expert 권장:** BCa는 현재까지 개발된 부트스트랩 신뢰 구간 중 가장 정교하고 일반적으로 가장 정확한 방법으로 널리 인정받고 있습니다. 특히 부트스트랩 분포가 비대칭적이거나 편향되어 있을 때 강력한 성능을 발휘합니다. 계산 비용이 더 들고 더 많은 부트스트랩 반복 횟수(예: 1000회 이상)가 필요하지만, **특별한 이유가 없다면 BCa 신뢰 구간을 우선적으로 고려하고 보고할 것을 강력히 권장합니다.**

### (5) Random Seed 설정 방법

Statistician님의 의견에 전적으로 동의하며, Random Seed (난수 시드) 설정은 부트스트랩 분석의 **재현성(reproducibility)**을 보장하는 데 절대적으로 필수적입니다.

SPSS에서 Random Seed를 설정하는 방법은 다음과 같습니다.

1.  **`SET SEED` 명령어 사용 (SPSS 세션 전체에 적용):**
    *   구문(Syntax) 창을 열고 다음 명령어를 실행합니다:
        ```sps
        SET SEED = 2023.
        ```
    *   `2023` 대신 어떤 양의 정수든 사용할 수 있습니다. 이 명령어는 현재 SPSS 세션 동안 난수 생성기의 초기값을 고정합니다. 따라서 이 명령어를 실행한 후 수행되는 모든 무작위 샘플링(부트스트랩 포함)은 동일한 시퀀스의 난수를 사용하게 되어 결과를 재현할 수 있습니다.
    *   **활용 팁:** 분석 스크립트의 시작 부분에 `SET SEED` 명령어를 포함시켜, 해당 스크립트를 실행할 때마다 동일한 부트스트랩 결과를 얻을 수 있도록 하는 것이 좋습니다.

2.  **`BOOTSTRAP` 프로시저 대화 상자 내 "Random Seed" 옵션 사용 (특정 분석에만 적용):**
    *   `Analyze` > `Regression` (또는 다른 분석 프로시저) > `Bootstrap...` 대화 상자로 진입합니다.
    *   `Bootstrap` 대화 상자에서 `Perform bootstrapping`을 체크한 후, `Random Seed` 필드에 특정 정수 값을 입력합니다.
    *   **활용 팁:** 이 방법은 특정 부트스트랩 분석에만 시드를 적용하고자 할 때 유용합니다. `SET SEED` 명령어가 전역적으로 적용되는 반면, 이 옵션은 해당 `BOOTSTRAP` 프로시저에만 영향을 미칩니다.

**SPSS_Expert 권장 사항:**
연구 보고서나 논문 작성 시에는 사용한 Random Seed 값을 명시하여, 다른 연구자들이 동일한 데이터와 설정으로 분석을 수행했을 때 동일한 부트스트랩 결과를 재현할 수 있도록 하는 것이 통계적 투명성과 재현성을 높이는 데 매우 중요합니다.

---

### 결론

SPSS Statistics의 `BOOTSTRAP` 프로시저는 Efron과 Tibshirani의 고전적인 부트스트랩 방법론을 충실히 구현하고 있으며, 단순 무작위 표본 추출을 기본으로 하되 계층화 옵션을 제공하여 다양한 연구 설계에 대응할 수 있습니다. Cox 회귀 및 로지스틱 회귀와 같은 복잡한 모델에서는 원데이터 재추출 방식을 사용하여 모델의 가정을 최소화하고 데이터의 본질적인 구조를 보존합니다. 1000회 반복은 일반적인 기본값이지만, 특히 BCa 신뢰 구간을 사용할 때는 더 많은 반복 횟수를 고려하는 것이 좋습니다. 마지막으로, `SET SEED` 명령어나 대화 상자 옵션을 통해 난수 시드를 설정하는 것은 분석의 재현성을 보장하는 핵심적인 단계입니다.

이러한 기능들은 SPSS 사용자가 통계적 추론의 견고성을 높이고, 모수적 가정에 대한 의존도를 줄이며, 복잡한 통계량의 표준 오차와 신뢰 구간을 효과적으로 추정할 수 있도록 지원합니다.

### 🔹 Reproducer

안녕하십니까, "재현성 검증자(Reproducibility Verifier)" 역할을 맡은 AI 에이전트 Reproducer입니다.

Statistician님과 SPSS_Expert님의 SPSS `BOOTSTRAP` 프로시저 알고리즘에 대한 심층적이고 전문적인 분석에 깊이 공감합니다. 두 분의 의견은 SPSS의 부트스트랩 작동 방식을 명확히 이해하는 데 큰 도움이 되었으며, 이를 바탕으로 SPSS 분석 결과를 Python에서 재현하기 위한 구체적이고 실행 가능한 방법을 설계하겠습니다.

재현성(Reproducibility)은 과학적 연구의 신뢰성을 보장하는 핵심 요소입니다. 특히 무작위 추출(random sampling)을 기반으로 하는 부트스트랩(bootstrap)과 같은 방법론에서는 난수 시드(random seed) 고정, 동일한 재표본 추출 전략, 통계 모델의 일치 등 여러 조건을 철저히 통제해야만 동일한 결과를 얻을 수 있습니다.

---

## SPSS Bootstrapping 결과 Python 재현을 위한 구체적 방법 설계

SPSS `BOOTSTRAP` 프로시저의 결과를 Python에서 재현하기 위한 절차는 다음과 같은 핵심 원칙을 따릅니다.

1.  **동일 데이터셋 사용 (Identical Dataset):** SPSS와 Python에서 분석에 사용되는 원본 데이터셋은 한 글자도 틀리지 않고 동일해야 합니다.
2.  **난수 시드 고정 (Fixed Random Seed):** 무작위 추출의 결과를 통제하기 위해 SPSS와 Python 모두에서 동일한 난수 시드를 설정해야 합니다.
3.  **동일한 재표본 추출 전략 (Identical Resampling Strategy):** 단순 무작위 복원 추출(Simple Random Sampling with replacement)인지, 계층화된 복원 추출(Stratified Sampling with replacement)인지 정확히 일치시켜야 합니다.
4.  **동일한 재표본 추출 단위 (Identical Resampling Unit):** SPSS가 원데이터(original data/cases)를 재추출하는 방식을 Python에서도 그대로 따라야 합니다.
5.  **동일한 반복 횟수 (Identical Number of Iterations):** 부트스트랩 표본 생성 횟수($B$)를 일치시켜야 합니다.
6.  **동일한 통계 모델 (Identical Statistical Model):** Cox 회귀분석 또는 로지스틱 회귀분석의 모델 사양(종속변수, 독립변수, 공변량 등)이 SPSS와 Python에서 완전히 동일해야 합니다.
7.  **동일한 신뢰 구간 계산 방식 (Identical Confidence Interval Calculation Method):** Percentile, BC, BCa 중 SPSS에서 선택한 방식과 동일한 알고리즘으로 Python에서 신뢰 구간을 계산해야 합니다.

아래에서는 조사 질문에 대한 답변을 바탕으로 각 재현 조건을 점검하고 Python 코드로 검증 가능한 절차를 제안합니다.

---

### (1) 기본 Resampling 방법: Simple Random Sampling vs. Stratified Sampling

Statistician님과 SPSS_Expert님의 의견에 전적으로 동의합니다. SPSS `BOOTSTRAP` 프로시저의 기본 재표본 추출 방법은 **복원 추출을 통한 단순 무작위 표본 추출(Simple Random Sampling with replacement)**입니다. 계층화된 부트스트랩(Stratified Bootstrap) 옵션도 제공됩니다.

**재현성 검증을 위한 접근:**

*   **SPSS 설정 확인:** `BOOTSTRAP` 대화 상자에서 `Simple random sampling`이 선택되었는지, 또는 `Stratification` 탭에서 특정 계층 변수가 지정되었는지 확인합니다.
*   **Python 구현:**
    *   **단순 무작위 표본 추출:** `numpy.random.choice` 함수를 사용하여 원본 데이터셋의 인덱스(index)를 복원 추출합니다.
    *   **계층화된 표본 추출:** 계층 변수(stratification variable)의 각 수준(level)별로 데이터를 분할한 후, 각 계층 내에서 독립적으로 `numpy.random.choice`를 사용하여 인덱스를 복원 추출하고, 이를 다시 결합하여 하나의 부트스트랩 표본을 만듭니다.

**Python 코드 예시 (Pseudo-code):**

```python
import numpy as np
import pandas as pd

# 1. 데이터 로드 (SPSS에서 사용한 동일한 데이터셋)
# 예: df = pd.read_csv('your_data.csv')

# 2. 난수 시드 고정 (SPSS와 동일한 시드 사용)
SPSS_SEED = 2023 # SPSS에서 SET SEED = 2023. 으로 설정했다고 가정
np.random.seed(SPSS_SEED)

# 원본 데이터셋의 인덱스
original_indices = np.arange(len(df))

# 부트스트랩 표본 생성 함수
def generate_bootstrap_sample(data, stratification_variable=None):
    if stratification_variable is None:
        # 단순 무작위 표본 추출 (Simple Random Sampling with replacement)
        sample_indices = np.random.choice(original_indices, size=len(data), replace=True)
        bootstrap_df = data.iloc[sample_indices]
    else:
        # 계층화된 표본 추출 (Stratified Sampling with replacement)
        bootstrap_parts = []
        for stratum_name, stratum_df in data.groupby(stratum_variable):
            stratum_indices = stratum_df.index.values # 계층별 원본 인덱스
            sample_stratum_indices = np.random.choice(stratum_indices, size=len(stratum_df), replace=True)
            bootstrap_parts.append(data.loc[sample_stratum_indices])
        bootstrap_df = pd.concat(bootstrap_parts)
        # 계층화된 부트스트랩의 경우, 원본 데이터의 순서가 아닌 계층별로 재추출되므로,
        # 인덱스를 재설정하여 분석의 일관성을 유지할 수 있습니다.
        bootstrap_df = bootstrap_df.reset_index(drop=True)
    return bootstrap_df

# 사용 예시:
# bootstrap_sample_df = generate_bootstrap_sample(df) # 단순 무작위
# bootstrap_sample_df_stratified = generate_bootstrap_sample(df, stratification_variable='Gender') # 계층화
```

### (2) Cox 회귀분석 및 로지스틱 회귀분석에서의 적용 (원데이터 재추출 vs. 잔차 재추출)

Statistician님과 SPSS_Expert님의 의견에 전적으로 동의합니다. SPSS `BOOTSTRAP` 프로시저는 Cox 회귀분석과 로지스틱 회귀분석을 포함한 대부분의 모델에서 **원데이터 재추출(resampling original data/cases) 방식**을 사용합니다. 잔차 재추출(resampling residuals)은 이러한 모델에 적합하지 않으며, SPSS에서도 해당 옵션을 제공하지 않습니다.

**재현성 검증을 위한 접근:**

*   **SPSS 설정 확인:** SPSS는 기본적으로 원데이터 재추출을 사용하므로, 별도로 설정할 필요는 없지만 이 점을 명확히 인지해야 합니다.
*   **Python 구현:** 위 `generate_bootstrap_sample` 함수를 사용하여 생성된 부트스트랩 데이터셋(`bootstrap_df`)에 대해 Cox 회귀분석 또는 로지스틱 회귀분석 모델을 적합(fit)합니다. 각 부트스트랩 표본은 원본 데이터셋의 완전한 행(row)을 복원 추출하여 구성됩니다.

**Python 코드 예시 (Pseudo-code):**

```python
# 로지스틱 회귀분석 (Logistic Regression) 예시
import statsmodels.api as sm

def fit_logistic_model(data, dependent_var, independent_vars):
    X = data[independent_vars]
    y = data[dependent_var]
    # 상수항 추가 (SPSS는 기본적으로 상수항을 포함)
    X = sm.add_constant(X, prepend=False)
    model = sm.Logit(y, X)
    results = model.fit(disp=0)

---

## Q2: SPSS vs Python 재현성: 동일 결과가 나오는 조건

### 🔹 Statistician

"Statistician"으로서, SPSS와 Python 간의 Bootstrapping 결과 재현성 문제에 대해 통계 방법론적 관점에서 심층적으로 분석하고 답변드리겠습니다.

Bootstrapping은 데이터의 기저 분포(underlying distribution)에 대한 가정을 최소화하면서 통계량의 표본 분포(sampling distribution)를 추정하는 강력한 비모수적(non-parametric) 재표본추출(resampling) 기법입니다. 그 핵심 원리는 관측된 표본(observed sample)을 모집단(population)의 대리(proxy)로 간주하고, 이 표본으로부터 반복적으로 재표본을 추출하여 관심 있는 통계량의 분포를 시뮬레이션하는 것입니다.

## 조사 질문에 대한 답변

### (1) 어떤 요소가 결과 차이를 만들 수 있는지 (난수 생성기 차이, resampling 단위, CI 방식, 알고리즘 구현 차이)

SPSS Bootstrapping 결과가 Python(예: lifelines, statsmodels, scikit-learn)으로 직접 계산한 결과와 '동일하게' 나오지 않는 주된 원인은 다음과 같습니다.

1.  **난수 생성기(Random Number Generator, RNG)의 차이:**
    *   **원리:** Bootstrapping은 무작위 재표본추출(random resampling)에 기반하므로, 난수 생성기의 역할이 결정적입니다. SPSS와 Python(주로 NumPy 기반)은 서로 다른 난수 생성기 알고리즘을 사용하거나, 동일한 알고리즘(예: Mersenne Twister)이라도 내부 구현 방식, 초기 시드(seed) 설정 방식, 그리고 난수 시퀀스(sequence)를 생성하는 세부 로직이 다를 수 있습니다.
    *   **영향:** 이로 인해 각 재표본추출 단계에서 추출되는 데이터 포인트의 구성이 달라지고, 결과적으로 수천 번의 반복(replications)을 통해 얻어지는 부트스트랩 분포(bootstrap distribution) 자체가 미묘하게 달라집니다. 이는 최종 추정치(estimate)와 신뢰구간(confidence interval)의 차이로 이어집니다.

2.  **리샘플링 단위(Resampling Unit) 및 전략의 차이:**
    *   **원리:** 대부분의 경우 개별 관측치(individual observation)를 리샘플링 단위로 사용하지만, 복잡한 데이터 구조(예: 시계열 데이터, 계층적 데이터, 종속 관측치)에서는 블록 부트스트랩(block bootstrap), 잔차 부트스트랩(residual bootstrap), 또는 모형 기반 부트스트랩(model-based bootstrap)과 같은 특수한 리샘플링 전략이 필요할 수 있습니다.
    *   **영향:** SPSS와 Python 라이브러리(예: `statsmodels`의 회귀 모델 부트스트랩)가 특정 모델이나 데이터 유형에 대해 기본적으로 채택하는 리샘플링 전략이 다를 수 있습니다. 예를 들어, 회귀 분석에서 잔차를 부트스트랩할지, 아니면 전체 관측치 쌍(X, Y)을 부트스트랩할지에 따라 결과가 달라집니다.

3.  **신뢰구간(Confidence Interval, CI) 추정 방식의 차이:**
    *   **원리:** 부트스트랩 신뢰구간은 여러 가지 방법으로 추정될 수 있습니다.
        *   **Percentile Method (백분위수 방법):** 부트스트랩 분포의 특정 백분위수(예: 2.5%와 97.5%)를 직접 사용합니다. 가장 간단하지만, 부트스트랩 분포가 편향되거나 비대칭일 때 정확도가 떨어질 수 있습니다.
        *   **Bias-Corrected and Accelerated (BCa) Method (편향 수정 가속 방법):** 부트스트랩 분포의 편향(bias)과 왜도(skewness)를 보정하여 더 정확한 신뢰구간을 제공합니다. 계산 과정이 복잡하며, 가속도(acceleration) 및 편향 보정 계수(bias correction factor) 추정이 필요합니다.
        *   **Basic Bootstrap (기본 부트스트랩):** 원본 표본 통계량과 부트스트랩 통계량의 차이를 이용합니다.
        *   **Studentized Bootstrap (스튜던트화 부트스트랩):** 부트스트랩 분포의 표준 오차를 사용하여 t-분포와 유사하게 신뢰구간

### 🔹 SPSS_Expert

안녕하십니까, SPSS_Expert입니다. IBM SPSS Statistics의 BOOTSTRAP 명령어, PROCESS macro, 그리고 Cox/로지스틱 회귀에서의 Bootstrap 구현 방식에 대한 깊이 있는 이해를 바탕으로, SPSS와 Python 간의 Bootstrapping 결과 재현성 문제에 대해 전문적으로 분석하고 답변드리겠습니다.

"Statistician" 에이전트님의 의견은 통계 방법론적 관점에서 매우 정확하며, 저는 그 의견에 전적으로 동의하며 SPSS의 구체적인 작동 방식과 Python과의 비교를 통해 더욱 심층적인 분석을 제공하고자 합니다.

---

## 조사 질문에 대한 답변: SPSS vs Python Bootstrapping 재현성

SPSS Bootstrapping 결과가 Python(예: lifelines, statsmodels, scikit-learn)으로 직접 계산한 결과와 '동일하게' 나오기 위한 조건은 매우 까다로우며, 현실적으로 '완전히 동일한 숫자'를 얻는 것은 거의 불가능에 가깝습니다. 그 이유는 다음과 같습니다.

### (1) 어떤 요소가 결과 차이를 만들 수 있는지 (난수 생성기 차이, resampling 단위, CI 방식, 알고리즘 구현 차이)

SPSS와 Python 간의 Bootstrapping 결과 차이를 유발하는 주요 요소들은 다음과 같습니다.

1.  **난수 생성기(Random Number Generator, RNG)의 차이:**
    *   **SPSS의 난수 생성:** IBM SPSS Statistics는 자체적인 난수 생성기를 사용합니다. SPSS Statistics 17 버전부터 기본 난수 생성기는 Mersenne Twister 알고리즘을 채택하고 있습니다. 그러나 동일한 Mersenne Twister 알고리즘을 사용하더라도, 각 소프트웨어의 내부 구현 방식, 초기 시드(seed)를 처리하는 방식, 그리고 난수 시퀀스를 생성하고 사용하는 세부 로직이 다를 수 있습니다. 예를 들어, SPSS는 `SET RNG = MT` 명령어를 통해 Mersenne Twister를 명시적으로 설정할 수 있으며, `SET SEED = <value>` 또는 `BOOTSTRAP /SEED = <value>`를 통해 시드를 설정합니다.
    *   **Python의 난수 생성:** Python의 통계 라이브러리(예: NumPy, SciPy)는 주로 `numpy.random` 모듈을 통해 난수를 생성합니다. `numpy.random` 또한 기본적으로 Mersenne Twister를 사용하지만, 그 구현 방식과 내부 상태 관리, 그리고 난수를 소비하는 방식이 SPSS와 다릅니다. `np.random.seed(<value>)` 또는 `np.random.default_rng(<value>)`를 사용하여 시드를 설정합니다.
    *   **영향:** 이처럼 서로 다른 RNG 구현은 각 재표본추출(resampling) 단계에서 추출되는 데이터 포인트의 순서와 구성이 달라지게 만듭니다. 1000회 반복(replications) 동안 생성되는 부트스트랩 표본(bootstrap samples)의 내용이 매번 달라지므로, 최종적으로 계산되는 통계량의 부트스트랩 분포(bootstrap distribution)가 미묘하게 달라지고, 이는 추정치와 신뢰구간의 차이로 이어집니다.

2.  **리샘플링 단위(Resampling Unit) 및 전략의 차이:**
    *   **SPSS의 리샘플링:**
        *   **`BOOTSTRAP` 명령어:** SPSS의 일반 `BOOTSTRAP` 명령어는 기본적으로 **개별 케이스(individual cases)**를 리샘플링 단위로 사용하여 단순 무작위 복원 추출(simple random sampling with replacement)을 수행합니다. 이는 대부분의 통계량(평균, 회귀계수 등)에 적용됩니다.
        *   **`PROCESS` macro (by Hayes):** PROCESS macro는 회귀 기반의 매개(mediation), 조절(moderation) 분석을 위한 강력한 도구이며, 자체적으로 부트스트랩 기능을 내장하고 있습니다. PROCESS macro 또한 기본적으로 개별 케이스를 리샘플링 단위로 사용합니다.
        *   **Cox/로지스틱 회귀에서의 Bootstrap:** `REGRESSION` (로지스틱 회귀) 또는 `COXREG` (Cox 회귀) 명령어 내의 `BOOTSTRAP` 하위 명령은 해당 모델의 계수에 대한 부트스트랩 신뢰구간을 제공하며, 이 역시 기본적으로 **개별 관측치(cases)**를 리샘플링합니다.
    *   **Python 라이브러리의 리샘플링:** Python 라이브러리(예: `statsmodels`의 회귀 모델, `lifelines`의 생존 분석)는 특정 모델이나 데이터 유형에 따라 다양한 리샘플링 전략을 채택할 수 있습니다. 예를 들어, 회귀 분석에서 잔차를 부트스트랩하는 잔차 부트스트랩(residual bootstrap)을 사용하거나, 관측치 쌍(X, Y)을 부트스트랩하는 쌍 부트스트랩(pairs bootstrap)을 사용할 수 있습니다. `lifelines`와 같은 생존 분석 라이브러리는 Censoring이 있는 데이터를 처리하기 위해 특정 리샘플링 로직을 가질 수 있습니다.
    *   **영향:** 리샘플링 단위나 전략이 다르면, 각 부트스트랩 표본의 구성이 근본적으로 달라지므로, 결과 통계량의 분포 또한 달라질 수밖에 없습니다.

3.  **신뢰구간(Confidence Interval, CI) 추정 방식의 차이:**
    *   **SPSS의 CI 방식:** SPSS `BOOTSTRAP` 명령어는 다양한 신뢰구간 추정 방식을 제공합니다.
        *   **Percentile (백분위수 방법):** 부트스트랩 분포의 특정 백분위수(예: 2.5%와 97.5%)를 직접 사용합니다.
        *   **BCa (Bias-Corrected and Accelerated Method, 편향 수정 가속 방법):** 부트스트랩 분포의 편향과 왜도(skewness)를 보정하여 더 정확한 신뢰구간을 제공합니다. 계산 과정이 복잡하며, 가속도(acceleration) 및 편향 보정 계수(bias correction factor) 추정이 필요합니다.
        *   **Basic (기본 부트스트랩):** 원본 표본 통계량과 부트스트랩 통계량의 차이를 이용합니다.
        *   **Standard Normal (표준 정규 근사):** 부트스트랩 분포의 평균과 표준 오차를 사용하여 정규 분포를 가정하고 신뢰구간을 계산합니다.
    *   **Python 라이브러리의 CI 방식:** Python 라이브러리 또한 Percentile, BCa 등을 포함한 다양한 CI 추정 방식을 제공합니다.
    *   **영향:** 동일한 부트스트랩 분포가 생성되었다 하더라도, 어떤 CI 추정 방식을 선택하느냐에 따라 신뢰구간의 하한과 상한이 달라집니다. 특히 BCa 방법은 계산 과정이 복잡하여, 구현 세부 사항에 따라 미세한 차이가 발생할 수 있습니다.

4.  **알고리즘 구현 및 수치 계산의 차이:**
    *   **기저 통계 모델 구현:** 부트스트랩은 재표본추출된 데이터에 대해 특정 통계 모델(예: 로지스틱 회귀, Cox 회귀)을 반복적으로 적합하는 과정입니다. SPSS와 Python 라이브러리(예: `statsmodels`, `scikit-learn`, `lifelines`)는 동일한 통계 모델이라도 최대 우도 추정(Maximum Likelihood Estimation, MLE)을 위한 최적화 알고리즘, 수렴 기준, 초기값 설정, 수치 정밀도(numerical precision) 등에서 미묘한 차이를 가질 수 있습니다.
    *   **결측치 처리:** 결측치(missing values)를 처리하는 방식(예: listwise deletion, pairwise deletion, imputation)이 다를 경우, 각 부트스트랩 표본의 유효한 관측치 수가 달라져 결과에 영향을 미칠 수 있습니다.
    *   **엣지 케이스 처리:** 로지스틱 회귀에서의 완전 분리(perfect separation)와 같은 엣지 케이스를 처리하는 방식도 다를 수 있습니다.
    *   **영향:** 이러한 미세한 구현 차이들은 1000회 이상의 반복 과정에서 누적되어 최종 부트스트랩 추정치와 신뢰구간에 차이를 발생시킬 수 있습니다.

### (2) 정확히 동일한 숫자를 기대하는 것이 통계적으로 올바른 기대인지

**아닙니다. 정확히 동일한 숫자를 기대하는 것은 통계적으로 올바른 기대가 아닙니다.**

부트스트랩은 본질적으로 확률적(stochastic) 재표본추출 과정에 기반합니다. 위에서 언급된 난수 생성기, 구현 방식, 수치 정밀도 등의 차이로 인해 두 개의 독립적인 소프트웨어 환경에서 '완전히 동일한' 부트스트랩 결과를 얻는 것은 현실적으로 불가능합니다.

통계 분석의 목표는 모집단에 대한

---

## Q3: 의학 논문에서의 SPSS Bootstrap 표준 사용법 및 검증 절차

### 🔹 Statistician

안녕하십니까, "Statistician"입니다. 통계 방법론자로서 의학 논문에서의 SPSS Bootstrap 표준 사용법 및 검증 절차에 대한 조사 질문에 전문적인 답변을 드리겠습니다. Bootstrapping의 통계적 원리, resampling 전략, confidence interval 추정 방식(percentile, BCa), bias correction 등 방법론적 정확성을 검증하며, Python 구현과의 수학적 동일성 판단을 염두에 두고 답변하겠습니다.

---

### **Bootstrapping의 통계적 원리 및 역할 (개요)**

Bootstrapping은 모분포에 대한 가정을 최소화하면서 표본 통계량의 표본 분포(sampling distribution)를 추정하는 강력한 비모수적(non-parametric) 재표본추출(resampling) 방법입니다. 원본 데이터에서 복원 추출(sampling with replacement)을 통해 다수의 가상 표본(bootstrap samples)을 생성하고, 각 가상 표본에서 관심 통계량(예: 회귀 계수, C-index, 오즈비 등)을 계산합니다. 이 통계량들의 분포를 통해 원본 통계량의 표준 오차(standard error), 신뢰 구간(confidence interval)을 추정하고, 편향(bias)을 보정할 수 있습니다. SPSS를 비롯한 대부분의 통계 소프트웨어에서 구현된 Bootstrap은 이러한 원리를 따르며, Python의 `scikit-learn`, `lifelines`, `statsmodels` 등에서도 동일한 수학적 원리에 기반한 구현이 가능합니다.

---

### **조사 질문에 대한 전문적 답변**

#### **(1) TRIPOD 지침이나 STROBE에서 bootstrap 내부 검증을 어떻게 요구하는지**

*   **TRIPOD (Transparent Reporting of a Multivariable Prediction Model for Individual Prognosis or Diagnosis) 지침:**
    *   TRIPOD 지침은 예후 예측 또는 진단 모델의 개발 및 검증 연구 보고의 투명성을 높이기 위해 고안되었습니다. 이 지침은 특정 통계 방법을 *명시적으로 요구*하기보다는, 모델의 **내부 검증(internal validation)**의 중요성을 매우 강조합니다.
    *   내부 검증의 목적은 모델이 개발 데이터에 과적합(overfitting)되었는지 평가하고, 동일한 모집단에서 추출된 새로운 데이터에 얼마나 잘 일반화될 수 있는지 추정하는 것입니다.
    *   TRIPOD는 내부 검증 방법으로 **재표본추출(resampling) 방법**을 명시적으로 언급하며, 그 예시로 **부트스트랩(bootstrap)**과 교차 검증(cross-validation)을 제시합니다. 특히, 부트스트랩은 모델의 성능 지표(예: C-index, R-squared)에 대한 낙관적 편향(optimism bias)을 추정하고 보정하는 데 효과적인 방법으로 권장됩니다.
    *   **요구 사항:** TRIPOD는 모델 개발 시 내부 검증 방법을 명확히 기술하고, 해당 방법을 통해 얻은 모델 성능 지표(예: 판별력(discrimination), 보정(calibration))와 그 신뢰 구간을 보고하도록 요구합니다. 부트스트랩을 사용했다면, 부트스트랩 반복 횟수(number of resamples), 사용된 신뢰 구간 추정 방식(예: percentile, BCa) 등을 상세히 기술해야 합니다.

*   **STROBE (STrengthening the Reporting of OBservational studies in Epidemiology) 지침:**
    *   STROBE 지침은 관찰 연구(observational studies)의 보고 품질을 향상시키기 위한 체크리스트입니다. STROBE는 주로 연구 설계, 데이터 수집, 통계 분석 및 결과 보고의 투명성과 완전성에 초점을 맞춥니다.
    *   TRIPOD와 달리, STROBE는 특정 통계 모델의 **내부 검증**에 대한 상세한 요구 사항을 직접적으로 다루지는 않습니다. STROBE는 주로 연구의 방법론적 엄격성과 편향(bias) 제어에 중점을 두지만, 예측 모델 자체의 과적합 문제를 직접적으로 다루는 지침은 아닙니다.
    *   **관련성:** 그러나 STROBE의 일반적인 원칙, 즉 통계 분석 방법의 투명한 보고와 결과의 신뢰성 강조는 부트스트랩과 같은 견고한 통계 방법을 사용하여 분석 결과의 불확실성을 정량화하고 보고하는 것을 지지합니다. 예를 들어, 회귀 분석 계수나 효과 크기(effect size)의 신뢰 구간을 추정할 때 부트스트랩을 사용했다면, 이는 STROBE의 "통계 분석" 섹션에서 명확히 보고되어야 합니다.

**결론:** TRIPOD는 예측 모델의 내부 검증을 위해 부트스트랩을 강력히 권장하며, 그 사용법과 결과를 상세히 보고하도록 요구합니다. STROBE는 직접적인 요구 사항은 없으나, 통계 분석의 투명성과 신뢰성 확보라는 광범위한 맥락에서 부트스트랩 사용을 지지합니다.

#### **(2) Harrell의 optimism correction과 bootstrap의 관계**

*   **낙관적 편향(Optimism Bias):** 예측 모델을 개발할 때, 모델은 훈련 데이터(training data)의 특이한 패턴까지 학습하여 실제로는 존재하지 않는 관계를 포착할 수 있습니다. 이로 인해 훈련 데이터에서 평가된 모델 성능은 새로운 데이터에서의 실제 성능보다 과대평가되는 경향이 있는데, 이를 낙관적 편향이라고 합니다.
*   **Harrell의 Optimism Correction:** Frank Harrell 교수는 예측 모델의 성능 지표(예: C-index, R-squared, Brier score)에 내재된 낙관적 편향을 정량화하고 보정하기 위해 부트스트랩 기반의 방법을 제안했습니다. 이 방법은 모델의 일반화 성능(generalizability)을 보다 현실적으로 추정하는 데 필수적입니다.
*   **작동 방식 (Bootstrapping을 이용한 Optimism Correction 절차):**
    1.  **원본 모델 성능 평가:** 원본 데이터셋($D_{orig}$)을 사용하여 예측 모델을 개발하고, 이 모델의 성능 지표($P_{orig}$)를 계산합니다. (예: C-index for survival models).
    2.  **부트스트랩 반복:** $B$번의 반복(예: 100~200회)을 수행합니다. 각 반복 $b$에 대해:
        *   **부트스트랩 표본 생성:** 원본 데이터셋 $D_{orig}$에서 복원 추출하여 부트스트랩 표본 $D_b^*$를 생성합니다.
        *   **모델 개발 및 성능 평가 (부트스트랩 표본 내):** $D_b^*$를 사용하여 새로운 모델을 개발하고, 이 모델의 성능을 $D_b^*$ 내에서 평가합니다 ($P_b^{apparent}$).
        *   **모델 성능 평가 (원본 데이터셋에 적용):** $D_b^*$에서 개발된 모델을 원본 데이터셋 $D_{orig}$에 적용하여 성능을 평가합니다 ($P_b^{test}$).
        *   **낙관적 편향 계산:** 각 반복에서 낙관적 편향($O_b$)은 $P_b^{apparent} - P_b^{test}$로 계산됩니다. 이는 부트스트랩 모델이 자체 데이터에서 얻은 성능과 원본 데이터에서 얻은 성능 간의 차이입니다.
    3.  **평균 낙관적 편향 계산:** 모든 $B$번의 반복에서 계산된 낙관적 편향($O_b$)의 평균($\bar{O}$)을 구합니다.
    4.  **편향 보정된 성능:** 원본 모델의 성능 지표($P_{orig}$)에서 평균 낙관적 편향($\bar{O}$)을 빼서 편향 보정된(optimism-corrected) 모델 성능을 얻습니다: $P_{corrected} = P_{orig} - \bar{O}$.

*   **SPSS 구현:** SPSS는 직접적으로 Harrell의 optimism correction 기능을 제공하지는 않습니다. 그러나 R의 `rms` 패키지나 Python의 `lifelines` (예: `ConcordanceStatistic.bootstrap_score`) 등에서는 이러한 기능을 구현하여 사용할 수 있습니다. SPSS 사용자는 수동으로 부트스트랩을 반복하고 각 단계의 결과를 기록하여 계산하거나, SPSS의 프로그래밍 기능(Python, R 확장)을 활용하여 구현할 수 있습니다.

#### **(3) SPSS bootstrap을 10-fold cross-validation과 병행할 때의 표준 보고 형식**

부트스트랩과 k-fold 교차 검증(cross-validation)은 모두 내부 검증을 위한 재표본추출 방법이지만, 주로 다른 목적이나 다른 단계에서 사용될 수 있습니다. 이 두 방법을 병행하는 것은 드물지만, 특정 상황에서는 상호 보완적으로 활용될 수 있습니다.

*   **일반적인 사용 시나리오:**
    *   **교차 검증:** 주로 모델의 하이퍼파라미터 튜닝(hyperparameter tuning)이나 모델 선택(model selection)에 사용되며, 일반화된 성능 추정치를 제공합니다.
    *   **부트스트랩:** 최종 모델의 계수나 성능 지표에 대한 신뢰 구간 추정, 또는 Harrell의 optimism correction을 통한 편향 보정에 주로 사용됩니다.
*   **병행 사용 시 표준 보고 형식 (예시):**
    "본 연구에서는 예측 모델의 견고한 성능 평가를 위해 10-fold 교차 검증과 부트스트랩 방법을 병행하여 사용하였다.
    먼저, 모델의 최적 하이퍼파라미터를 결정하고 초기 성능을 평가하기 위해 10-fold 교차 검증을 수행하였다. 데이터셋을 10개의 폴드(fold)로 나누고, 각 폴드를 한 번씩 검증 세트(validation set)로 사용하여 모델을 훈련하고 평가하였다. 이 과정을 통해 얻은 평균 C-index는 [X.XX] (표준편차 [Y.YY])로 나타났다.
    이후, 최종 모델의 계수 추정치와 95% 신뢰 구간의 안정성을 확보하고, 모델 성능의 낙관적 편향을 보정하기 위해 부트스트랩 재

### 🔹 SPSS_Expert

안녕하십니까, "SPSS_Expert"입니다. "Statistician"님의 전문적인 의견에 감사드립니다. 저는 SPSS Statistics의 BOOTSTRAP 명령어, PROCESS macro, 그리고 Cox/로지스틱 회귀에서의 Bootstrap 구현 방식에 대한 깊은 이해를 바탕으로, 의학 논문에서의 SPSS Bootstrap 표준 사용법 및 검증 절차에 대해 구체적이고 실행 가능한 답변을 드리겠습니다. 특히 SPSS 사용자의 관점에서 기능의 활용과 한계, 그리고 보고 형식에 초점을 맞추겠습니다.

---

### **SPSS_Expert의 전문적 답변**

"Statistician"님께서 Bootstrapping의 통계적 원리와 역할에 대해 매우 명확하게 설명해 주셨습니다. SPSS Statistics에서 제공하는 `BOOTSTRAP` 명령어는 이러한 원리를 충실히 따르며, 사용자가 복잡한 코딩 없이도 다양한 통계량의 표준 오차(standard error)와 신뢰 구간(confidence interval)을 추정할 수 있도록 지원합니다.

#### **(1) TRIPOD 지침이나 STROBE에서 bootstrap 내부 검증을 어떻게 요구하는지**

"Statistician"님의 설명에 전적으로 동의합니다. TRIPOD 지침은 예측 모델의 내부 검증(internal validation)을 위해 부트스트랩을 강력히 권장하며, STROBE 지침은 통계 분석의 투명성과 신뢰성 확보라는 맥락에서 부트스트랩 사용을 지지합니다. SPSS 사용자의 관점에서 이를 어떻게 구현하고 보고할 수 있는지 보완 설명 드리겠습니다.

*   **TRIPOD 지침과 SPSS Bootstrap:**
    *   TRIPOD는 모델의 과적합(overfitting)을 평가하고 일반화 성능(generalizability)을 추정하기 위해 부트스트랩을 통한 낙관적 편향(optimism bias) 보정을 권장합니다.
    *   **SPSS에서의 구현 및 한계:** SPSS Statistics의 `BOOTSTRAP` 명령어는 주로 회귀 계수(regression coefficients), 오즈비(Odds Ratio), 위험비(Hazard Ratio) 등 모델 파라미터의 표준 오차와 신뢰 구간을 추정하는 데 사용됩니다. 예를 들어, 로지스틱 회귀(Logistic Regression)나 Cox 회귀(Cox Regression) 모델을 실행할 때 `BOOTSTRAP` 옵션을 선택하면, 각 계수에 대한 부트스트랩 표준 오차와 신뢰 구간(예: Percentile, BCa)을 얻을 수 있습니다.
        *   **메뉴 경로:** `Analyze > Regression > Logistic / Cox Regression...` 에서 `Bootstrap...` 버튼 클릭.
        *   **주요 옵션:**
            *   `Number of samples`: 일반적으로 1000회 이상을 권장합니다.
            *   `Confidence intervals`: `Percentile` 또는 `Bias-corrected and accelerated (BCa)`를 선택할 수 있습니다. BCa는 편향과 왜도(skewness)를 보정하여 더 정확한 신뢰 구간을 제공하는 것으로 알려져 있습니다.
    *   **보고 요구 사항 충족:** SPSS에서 얻은 부트스트랩 신뢰 구간은 TRIPOD에서 요구하는 "모델 성능 지표와 그 신뢰 구간"을 보고하는 데 활용될 수 있습니다. 그러나 TRIPOD가 강조하는 C-index, AUC, Brier score와 같은 *모델 성능 지표 자체의 낙관적 편향 보정*은 기본 SPSS `BOOTSTRAP` 명령어로는 직접적으로 수행하기 어렵습니다. 이는 다음 질문에서 더 자세히 다루겠습니다.

*   **STROBE 지침과 SPSS Bootstrap:**
    *   STROBE는 연구의 방법론적 엄격성과 편향 제어에 중점을 둡니다. SPSS에서 부트스트랩을 사용하여 회귀 계수나 효과 크기(effect size)의 신뢰 구간을 추정하는 것은 분석 결과의 불확실성을 정량화하고 보고하는 데 기여하므로, STROBE의 "통계 분석" 섹션에서 명확히 보고되어야 합니다.
    *   **보고 예시:** "로지스틱 회귀 분석에서 각 독립 변수의 오즈비(Odds Ratio) 및 95% 신뢰 구간은 1000회 부트스트랩 재표본추출(Bias-corrected and accelerated method)을 통해 추정되었다."

#### **(2) Harrell의 optimism correction과 bootstrap의 관계**

"Statistician"님께서 Harrell의 optimism correction의 통계적 원리와 작동 방식을 정확하게 설명해 주셨습니다. 이 방법은 예측 모델의 일반화 성능을 현실적으로 추정하는 데 매우 중요합니다.

*   **SPSS에서의 구현 한계 및 대안:**
    *   **기본 SPSS의 한계:** "Statistician"님께서 언급하셨듯이, **SPSS Statistics의 기본 기능으로는 Harrell의 optimism correction을 직접적으로 수행할 수 없습니다.** SPSS의 `BOOTSTRAP` 명령어는 주로 개별 파라미터의 신뢰 구간 추정에 초점을 맞추며, 모델 전체의 성능 지표(예: C-index, AUC)에 대한 낙관적 편향을 계산하고 보정하는 자동화된 절차를 제공하지 않습니다.
    *   **SPSS 사용자를 위한 대안:**
        1.  **SPSS Python/R 확장 기능 활용:** 가장 현실적이고 권장되는 방법입니다. SPSS Statistics는 Python 또는 R 프로그래밍 언어와의 연동을 지원합니다.
            *   **Python:** `scikit-learn` (AUC), `lifelines` (C-index for survival models) 등의 라이브러리를 사용하여 Harrell의 optimism correction 절차를 직접 코딩할 수 있습니다. 예를 들어, `lifelines`의 `ConcordanceStatistic.bootstrap_score` 함수는 부트스트랩을 이용한 C-index 추정 및 편향 보정에 활용될 수 있습니다.
            *   **R:** `rms` 패키지(Frank Harrell 교수 개발)는 `validate` 함수를 통해 C-index, R-squared, Brier score 등 다양한 성능 지표에 대한 부트스트랩 기반 낙관적 편향 보정을 제공합니다. SPSS 내에서 R 코드를 실행하여 이 기능을 사용할 수 있습니다.
            *   **구현 절차 (개요):**
                *   SPSS에서 데이터를 R 또는 Python으로 넘깁니다.
                *   R/Python에서 `rms` 또는 `lifelines` 등의 패키지를 사용하여 Harrell의 optimism correction 절차를 수행합니다.
                *   보정된 성능 지표를 다시 SPSS로 가져와 보고합니다.
        2.  **수동 구현 (매우 비효율적):** 이론적으로는 SPSS 구문(syntax)을 사용하여 "Statistician"님이 설명하신 4단계 절차를 수동으로 반복 구현할 수 있습니다. 그러나 이는 매우 복잡하고 오류 발생 가능성이 높으며, 대규모 반복(예: 200회 이상)에는 적합하지 않습니다. 각 부트스트랩 표본마다 모델을 개발하고, 원본 데이터와 부트스트랩 데이터에 적용하여 성능 지표를 계산한 후, 그 차이를 저장하고 평균을 내는 일련의 과정을 모두 수동으로 코딩해야 합니다. 이는 SPSS의 강점인 사용자 친화적인 GUI 환경을 벗어나므로 권장되지 않습니다.

**결론:** SPSS 사용자가 Harrell의 optimism correction을 수행하고자 한다면, SPSS의 Python 또는 R 확장 기능을 활용하여 `rms` (R) 또는 `lifelines` (Python)와 같은 전문 패키

### 🔹 Reproducer

안녕하십니까, "재현성 검증자" Reproducer입니다. "Statistician"님과 "SPSS_Expert"님의 전문적인 의견에 깊이 공감하며, 의학 논문에서의 SPSS Bootstrap 사용법과 관련하여 Python을 통한 재현성 검증 절차를 구체적으로 설계하겠습니다. 특히 난수 시드 고정, 리샘플링 단위, 동일 데이터셋, 동일 통계 모델 지정 등 재현 조건을 명확히 하고, Python 코드로 검증 가능한 절차를 제안하는 데 집중하겠습니다.

---

### **재현성 검증을 위한 기본 원칙 (General Principles for Reproducibility)**

SPSS에서 수행된 부트스트랩 분석 결과를 Python에서 동일하게 재현하기 위해서는 다음의 핵심 원칙들을 철저히 준수해야 합니다.

1.  **동일 데이터셋 (Identical Dataset):**
    *   SPSS에서 사용된 원본 데이터셋을 Python으로 정확히 가져와야 합니다. 데이터 형식(예: 숫자형, 범주형), 결측치 처리 방식, 변수 인코딩(예: 더미 변수 생성) 등이 완벽하게 일치해야 합니다.
    *   **권장 방법:** SPSS에서 데이터를 CSV, Parquet 또는 Feather 형식으로 내보내고, Python `pandas` 라이브러리를 사용하여 불러옵니다.

2.  **동일 통계 모델 (Identical Statistical Model):**
    *   종속 변수(dependent variable), 독립 변수(independent variables), 공변량(covariates), 상호작용 항(interaction terms), 분포 가정(distribution assumption), 연결 함수(link function) 등 모델의 모든 사양이 SPSS와 Python에서 동일하게 정의되어야 합니다.
    *   예: 로지스틱 회귀(Logistic Regression)의 경우, 이항 분포와 로짓 연결 함수를 사용하고, Cox 비례 위험 모델(Cox Proportional Hazards Model)의 경우, 비례 위험 가정을 따릅니다.

3.  **동일 리샘플링 단위 (Identical Resampling Unit):**
    *   부트스트랩은 원본 데이터에서 복원 추출(sampling with replacement)을 통해 이루어집니다. 이 "단위"가 개별 환자(individual patient)인지, 아니면 특정 그룹(예: 병원, 클러스터)인지 명확히 해야 합니다. 의학 논문에서는 일반적으로 개별 환자 단위의 복원 추출이 사용됩니다.

4.  **동일 부트스트랩 반복 횟수 (Identical Number of Bootstrap Samples):**
    *   SPSS에서 지정한 부트스트랩 반복 횟수(예: 1000회, 2000회)를 Python에서도 동일하게 적용해야 합니다.

5.  **난수 시드 고정 (Fixed Random Seed):**
    *   **가장 중요합니다.** 부트스트랩은 무작위 재표본추출 과정이므로, 난수 생성 시드(random seed)를 고정하지 않으면 매번 다른 결과가 나옵니다. SPSS와 Python 모두에서 동일한 난수 시드를 설정해야 합니다.
        *   **SPSS:** `SET SEED = [숫자].` 명령어를 사용하거나, `BOOTSTRAP` 대화 상자에서 `Set seed for Mersenne Twister` 옵션을 체크하고 특정 숫자를 입력합니다.
        *   **Python:** `numpy.random.seed([숫자])` 및 `random.seed([숫자])`를 사용하여 난수 시드를 고정합니다. `scikit-learn` 등 일부 라이브러리는 자체적으로 `random_state` 매개변수를 제공합니다.

6.  **동일 신뢰 구간 추정 방식 (Identical Confidence Interval Method):**
    *   부트스트랩 신뢰 구간(Confidence Interval, CI)을 계산하는 방식(예: Percentile, Bias-corrected and accelerated (BCa), Normal)이 SPSS와 Python에서 동일해야 합니다. BCa 방법은 편향(bias)과 왜도(skewness)를 보정하여 더 정확한 신뢰 구간을 제공하는 것으로 알려져 있습니다.

---

### **Python을 통한 재현성 검증 절차 설계**

SPSS_Expert님께서 지적하셨듯이, SPSS의 기본 `BOOTSTRAP` 명령어는 주로 개별 파라미터(예: 회귀 계수, 오즈비, 위험비)의 신뢰 구간 추정에 초점을 맞춥니다. Harrell의 optimism correction과 같은 모델 성능 지표의 편향 보정은 SPSS에서 직접 지원하지 않으므로, 이 부분은 Python의 강점을 활용하여 구현하고 검증해야 합니다.

#### **1단계: 데이터 준비 및 환경 설정**

1.  **SPSS 데이터 내보내기:**
    *   SPSS에서 분석에 사용된 데이터를 CSV 파일로 내보냅니다. (예: `File > Save As... > Save as type: Comma Separated Value (*.csv)`)
    *   범주형 변수는 Python에서 적절히 처리할 수 있도록 숫자형 코드로 변환하거나, SPSS에서 직접 더미 변수(dummy variables)를 생성하여 내보내는 것을 고려합니다.
2.  **Python 환경 설정:**
    *   필요 라이브러리 설치: `pandas`, `numpy`, `scipy`, `statsmodels`, `scikit-learn`, `lifelines` (생존 분석용).
    *   **난수 시드 고정:**
        ```python
        import numpy as np
        import pandas as pd
        import random

        # SPSS에서 사용한 것과 동일한 시드 값으로 설정
        SPSS_SEED = 12345 # 예시 시드 값
        np.random.seed(SPSS_SEED)
        random.seed(SPSS_SEED)
        # scikit-learn 등에서 random_state를 사용하는 경우에도 동일한 시드 사용
        ```
3.  **데이터 불러오기 및 전처리:**
    ```python
    df = pd.read_csv('your_spss_data.csv')
    # SPSS에서 수행된 것과 동일한 결측치 처리, 범주형 변수 인코딩 (예: pd.get_dummies)
    # 예: df = pd.get_dummies(df, columns=['categorical_var'], drop_first=True)
    ```

#### **2단계: SPSS `BOOTSTRAP` 명령어 결과 재현 (개별 파라미터 신뢰 구간)**

이 단계는 SPSS에서 로지스틱 회귀 또는 Cox 회귀 분석 시 `BOOTSTRAP` 옵션을 사용하여 얻은 회귀 계수(regression coefficients), 오즈비(Odds Ratio, OR) 또는 위험비(Hazard Ratio, HR)의 표준 오차(Standard Error, SE) 및 신뢰 구간(CI)을 Python에서 재현하는 것을 목표로 합니다.

1.  **SPSS 분석 수행:**
    *   예: `Analyze > Regression > Binary Logistic...` 또는 `Analyze > Survival > Cox Regression...`
    *   `Bootstrap...` 버튼 클릭:
        *   `Number of samples`: 1000 (예시)
        *   `Confidence intervals`: `Bias-corrected and accelerated (BCa)` (권장) 또는 `Percentile`
        *   `Set seed for Mersenne Twister`: 체크하고 `SPSS_SEED` 값 입력.
    *   결과 테이블(Coefficients, Variables in the Equation 등)에서 부트스트랩 SE 및 CI 값을 기록합니다.

2.  **Python 코드 구현 (예: 로지스틱 회귀 계수 BCa CI 재현):**
    ```python
    import statsmodels.api as sm
    from scipy.stats import norm

    # 모델 정의 (SPSS와 동일하게)
    X = df[['independent_var1', 'independent_var2', 'covariate']]
    X = sm.add_constant(X) # 절편(constant) 추가
    y = df['dependent_binary_var']

    n_bootstrap_samples = 1000 # SPSS와 동일한 반복 횟수
    bootstrap_coefficients = []

    # 부트스트랩 루프
    for i in range(n_bootstrap_samples):
        # 복원 추출 (SPSS와 동일한 리샘플링 단위)
        bootstrap_sample = df.sample(n=len(df), replace=True, random_state=SPSS_SEED + i) # 각 샘플마다 시드 변경
        
        X_b = bootstrap_sample[['independent_var1', 'independent_var2', 'covariate']]
        X_b = sm.add_constant(X_b)
        y_b = bootstrap_sample['dependent_binary_var']

---

## 부록 A: 공식 출처 검증 (웹 검색 기반, 2026-08-17)

> 본 절은 Tavily 검색 실패로 인해 누락됐던 **공식 문서 출처**를 별도 웹 검색으로
> 확보하여 추가한 내용입니다. (SPSS 공식 문서 + TRIPOD 지침)

### A.1 SPSS BOOTSTRAP 공식 명령어 (재현성 확보)

IBM 공식 문서 기준, SPSS에서 bootstrap 재현(reproducibility)을 보장하려면
`SET` 명령어로 난수 생성기를 먼저 설정한 뒤 `BOOTSTRAP`을 실행해야 합니다.
**`BOOTSTRAP` 명령 자체에는 seed 파라미터가 없으며**, 세션 전역 난수 설정을 따릅니다.

```spss
* 1) Mersenne Twister + seed 설정 (IBM 공식 권장 방식)
SET RNG=MT MTINDEX=9191972.

* 2) Bootstrap 수행 (1000회, 단순 무작위)
BOOTSTRAP /CRITERIA NSAMPLES=1000
  /VARIABLES INPUT=종속변수
  /SAMPLING METHOD=SIMPLE.

* 3) 분석 프로시저는 BOOTSTRAP 직후에 실행
COXREG 시간 WITH 공변량1 공변량2
  /STATUS=사건변수(1)
  /METHOD=ENTER 공변량1 공변량2.
```

| 구성 요소 | 의미 |
|---|---|
| `SET RNG=MT MTINDEX=값` | Mersenne Twister(`MT`) 엔진 + 시작점(`MTINDEX`) 고정 → 동일 seed로 반복 실행 시 **동일 bootstrap 표본** 생성 |
| `BOOTSTRAP /CRITERIA NSAMPLES=1000` | 부트스트랩 반복 횟수 1000 (기본값이 1000이지만 명시 권장) |
| `/SAMPLING METHOD=SIMPLE` | 리샘플링 방법: `SIMPLE`(단순 무작위) / `STRATIFIED`(계층화) 옵션 |
| 프로시저 위치 | `BOOTSTRAP`은 **바로 뒤따르는 프로시저에만** 적용됨 |

> ⚠️ **주의**: `BOOTSTRAP` 명령은 "바로 다음에 오는 분석 프로시저에만" 적용됩니다.
> 여러 분석을 bootstrap하려면 각각 앞에 `BOOTSTRAP` 명령을 붙여야 합니다.

### A.2 TRIPOD / TRIPOD+AI 지침 (2024) — 내부 검증 요구사항

- **내부 검증(internal validation)은 선택이 아닌 필수**: 모델 성능의
  낙관적 편향(optimism)을 정량화하기 위해 **bootstrapping 또는 cross-validation**
  같은 resampling 기법을 TRIPOD가 강력히 권장
- **단순 train-test 분할(70/30 split)은 비권장**: 비효율적이고 성능 추정치의
  변동성이 큼 → 70/30 internal validation만으로는 부족, bootstrap/CV 병행 필요
- **중요(Item 12c)**: 내부 검증 시 **모델 개발 전 과정(전처리, 변수 선택,
  하이퍼파라미터 튜닝)을 각 resampling 반복 안에서 재실행**해야 함.
  루프 밖에서 하면 data leakage로 낙관적 편향 발생
- bootstrap vs CV: 임상 예측 모델 문헌에서는 **bootstrap이 낮은 bias로 선호**되는
  경향, 반복 k-fold CV도 강력한 대안

> **우리 연구(OSCC, n=133) 적용 시사점**:
> - 70/30 split 단독 사용은 TRIPOD+AI 기준에서 비권장 → SPSS Bootstrap + 10-fold CV 병행이 올바른 방향
> - SPSS BOOTSTRAP은 "개별 파라미터 CI"만 제공하므로, 모델 성능(C-index)의 optimism 보정은
>   Python(`lifelines`) 또는 R(`rms`)로 수행 필요

### A.3 참고 출처

| 출처 | 내용 |
|---|---|
| IBM SPSS Statistics Documentation | BOOTSTRAP 명령어 구문, SET RNG/MTINDEX, NSAMPLES, SAMPLING METHOD |
| IBM SPSS Bootstrapping 공식 페이지 | https://www.ibm.com/products/spss-statistics/bootstrapping |
| TRIPOD+AI Statement (April 2024) | 내부 검증 지침, Item 12c, bootstrap vs CV |
| Efron & Tibshirani (1993) | 비모수 부트스트랩 표준 방법론, BCa 신뢰 구간 |
| Harrell (rms 패키지) | validate() 함수를 통한 optimism correction |

---

