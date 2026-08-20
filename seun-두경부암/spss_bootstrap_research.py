#!/usr/bin/env python3
"""
SPSS Bootstrapping 방법론 조사 — 멀티 에이전트 시스템
=====================================================
목적: "SPSS Bootstrapping을 사용해도 정세운 선생이 Python으로 직접 분석한 결과와
동일하게 나오는가?" 라는 질문에 답하기 위한 방법론 조사.

multi_agent.py (사업계획서 협업) 구조를 SPSS 조사에 맞게 각색.

Agents:
  - Statistician (통계 방법론자)  → Google Gemini 2.5 Flash (API)
  - SPSS_Expert (SPSS 전문가)    → Ollama Mistral (로컬)
  - Reproducer (재현성 검증자)   → Ollama Llama 3.1 (로컬)
  - Researcher (문헌 조사원)     → Tavily Search (API)
"""

import os
import json
import time
import sys
from datetime import datetime
from typing import Optional

try:
    from dotenv import load_dotenv
except ImportError:
    print("❌ python-dotenv not installed. Run: pip3 install python-dotenv")
    sys.exit(1)

load_dotenv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '.env'))

AGENTS = [
    {
        "name": "Statistician",
        "role": "통계 방법론자",
        "description": "Bootstrapping의 통계적 원리, resampling 전략, "
                       "confidence interval 추정 방식(percentile, BCa), "
                       "bias correction 등 방법론적 정확성을 검증합니다. "
                       "Python(scikit-learn, lifelines, statsmodels) 구현과의 "
                       "수학적 동일성 판단을 담당합니다.",
        "provider": "gemini",
        "model": "gemini-2.5-flash",
        "api_key_env": "GEMINI_API_KEY",
        "color": "\033[94m",
    },
    {
        "name": "SPSS_Expert",
        "role": "SPSS 전문가",
        "description": "IBM SPSS Statistics의 BOOTSTRAP 명령어, "
                       "PROCESS macro, Cox 회귀/로지스틱 회귀에서의 "
                       "Bootstrap 구현 방식, 기본 옵션(simple random vs "
                       "stratified), seed 설정, 1000회 resampling의 실제 "
                       "동작을 전문적으로 분석합니다.",
        "provider": "gemini",
        "model": "gemini-2.5-flash",
        "api_key_env": "GEMINI_API_KEY",
        "color": "\033[92m",
    },
    {
        "name": "Reproducer",
        "role": "재현성 검증자",
        "description": "SPSS 결과가 Python 분석 결과와 동일하게 나오는지 "
                       "재현하기 위한 구체적 방법을 설계합니다. "
                       "난수 시드 고정, resampling 단위, 동일 데이터셋, "
                       "동일 통계 모델 지정 등 재현 조건을 점검하고 "
                       "Python 코드로 검증 가능한 절차를 제안합니다.",
        "provider": "gemini",
        "model": "gemini-2.5-flash",
        "api_key_env": "GEMINI_API_KEY",
        "color": "\033[93m",
    },
    {
        "name": "Researcher",
        "role": "문헌 조사원",
        "description": "웹 검색 기반으로 SPSS Bootstrapping의 공식 문서, "
                       "IBM 기술 자료, 의학 논문에서의 사용 사례, "
                       "Bootstrap vs Cross-validation 비교 연구, "
                       "SPSS와 R/Python 결과 일치성 관련 문헌을 조사합니다.",
        "provider": "tavily",
        "model": "search",
        "api_key_env": "TAVILY_API_KEY",
        "color": "\033[95m",
    },
]

RESET = "\033[0m"


def call_gemini(prompt: str, model: str = "gemini-2.5-flash") -> Optional[str]:
    import requests

    api_key = os.getenv("GEMINI_API_KEY")
    if not api_key:
        print("❌ GEMINI_API_KEY not set in .env file")
        return None

    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent?key={api_key}"
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {
            "temperature": 0.3,
            "maxOutputTokens": 4096,
        },
    }

    try:
        resp = requests.post(url, json=payload, timeout=90)
        resp.raise_for_status()
        data = resp.json()
        candidates = data.get("candidates", [])
        if candidates:
            parts = candidates[0].get("content", {}).get("parts", [])
            if parts:
                return parts[0].get("text", "")
        return None
    except Exception as e:
        print(f"❌ Gemini API error: {e}")
        return None


def call_ollama(prompt: str, model: str = "llama3.1") -> Optional[str]:
    import requests

    payload = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "stream": False,
        "options": {"temperature": 0.3},
    }

    try:
        resp = requests.post(
            "http://localhost:11434/api/chat",
            json=payload,
            timeout=300,
        )
        resp.raise_for_status()
        data = resp.json()
        return data["message"]["content"]
    except requests.exceptions.ConnectionError:
        print("❌ Ollama not running. Start it: ollama serve")
        return None
    except Exception as e:
        print(f"❌ Ollama error: {e}")
        return None


def call_tavily(prompt: str) -> Optional[str]:
    api_key = os.getenv("TAVILY_API_KEY")
    if not api_key:
        print("❌ TAVILY_API_KEY not set in .env file")
        return None

    import requests

    payload = {
        "api_key": api_key,
        "query": prompt,
        "search_depth": "advanced",
        "include_answer": True,
        "include_raw_content": False,
        "max_results": 6,
    }

    try:
        resp = requests.post(
            "https://api.tavily.com/search",
            json=payload,
            timeout=90,
        )
        resp.raise_for_status()
        data = resp.json()

        answer = data.get("answer", "")
        results = data.get("results", [])

        output = f"[Tavily Search Results]\n\n"
        if answer:
            output += f"📌 요약:\n{answer}\n\n"
        if results:
            output += "📄 검색 결과:\n"
            for i, r in enumerate(results[:6], 1):
                title = r.get("title", "Untitled")
                url = r.get("url", "")
                content = r.get("content", "")
                output += f"\n{i}. {title}\n   URL: {url}\n   {content[:600]}\n"
        return output
    except Exception as e:
        print(f"❌ Tavily API error: {e}")
        return None


def call_agent(agent: dict, prompt: str) -> Optional[str]:
    provider = agent["provider"]
    model = agent["model"]

    print(f"\n{agent['color']}  ┌─ [{agent['name']}] ({agent['role']}) ─────────────{RESET}")
    print(f"{agent['color']}  │  모델: {model}{RESET}")
    print(f"{agent['color']}  │  호출 중...{RESET}")

    if provider == "gemini":
        return call_gemini(prompt, model)
    elif provider == "ollama":
        return call_ollama(prompt, model)
    elif provider == "tavily":
        return call_tavily(prompt)
    else:
        print(f"❌ Unknown provider: {provider}")
        return None


RESEARCH_QUESTIONS = [
    {
        "id": "Q1",
        "title": "SPSS Bootstrapping의 알고리즘: simple random vs stratified, resampling 전략",
        "prompt": (
            "SPSS Statistics의 BOOTSTRAP 프로시저가 실제로 어떤 알고리즘을 사용하는지 "
            "상세히 분석하세요. (1) 기본 resampling 방법이 simple random sampling인지 "
            "stratified sampling인지, (2) Cox 회귀분석과 로지스틱 회귀분석에서 bootstrap이 "
            "어떻게 적용되는지(원데이터 재추출 vs 잔차 재추출), (3) 1000회 resampling의 "
            "구체적 의미, (4) bootstrap confidence interval 계산 방식(percentile, BCa, "
            "bias-corrected), (5) random seed 설정 방법을 포함하세요."
        ),
    },
    {
        "id": "Q2",
        "title": "SPSS vs Python 재현성: 동일 결과가 나오는 조건",
        "prompt": (
            "SPSS Bootstrapping 결과가 Python(예: lifelines, statsmodels, scikit-learn)으로 "
            "직접 계산한 결과와 '동일하게' 나오기 위한 조건을 분석하세요. "
            "(1) 어떤 요소가 결과 차이를 만들 수 있는지(난수 생성기 차이, resampling 단위, "
            "CI 방식, 알고리즘 구현 차이), (2) 정확히 동일한 숫자를 기대하는 것이 "
            "통계적으로 올바른 기대인지, (3) '동일한 결과'가 아니라 '동일한 통계적 "
            "결론(유의성/방향/대략적 크기)'을 검증하는 것이 올바른 기준인지, "
            "(4) 완전 재현을 원한다면 난수 시드를 어떻게 맞춰야 하는지 판단하세요."
        ),
    },
    {
        "id": "Q3",
        "title": "의학 논문에서의 SPSS Bootstrap 표준 사용법 및 검증 절차",
        "prompt": (
            "의학/임상 논문(특히 생존분석, 예후 예측 모델)에서 SPSS Bootstrapping을 "
            "표준적으로 어떻게 사용하고 보고하는지 조사하세요. (1) TRIPOD 지침이나 "
            "STROBE에서 bootstrap 내부 검증을 어떻게 요구하는지, (2) Harrell의 "
            "optimism correction과 bootstrap의 관계, (3) SPSS bootstrap을 "
            "10-fold cross-validation과 병행할 때의 표준 보고 형식, "
            "(4) 소표본(n=133)에서 bootstrap 사용의 적절성과 한계를 포함하세요."
        ),
    },
]

SYSTEM_PROMPT_TEMPLATE = """당신은 "{role}" 역할을 맡은 AI 에이전트 "{name}"입니다.
{description}

## 조사 주제
{question_title}

## 조사 질문
{question_prompt}

## 지금까지 다른 에이전트들의 의견
{context}

당신의 역할에 맞춰 위 조사 질문에 대한 전문적 답변을 작성해주세요.
- 이전 에이전트의 의견에 동의/반박/보완할 수 있습니다.
- 구체적이고 실행 가능한 내용을 제시하세요.
- 한국어로 작성하되, 전문 용어는 영어를 병기해주세요.
- 인용이 필요한 경우 출처를 명시하세요.
"""


def run_multi_agent_research():
    print(f"\n{'='*60}")
    print(f"  🔬 SPSS Bootstrapping 방법론 조사 — 멀티 에이전트")
    print(f"{'='*60}")
    print(f"  에이전트: {len(AGENTS)}명")
    print(f"  조사 질문: {len(RESEARCH_QUESTIONS)}개")
    print(f"\n{'='*60}\n")

    final_report = {
        "metadata": {
            "topic": "SPSS Bootstrapping 방법론 조사",
            "goal": "SPSS Bootstrapping 사용 시 Python 직접 분석 결과와 동일하게 나오는지 검증",
            "generated_at": datetime.now().isoformat(),
            "agents": [a["name"] for a in AGENTS],
        },
        "questions": {},
    }

    for question in RESEARCH_QUESTIONS:
        qid = question["id"]
        title = question["title"]
        print(f"\n{'='*60}")
        print(f"  📋 {qid}: {title}")
        print(f"{'='*60}")

        context = f"조사 주제: {title}\n\n"
        q_outputs = []

        for agent in AGENTS:
            prompt = SYSTEM_PROMPT_TEMPLATE.format(
                name=agent["name"],
                role=agent["role"],
                description=agent["description"],
                question_title=title,
                question_prompt=question["prompt"],
                context=context,
            )

            result = call_agent(agent, prompt)

            if result:
                q_outputs.append({"agent": agent["name"], "content": result})
                context += f"\n--- [{agent['name']}]의 의견 ---\n{result}\n"
                print(f"\n{agent['color']}  └─ 응답 완료 ({len(result)}자){RESET}")
            else:
                print(f"\n{agent['color']}  └─ ⚠️  응답 실패{RESET}")

            time.sleep(0.5)

        final_report["questions"][qid] = {
            "title": title,
            "outputs": q_outputs,
        }

    return final_report


def save_report(report: dict):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    json_path = f"spss_bootstrap_research_{timestamp}.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    print(f"\n✅ JSON 저장 완료: {json_path}")

    md_path = f"spss_bootstrap_research_{timestamp}.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# SPSS Bootstrapping 방법론 조사 결과\n\n")
        f.write(f"> 조사 목적: SPSS Bootstrapping 사용 시 Python 직접 분석 결과와 "
                f"동일하게 나오는지 검증\n")
        f.write(f"> 생성일: {report['metadata']['generated_at']}\n")
        f.write(f"> 참여 에이전트: {', '.join(report['metadata']['agents'])}\n\n")
        f.write("---\n\n")

        for qid, qdata in report["questions"].items():
            f.write(f"## {qid}: {qdata['title']}\n\n")
            for output in qdata["outputs"]:
                f.write(f"### 🔹 {output['agent']}\n\n")
                f.write(f"{output['content']}\n\n")
            f.write("---\n\n")

    print(f"✅ Markdown 저장 완료: {md_path}")


if __name__ == "__main__":
    print("=" * 60)
    print("  🔬 SPSS Bootstrapping 방법론 조사")
    print("=" * 60)
    print()
    print("  에이전트 구성:")
    for a in AGENTS:
        provider_icon = {"gemini": "🔷", "ollama": "🟢", "tavily": "🟠"}
        print(f"    {provider_icon.get(a['provider'], '⚪')} [{a['name']}] {a['role']} → {a['provider']}/{a['model']}")
    print()

    report = run_multi_agent_research()

    print(f"\n{'='*60}")
    print("  ✅ 모든 조사 완료!")
    print(f"{'='*60}")

    save_report(report)

    print(f"\n{'='*60}")
    print("  🎉 SPSS Bootstrapping 조사 완료!")
    print(f"{'='*60}")
