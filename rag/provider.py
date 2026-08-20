import os
import time
from pathlib import Path

from dotenv import load_dotenv
from openai import APIStatusError, BadRequestError, OpenAI, RateLimitError

REPO_ROOT = Path(__file__).resolve().parent.parent
load_dotenv(REPO_ROOT / ".env")

DEFAULT_BASE_URL = "https://integrate.api.nvidia.com/v1"
DEFAULT_MODEL = "nvidia/nemotron-3-ultra-550b-a55b"
GLM_MODEL = "z-ai/glm-5.2"
LLAMA_MODEL = "meta/llama-3.3-70b-instruct"


# The project's NVIDIA API is rate-limited to 40 calls/minute (RAG_ARCHITECTURE.md SS5) -- this
# is the floor that protects against ever going over it, no matter which caller forgets to
# count its own calls. It does NOT replace the "one generation call per result set, not per
# cell line" design rule -- that discipline belongs to each calling module.
#
# Disclosed limit of this guard: it is a per-*process* sliding window, so it cannot see calls
# made by a different process sharing the same API key -- observed live during this project's
# own testing (several separate `pytest` invocations in quick succession each stayed under 40
# within themselves, but collectively tripped the server's real 429). The RETRY_ON_RATE_LIMIT
# backoff below is what actually recovers from that case; a production deployment with multiple
# worker processes would need a shared counter (Redis or similar), not this in-memory list.
CALLS_PER_MINUTE = 40
RATE_LIMIT_RETRY_DELAYS = (5, 15, 30)  # seconds -- exhausts in under a minute, matching the window

_client = None
_call_timestamps = []  # sliding one-minute window of chat() call times, oldest first


class ProviderRequestError(RuntimeError):
    """Raised when the LLM API rejects a request outright (e.g. context_length_exceeded) --
    distinct from RateLimitError (retried below) since retrying an oversized prompt unchanged
    would just fail identically."""


def _get_client():
    global _client
    if _client is None:
        api_key = os.environ.get("NVIDIA_API_KEY")
        if not api_key:
            raise RuntimeError(
                "NVIDIA_API_KEY is not set. Copy .env.example to .env at the repo root and "
                "fill in a real key -- see rag/_.md's Secrets section."
            )
        base_url = os.environ.get("NVIDIA_BASE_URL", DEFAULT_BASE_URL)
        _client = OpenAI(api_key=api_key, base_url=base_url)
    return _client


def _respect_rate_limit():
    now = time.monotonic()
    while _call_timestamps and now - _call_timestamps[0] > 60:
        _call_timestamps.pop(0)
    if len(_call_timestamps) >= CALLS_PER_MINUTE:
        sleep_for = 60 - (now - _call_timestamps[0])
        if sleep_for > 0:
            time.sleep(sleep_for)
    _call_timestamps.append(time.monotonic())


def chat(messages, tools=None, model=None, temperature=0.0, response_format=None):
    client = _get_client()
    kwargs = {
        "model": model or os.environ.get("NVIDIA_MODEL", DEFAULT_MODEL),
        "messages": messages,
        "temperature": temperature,
    }
    if tools:
        kwargs["tools"] = tools
    if response_format:
        kwargs["response_format"] = response_format

    for delay in (*RATE_LIMIT_RETRY_DELAYS, None):
        _respect_rate_limit()
        try:
            completion = client.chat.completions.create(**kwargs)
            return completion.choices[0].message
        except RateLimitError:
            if delay is None:
                raise
            time.sleep(delay)
        except (BadRequestError, APIStatusError) as exc:
            # Not retryable -- most commonly context_length_exceeded from an oversized prompt.
            # Retrying the same request would fail identically, so fail closed immediately with
            # a message callers (narrator.py) can surface instead of an opaque SDK traceback.
            raise ProviderRequestError(f"LLM request rejected by the API: {exc}") from exc
