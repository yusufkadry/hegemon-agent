"""Minimal OpenAI Responses API client (standard library only).

- One request per call, with a hard timeout. The only retries are one visible
  retry for rate limits / server errors, printed to the terminal.
- Key, quota, model-access and timeout problems raise ProviderError, which the
  agent reports as an AI-provider problem, never as "dataset failed".
"""
from getpass import getpass
import json
import os
import socket
import time
import urllib.error
import urllib.request

from . import ui

API_URL = os.environ.get("HEGEMON_OPENAI_URL", "https://api.openai.com/v1/responses")


class ProviderError(Exception):
    def __init__(self, kind, message):
        Exception.__init__(self, message)
        self.kind = kind


def load_key(cfg, prompt=True):
    key = (os.environ.get("OPENAI_API_KEY") or "").strip()
    if key:
        return key
    if cfg.key_file.is_file():
        key = cfg.key_file.read_text().strip()
        if key:
            return key
    if not prompt:
        return None
    print("\n   OpenAI API key needed (input hidden). It is saved to {} (mode 600).".format(cfg.key_file))
    key = getpass("   OpenAI API key: ").strip()
    if key:
        save_key(cfg, key)
    return key or None


def save_key(cfg, key):
    cfg.key_file.parent.mkdir(parents=True, exist_ok=True)
    os.chmod(str(cfg.key_file.parent), 0o700)
    tmp = cfg.key_file.with_suffix(".tmp")
    tmp.write_text(key.strip() + "\n")
    os.chmod(str(tmp), 0o600)
    os.replace(str(tmp), str(cfg.key_file))


def _error_from_http(code, body):
    try:
        err = json.loads(body).get("error") or {}
    except ValueError:
        err = {}
    msg = err.get("message") or body[:300] or "no details"
    etype = (err.get("code") or err.get("type") or "")
    if code == 401:
        return ProviderError("auth", "OpenAI rejected the API key (401). Run: hegemon key")
    if code == 429 and etype == "insufficient_quota":
        return ProviderError("quota", "OpenAI account is out of credit/quota (429 insufficient_quota): " + msg)
    if code == 429:
        return ProviderError("rate", "OpenAI rate limit (429): " + msg)
    if code in (403, 404) or "model" in etype:
        return ProviderError("model", "OpenAI refused the model ({}): {}. Set HEGEMON_MODEL to a model your key can use.".format(code, msg))
    if code >= 500:
        return ProviderError("server", "OpenAI server error {}: {}".format(code, msg))
    return ProviderError("bad_request", "OpenAI request rejected ({}): {}".format(code, msg))


class LLM(object):
    def __init__(self, key, cfg):
        self.key = key
        self.cfg = cfg
        self.calls = 0
        self.input_tokens = 0
        self.output_tokens = 0

    def _post(self, payload, timeout):
        data = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(API_URL, data=data, headers={
            "Authorization": "Bearer " + self.key,
            "Content-Type": "application/json",
        })
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as err:
            raise _error_from_http(err.code, err.read().decode("utf-8", "replace"))
        except socket.timeout:
            raise ProviderError("timeout", "No answer from OpenAI within {:.0f}s (HEGEMON_AI_TIMEOUT).".format(timeout))
        except urllib.error.URLError as err:
            if isinstance(getattr(err, "reason", None), socket.timeout):
                raise ProviderError("timeout", "No answer from OpenAI within {:.0f}s.".format(timeout))
            raise ProviderError("network", "Could not reach OpenAI: {}".format(err.reason))

    def ask(self, instructions, text, model=None, effort=None, max_output_tokens=32000,
            timeout=None, label="AI"):
        payload = {
            "model": model or self.cfg.model,
            "instructions": instructions,
            "input": text,
            "max_output_tokens": max_output_tokens,
        }
        effort = self.cfg.reasoning if effort is None else effort
        if effort:
            payload["reasoning"] = {"effort": effort}
        timeout = timeout or self.cfg.ai_timeout
        retried = False
        while True:
            try:
                with ui.Heartbeat("waiting on {} ({})".format(label, payload["model"])) as hb:
                    data = self._post(payload, timeout)
                break
            except ProviderError as err:
                if err.kind == "bad_request" and "reasoning" in str(err) and "reasoning" in payload:
                    ui.warn("model does not take a reasoning setting; retrying without it")
                    payload.pop("reasoning")
                    continue
                if err.kind in ("rate", "server") and not retried:
                    retried = True
                    ui.warn(str(err) + " -- retrying once in 20s")
                    time.sleep(20)
                    continue
                raise
        self.calls += 1
        usage = data.get("usage") or {}
        self.input_tokens += int(usage.get("input_tokens") or 0)
        self.output_tokens += int(usage.get("output_tokens") or 0)
        text_out = []
        for item in data.get("output") or []:
            if item.get("type") == "message":
                for part in item.get("content") or []:
                    if part.get("type") in ("output_text", "text"):
                        text_out.append(part.get("text") or "")
        result = "".join(text_out).strip()
        ui.info("{} answered in {:.0f}s ({} in / {} out tokens)".format(
            label, hb.elapsed, usage.get("input_tokens", "?"), usage.get("output_tokens", "?")))
        if data.get("status") == "incomplete" and not result:
            reason = (data.get("incomplete_details") or {}).get("reason", "unknown")
            raise ProviderError("incomplete", "OpenAI stopped before answering (reason: {}).".format(reason))
        if not result:
            raise ProviderError("empty", "OpenAI returned no text.")
        return result

    def ping(self, model):
        """A real 1-request inference check (proves key, quota and model access)."""
        self._post({"model": model, "input": "Reply with OK.", "max_output_tokens": 64}, timeout=120)
