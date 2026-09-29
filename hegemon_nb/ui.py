import sys
import threading
import time

_TTY = sys.stdout.isatty()


def _c(code, text):
    return "\033[{}m{}\033[0m".format(code, text) if _TTY else text


def step(msg):
    print(_c("1", "\n== " + msg), flush=True)


def info(msg):
    print("   " + msg, flush=True)


def ok(msg):
    print(_c("32", "   OK ") + msg, flush=True)


def warn(msg):
    print(_c("33", "   WARNING ") + msg, flush=True)


def fail(msg):
    print(_c("31", "   FAILED ") + msg, flush=True)


def ask(prompt, default=""):
    try:
        answer = input("   " + prompt).strip()
    except EOFError:
        return default
    return answer or default


def yes_no(prompt, default=False, assume=None):
    if assume is not None:
        info("{} {}".format(prompt, "yes (--yes)" if assume else "no"))
        return assume
    if not sys.stdin.isatty():
        return default
    suffix = " [Y/n] " if default else " [y/N] "
    answer = ask(prompt + suffix).lower()
    if not answer:
        return default
    return answer.startswith("y")


class Heartbeat(object):
    """Prints 'still working' lines while something slow runs, so a wait never looks like a hang."""

    def __init__(self, label, every=20):
        self.label = label
        self.every = every
        self._stop = threading.Event()
        self._thread = None
        self.start_time = None

    def __enter__(self):
        self.start_time = time.time()
        self._thread = threading.Thread(target=self._run)
        self._thread.daemon = True
        self._thread.start()
        return self

    def _run(self):
        while not self._stop.wait(self.every):
            info("... {} ({:.0f}s)".format(self.label, time.time() - self.start_time))

    def __exit__(self, *exc):
        self._stop.set()
        if self._thread:
            self._thread.join(1)
        return False

    @property
    def elapsed(self):
        return time.time() - self.start_time


def duration(seconds):
    seconds = max(0, seconds)
    return "{:.0f}s".format(seconds) if seconds < 60 else "{}m{:02.0f}s".format(int(seconds // 60), seconds % 60)


class Clock(object):
    """Times each step of a build, so slow steps are visible instead of guessed."""

    def __init__(self):
        self.start = self.last = time.time()
        self.laps = []

    def lap(self, name):
        now = time.time()
        self.laps.append((name, now - self.last))
        self.last = now

    def total(self):
        return time.time() - self.start

    def text(self):
        return " | ".join(["{} {}".format(n, duration(s)) for n, s in self.laps] + ["total " + duration(self.total())])
