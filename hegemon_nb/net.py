"""HTTP with the standard library only: explicit timeouts, bounded retries, size checks."""
import json
import os
import socket
import threading
import time
import urllib.error
import urllib.request

UA = "hegemon-agent/1.0 (BooleanLab dataset builder)"


class NetError(Exception):
    pass


def request(url, data=None, headers=None, timeout=60, retries=3, method=None, max_bytes=None):
    hdrs = {"User-Agent": UA}
    hdrs.update(headers or {})
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, data=data, headers=hdrs, method=method)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                body = resp.read() if max_bytes is None else resp.read(max_bytes)
                return resp.status, dict(resp.headers), body
        except urllib.error.HTTPError as err:
            last = "HTTP {} for {}".format(err.code, url)
            if err.code not in (429, 500, 502, 503, 504) or attempt == retries:
                raise NetError(last)
        except (urllib.error.URLError, socket.timeout, ConnectionError, OSError) as err:
            last = "{} for {}".format(getattr(err, "reason", err), url)
            if attempt == retries:
                raise NetError(last)
        time.sleep(min(30, 2 * (2 ** attempt)))
    raise NetError(last or url)


def get_json(url, **kw):
    return json.loads(request(url, **kw)[2].decode("utf-8"))


def get_text(url, **kw):
    return request(url, **kw)[2].decode("utf-8", errors="replace")


PARTS = int(os.environ.get("HEGEMON_DOWNLOAD_PARTS", "4"))
SPLIT_ABOVE = 16 * 1024 * 1024
_HEADS = {}   # url -> (size, accepts_ranges): a server is never asked the same thing twice in one run
_HEADS_LOCK = threading.Lock()


def head(url, timeout=30):
    """(size, accepts_ranges) from one HEAD request; (None, False) when the server will not say."""
    with _HEADS_LOCK:
        if url in _HEADS:
            return _HEADS[url]
    try:
        _, headers, _ = request(url, method="HEAD", timeout=timeout, retries=1)
    except NetError:
        return None, False
    lower = {k.lower(): v for k, v in headers.items()}
    if "text/html" in str(lower.get("content-type", "")).lower() and \
            not url.split("?", 1)[0].lower().endswith((".html", ".htm")):
        return None, False   # that is an error page's size, not the file's
    size = lower.get("content-length")
    size = int(size) if size and str(size).isdigit() else None
    info = (size, "bytes" in str(lower.get("accept-ranges", "")).lower())
    if size is not None:
        with _HEADS_LOCK:
            _HEADS[url] = info
    return info


def remote_size(url, timeout=30):
    return head(url, timeout)[0]


# What a file's first bytes must be, judged by its name. A server that is throttling or failing
# often answers "200 OK" with an HTML error page; saved under a .gz name, that page is not data.
_MAGIC = (
    ((".gz", ".tgz", ".bgz"), (b"\x1f\x8b",)),
    ((".zip", ".xlsx", ".xlsm", ".docx"), (b"PK\x03\x04", b"PK\x05\x06")),
    ((".bz2",), (b"BZh",)),
    ((".xz",), (b"\xfd7zXZ\x00",)),
    ((".h5", ".h5ad", ".hdf5", ".loom"), (b"\x89HDF\r\n\x1a\n",)),
)
_PAGE_STARTS = (b"<!doctype", b"<html", b"<head", b"<body", b"<?xml", b"<!--")


def _looks_like_page(head_bytes):
    return head_bytes.lstrip(b"\xef\xbb\xbf \t\r\n")[:20].lower().startswith(_PAGE_STARTS)


def content_problem(path, name=None):
    """None when the file's bytes can be what its name says it is; otherwise what is wrong."""
    name = (name or os.path.basename(str(path))).lower()
    try:
        with open(str(path), "rb") as fh:
            first = fh.read(4096)
    except OSError as err:
        return "unreadable ({})".format(err)
    page = _looks_like_page(first[:512])
    said = ""
    if page:   # what the server said instead, so the reason is visible, not guessed
        import re
        title = re.search(rb"(?is)<title[^>]*>(.*?)</title>", first)
        if title:
            said = ": " + re.sub(r"\s+", " ", title.group(1).decode("utf-8", "replace")).strip()[:150]
    for exts, magics in _MAGIC:
        if name.endswith(exts):
            if first.startswith(magics):
                return None
            what = "an HTML page (the server sent an error page{})".format(said) if page else \
                "starts with {!r}".format(first[:12])
            return "is not {} data: {}".format(exts[0].lstrip("."), what)
    if name.endswith((".html", ".htm", ".xml", ".svg")):
        return None
    if page:
        return "is an HTML page, not the file (the server sent an error page{})".format(said)
    return None


def _download_ranged(url, dest, size, parts, progress, retries, timeout):
    """Fetch [start, end] byte ranges on separate connections into one pre-sized file.

    Repositories cap each connection, not the total, so four ranges arrive about four
    times faster. Each range resumes from where it stopped if its connection drops.
    """
    import threading
    tmp = dest + ".part"
    with open(tmp, "wb") as fh:
        fh.truncate(size)
    step = -(-size // parts)
    ranges = [(i * step, min(size, (i + 1) * step) - 1) for i in range(parts) if i * step < size]
    done = [0] * len(ranges)
    errors = []
    lock = threading.Lock()

    def worker(k, start, end):
        pos = start
        for attempt in range(retries + 1):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": UA, "Range": "bytes={}-{}".format(pos, end)})
                with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "r+b") as out:
                    if resp.status != 206:
                        raise NetError("server ignored the byte range (HTTP {})".format(resp.status))
                    out.seek(pos)
                    while pos <= end:
                        chunk = resp.read(min(1 << 20, end - pos + 1))
                        if not chunk:
                            break
                        out.write(chunk)
                        pos += len(chunk)
                        with lock:
                            done[k] = pos - start
                if pos > end:
                    return
            except NetError as err:
                with lock:
                    errors.append(str(err))
                return
            except (urllib.error.URLError, socket.timeout, ConnectionError, OSError) as err:
                if attempt == retries:
                    with lock:
                        errors.append("{} for {}".format(getattr(err, "reason", err), url))
                    return
                time.sleep(min(30, 2 * (2 ** attempt)))
        with lock:
            errors.append("range {}-{} incomplete".format(start, end))

    threads = [threading.Thread(target=worker, args=(k, s, e)) for k, (s, e) in enumerate(ranges)]
    for th in threads:
        th.daemon = True
        th.start()
    last_report = time.time()
    while any(th.is_alive() for th in threads):
        time.sleep(0.5)
        if progress and time.time() - last_report > 15:
            with lock:
                progress(sum(done), size)
            last_report = time.time()
    for th in threads:
        th.join()
    if errors or sum(done) != size:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise NetError("parallel download failed: " + (errors[0] if errors else "{} of {} bytes".format(sum(done), size)))
    os.replace(tmp, dest)
    return size


def download(url, dest, max_bytes=None, progress=None, retries=3, timeout=120, size_hint=None):
    """Download to dest via dest.part; fail if the byte count doesn't match Content-Length,
    or if the bytes are not what the file's name says (an error page saved as .gz). Both are retried.

    Files over 16 MB on servers that accept byte ranges are fetched as HEGEMON_DOWNLOAD_PARTS
    (default 4) ranges at once; anything else, or any failure of that, uses one stream.
    size_hint: a size already known (from a listing), so small files skip the HEAD request.
    """
    dest = str(dest)
    name = os.path.basename(dest)
    if PARTS > 1 and (size_hint is None or size_hint >= SPLIT_ABOVE // 2):
        size, ranges = head(url)
        if ranges and size and size >= SPLIT_ABOVE and (max_bytes is None or size <= max_bytes):
            try:
                _download_ranged(url, dest, size, PARTS, progress, retries, timeout)
                problem = content_problem(dest, name)
                if problem is None:
                    return size
                os.unlink(dest)
            except NetError:
                pass  # fall back to one stream below
    tmp = dest + ".part"
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as out:
                total = resp.headers.get("Content-Length")
                total = int(total) if total else None
                if max_bytes and total and total > max_bytes:
                    raise NetError("{} is {:.1f} GB, over the limit".format(url, total / 1e9))
                written, last_report = 0, time.time()
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
                    written += len(chunk)
                    if max_bytes and written > max_bytes:
                        raise NetError("{} exceeded the size limit while downloading".format(url))
                    if progress and time.time() - last_report > 15:
                        progress(written, total)
                        last_report = time.time()
            if total is not None and written != total:
                raise NetError("incomplete download ({} of {} bytes)".format(written, total))
            problem = content_problem(tmp, name)
            if problem:
                raise NetError("{} {}".format(name, problem))
            os.replace(tmp, dest)
            return written
        except NetError as err:
            last = str(err)
            if "limit" in last:
                break
        except urllib.error.HTTPError as err:
            last = "HTTP {} for {}".format(err.code, url)
            if err.code not in (429, 500, 502, 503, 504):
                break
        except (urllib.error.URLError, socket.timeout, ConnectionError, OSError) as err:
            last = "{} for {}".format(getattr(err, "reason", err), url)
        if attempt < retries:
            time.sleep(min(60, 3 * (2 ** attempt)))
    if os.path.exists(tmp):
        os.unlink(tmp)
    raise NetError(last or "download failed: " + url)
