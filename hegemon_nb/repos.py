"""Search GEO + ArrayExpress and download a dataset's files.

GEO: the family SOFT file (what GEOparse.get_GEO downloads in the notebook), plus
processed supplementary files (GSE_RAW.tar extracted to GSE###_RAW/, like the
notebook asks) and NCBI-computed RNA-seq counts when GEO provides them.
ArrayExpress: SDRF/IDF and processed files (zips extracted, like counts/ in the
E-MTAB notebook).
"""
from collections import Counter, OrderedDict
import html
import os
import re
import shutil
import tarfile
import time
import zipfile
import zlib
from pathlib import Path
from urllib.parse import quote, urlencode, urljoin

from . import net, ui
from .geo_soft import read_geo_soft

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BIOSTUDIES = "https://www.ebi.ac.uk/biostudies/api/v1"
RAW_EXT = (".cel", ".cel.gz", ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".bai", ".sra", ".cram",
           ".bw", ".bigwig", ".wig", ".wig.gz", ".bed", ".bed.gz", ".idat", ".idat.gz", ".chp", ".chp.gz")
RAW_TYPES = {"CEL", "IDAT", "BAM", "BW", "BIGWIG", "BED", "WIG", "FASTQ", "SRA", "CHP"}


class FilesTooLarge(Exception):
    """Files the dataset needs were refused by the size limits: the data exists, the limits said no."""


class NotHosted(Exception):
    """The metadata names a processed file the repository does not actually serve."""


class SingleCell(Exception):
    """Single-cell study: no single-cell blueprint yet."""


SINGLE_CELL_WORDS = re.compile(
    r"single[- ]cell|single[- ]nucle|\bscrna|\bsnrna|\bscatac|cite-?seq|drop-?seq|10x genomics", re.I)
SINGLE_CELL_FILES = re.compile(
    r"(barcodes\.tsv|matrix\.mtx|features\.tsv|genes\.tsv|\.h5ad$|\.loom$|feature_bc_matrix|_cell_metadata)", re.I)


def looks_single_cell(title="", design="", files=(), extra=""):
    """Strong signals only: the title/design says so, or the files are cell-level matrices.
    A summary that merely cites single-cell work does not count (bulk studies often do)."""
    if SINGLE_CELL_WORDS.search(title or "") or SINGLE_CELL_WORDS.search(design or ""):
        return True
    if SINGLE_CELL_WORDS.search(extra or ""):
        return True
    return any(SINGLE_CELL_FILES.search(str(f)) for f in files)


def _probe(url):
    """(exists, size). exists is False only on a definite 404/410; anything unclear is 'try it'."""
    try:
        _, headers, _ = net.request(url, method="HEAD", timeout=20, retries=1)
        size = headers.get("Content-Length") or headers.get("content-length")
        return True, (int(size) if size and str(size).isdigit() else None)
    except net.NetError as err:
        if "HTTP 404" in str(err) or "HTTP 410" in str(err):
            return False, None
        return True, None


class TooFewSamples(Exception):
    def __init__(self, count, minimum):
        Exception.__init__(self, "{} samples, fewer than the minimum of {}".format(count, minimum))
        self.count, self.minimum = count, minimum


class TooManySamples(Exception):
    def __init__(self, count, maximum):
        Exception.__init__(self, "{} samples, more than --max-samples {}".format(count, maximum))
        self.count, self.maximum = count, maximum


class Dataset(object):
    def __init__(self, accession, repository, title="", summary="", kind="", organism="",
                 n_samples=None, platforms=None, note=""):
        self.accession = accession
        self.repository = repository  # "GEO" or "ArrayExpress"
        self.title = title
        self.summary = summary
        self.kind = kind
        self.organism = organism
        self.n_samples = n_samples
        self.platforms = platforms or []
        self.note = note
        self.catalog = OrderedDict()   # name -> CatalogEntry: every file the study offers (see catalog_add)
        self.failures = []             # downloads that did not work, shown to the AI

    def line(self):
        n = "{} samples".format(self.n_samples) if self.n_samples else ""
        bits = [b for b in (self.kind, self.organism, n) if b]
        return "{}  {}\n        {}".format(self.accession, self.title[:110], " | ".join(bits))


def detect_repository(accession, allow_other=False):
    """GEO and ArrayExpress are searched and downloaded automatically.

    Any other dataset is 'Other': you supply its files (--source) or its URLs (--url)
    and the notebook-plus-AI conversion treats it the same way.
    """
    acc = accession.strip().upper()
    if re.match(r"^GSE\d+$", acc):
        return "GEO"
    if re.match(r"^E-[A-Z]{4}-\d+$", acc):
        return "ArrayExpress"
    if allow_other:
        return "Other"
    raise ValueError(
        "{} is not a GEO series (GSE...) or an ArrayExpress accession (E-XXXX-...). Any other dataset works too, "
        "but you have to give it its files: --url <file url> (repeatable) or --source <folder>.".format(accession))


# ------------------------------------------------------------------ search ----

def _eutils(endpoint, params):
    params = dict(params, tool="hegemon-agent")
    time.sleep(0.35)  # NCBI allows ~3 requests/second without a key
    return net.get_json(EUTILS + endpoint + "?" + urlencode(params), timeout=30)


def _geo_summaries(uids):
    out = []
    for start in range(0, len(uids), 50):
        result = _eutils("esummary.fcgi", {"db": "gds", "id": ",".join(uids[start:start + 50]),
                                           "retmode": "json"}).get("result") or {}
        for uid in result.get("uids") or []:
            item = result.get(uid) or {}
            acc = str(item.get("accession") or "").upper()
            if not acc.startswith("GSE"):
                continue
            gpl = str(item.get("gpl") or "")
            out.append(Dataset(
                acc, "GEO", item.get("title") or acc, item.get("summary") or "",
                item.get("gdstype") or "", item.get("taxon") or "", item.get("n_samples"),
                ["GPL" + p for p in gpl.split(";") if p],
                note="no processed files listed" if not item.get("suppfile") and
                "sequencing" in (item.get("gdstype") or "") else ""))
    return out


def search_geo(query, limit=40):
    term = "({}) AND gse[ETYP]".format(query)
    ids = _eutils("esearch.fcgi", {"db": "gds", "term": term, "retmax": limit, "retmode": "json",
                                   "sort": "relevance"}).get("esearchresult", {}).get("idlist", [])
    found = _geo_summaries(ids) if ids else []
    # Only gene-expression series (the notebooks handle array + sequencing expression).
    return [d for d in found if "expression profiling" in d.kind.lower()]


def search_arrayexpress(query, limit=40):
    url = BIOSTUDIES + "/search?" + urlencode({"query": query, "collection": "arrayexpress",
                                               "pageSize": limit, "page": 1})
    out = []
    for hit in net.get_json(url, timeout=30).get("hits") or []:
        acc = str(hit.get("accession") or "").upper()
        if not acc.startswith("E-") or acc.startswith("E-GEOD-"):  # E-GEOD = copy of a GEO series
            continue
        out.append(Dataset(acc, "ArrayExpress", hit.get("title") or acc, hit.get("content") or "",
                           str(hit.get("type") or ""), ""))
    return out


# --------------------------------------------------------------- describe ----

def describe(accession):
    repo = detect_repository(accession)
    acc = accession.strip().upper()
    if repo == "GEO":
        try:
            ids = _eutils("esearch.fcgi", {"db": "gds", "term": "{}[ACCN] AND gse[ETYP]".format(acc),
                                           "retmode": "json"}).get("esearchresult", {}).get("idlist", [])
            hits = [d for d in _geo_summaries(ids) if d.accession == acc]
            if hits:
                return hits[0]
        except net.NetError as err:
            ui.warn("GEO summary lookup failed ({}); continuing with the SOFT file".format(err))
        return Dataset(acc, "GEO", acc)
    try:
        study = net.get_json("{}/studies/{}".format(BIOSTUDIES, quote(acc)), timeout=30)
        title = _attr(study, "Title") or _attr(study.get("section") or {}, "Title") or acc
        desc = _attr(study.get("section") or {}, "Description") or ""
        return Dataset(acc, "ArrayExpress", title, desc, _attr(study.get("section") or {}, "Study type") or "")
    except net.NetError as err:
        ui.warn("ArrayExpress lookup failed ({})".format(err))
        return Dataset(acc, "ArrayExpress", acc)


def _attr(node, name):
    for a in node.get("attributes") or []:
        if str(a.get("name", "")).lower() == name.lower():
            return a.get("value")
    return None


# --------------------------------------------------------------- download ----

class Budget(object):
    def __init__(self, cfg, assume_yes):
        import threading
        self.cfg = cfg
        self.used = 0
        self.assume_yes = assume_yes
        self.skipped = []
        self.lock = threading.Lock()

    def allow(self, name, size):
        with self.lock:
            gb = (size or 0) / 1e9
            over = gb > self.cfg.max_file_gb or (self.used + (size or 0)) / 1e9 > self.cfg.max_total_gb
            if over and not ui.yes_no("{} is {:.1f} GB (limits: {:.0f} GB/file, {:.0f} GB total). Download it?".format(
                    name, gb, self.cfg.max_file_gb, self.cfg.max_total_gb), False, self.assume_yes):
                self.skipped.append(name)
                return False
            return True

    def add(self, size):
        with self.lock:
            self.used += size or 0


WORKERS = int(os.environ.get("HEGEMON_DOWNLOAD_WORKERS", "4"))


def _parallel(items, fn, workers=None):
    """Run fn(item) for every item, a few at a time. Returns [(item, result, error)] in input order."""
    from concurrent.futures import ThreadPoolExecutor
    workers = max(1, min(workers or WORKERS, len(items) or 1))
    if workers == 1:
        out = []
        for item in items:
            try:
                out.append((item, fn(item), None))
            except Exception as err:
                out.append((item, None, err))
        return out
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(fn, item) for item in items]
        out = []
        for item, fut in zip(items, futures):
            try:
                out.append((item, fut.result(), None))
            except Exception as err:
                out.append((item, None, err))
        return out


def _fetch(url, dest, budget, size=None, quiet=False):
    dest = Path(dest)
    if dest.is_file() and (size is None or dest.stat().st_size == size):
        if not quiet:
            ui.info("have {} ({:.1f} MB)".format(dest.name, dest.stat().st_size / 1e6))
        return dest
    size = size if size is not None else net.remote_size(url)
    if not budget.allow(dest.name, size):
        return None
    free = shutil.disk_usage(str(dest.parent)).free
    if size and size + 2e9 > free:
        raise net.NetError("Not enough disk space for {} ({:.1f} GB free)".format(dest.name, free / 1e9))
    if not quiet:
        ui.info("downloading {}{}".format(dest.name, " ({:.1f} MB)".format(size / 1e6) if size else ""))
    cap = int(budget.cfg.max_file_gb * 1e9) if size is None else None

    def report(done, total):
        ui.info("  {} {:.0f} MB{}".format(dest.name, done / 1e6, " of {:.0f} MB".format(total / 1e6) if total else ""))

    written = net.download(url, dest, max_bytes=cap, progress=report, size_hint=size)
    budget.add(written)
    return dest


def _safe_under(root, relative):
    """Resolve a repository-supplied relative path inside root, or None if it escapes."""
    root = Path(root).resolve()
    cleaned = str(relative).replace("\\", "/").lstrip("/")
    if not cleaned or ".." in cleaned.split("/"):
        return None
    dest = (root / cleaned).resolve()
    if not str(dest).startswith(str(root) + os.sep):
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    return dest


def _safe_extract(archive, target):
    """Extract tar/zip without absolute paths, '..', links or devices."""
    target = Path(target).resolve()
    target.mkdir(parents=True, exist_ok=True)
    count = 0
    if zipfile.is_zipfile(str(archive)):
        with zipfile.ZipFile(str(archive)) as zf:
            for info in zf.infolist():
                dest = (target / info.filename).resolve()
                if not str(dest).startswith(str(target) + os.sep) or info.is_dir():
                    continue
                dest.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(info) as src, open(str(dest), "wb") as out:
                    shutil.copyfileobj(src, out)
                count += 1
        return count
    with tarfile.open(str(archive)) as tf:
        for member in tf.getmembers():
            dest = (target / member.name).resolve()
            if not member.isfile() or not str(dest).startswith(str(target) + os.sep):
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            src = tf.extractfile(member)
            with open(str(dest), "wb") as out:
                shutil.copyfileobj(src, out)
            count += 1
    return count


def _is_raw(name):
    low = name.lower()
    return low.endswith(RAW_EXT)


def geo_series_url(acc):
    num = acc[3:]
    stub = "GSE{}nnn".format(num[:-3]) if len(num) > 3 else "GSEnnn"
    return "https://ftp.ncbi.nlm.nih.gov/geo/series/{}/{}/".format(stub, acc)


def _listing(url):
    """Parse an NCBI FTP-over-HTTPS directory page into [(name, url, approx_bytes)]."""
    page = net.get_text(url, timeout=60)  # NetError propagates: "couldn't list" != "no files"
    out, seen = [], set()
    for m in re.finditer(r'<a href="([^"?/][^"]*)">[^<]*</a>\s+[\d-]+\s+[\d:]+\s+([\d.]+[KMGT]?)', page):
        name, size = html.unescape(m.group(1)), m.group(2)
        mult = {"K": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}.get(size[-1:], 1)
        try:
            approx = float(size.rstrip("KMGT")) * mult
        except ValueError:
            approx = None
        seen.add(name)
        out.append((name, urljoin(url, name), approx))
    for m in re.finditer(r'<a href="([^"?/#][^"]*)"', page):  # fallback if the page layout changes
        name = html.unescape(m.group(1))
        if name not in seen and not name.endswith("/") and "://" not in name:
            seen.add(name)
            out.append((name, urljoin(url, name), None))
    return out


def _raw_tar_contents(filelist_text):
    types = Counter()
    for line in filelist_text.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 5 and parts[0] == "File":
            types[parts[4].strip().upper()] += 1
    return types


def _ncbi_counts_links(acc):
    try:
        page = net.get_text("https://www.ncbi.nlm.nih.gov/geo/download/?acc=" + acc, timeout=60, retries=1)
    except net.NetError:
        return []
    links = set()
    for m in re.finditer(r'href="([^"]*type=rnaseq_counts[^"]*)"', page):
        link = html.unescape(m.group(1))
        if "raw_counts" in link or "annot" in link:
            links.add(urljoin("https://www.ncbi.nlm.nih.gov", link))
    return sorted(links)


def _https(url):
    """GEO's SOFT files give ftp:// URLs; the same paths serve over https."""
    url = str(url or "").strip()
    for prefix in ("ftp://ftp.ncbi.nlm.nih.gov/", "ftp://ftp.ebi.ac.uk/"):
        if url.startswith(prefix):
            return "https://" + url[len("ftp://"):]
    return url


def geo_sample_files(gse, platform=None):
    """Every GSM's own supplementary files, from the SOFT metadata.

    GEO serves each of these individually, so the GSE_RAW.tar bundle is never needed:
    fetching them into GSE###_RAW/ gives the exact layout the notebook untars to.
    """
    out = []
    for name, gsm in gse.gsms.items():
        if platform and (gsm.metadata.get("platform_id") or [""])[0] != platform:
            continue
        for key, values in gsm.metadata.items():
            if not key.startswith("supplementary_file"):
                continue
            for value in values:
                url = _https(value)
                if not url.lower().startswith("http") or url.upper().endswith("NONE"):
                    continue
                out.append((name, url.rsplit("/", 1)[-1], url))
    return out


def _sum_sizes(items, cap=40):
    """Ask for sizes only when there are few enough files that HEADs are cheap."""
    if len(items) > cap:
        return None
    total = 0
    for item in items:
        size = net.remote_size(item[-1])
        if size is None:
            return None
        total += size
    return total


def _show_plan(label, files, total):
    size = " (~{:.0f} MB)".format(total / 1e6) if total else ""
    ui.info("using {}: {} file(s){}".format(label, len(files), size))


class Incomplete(net.NetError):
    """Some files of a source could not be downloaded even one at a time."""


def _download_all(files, source, budget, subdir=None, quiet=None, require_all=False, errors_out=None):
    """files: [(name, url, size)]. Several at once when there are several. Returns names written.

    Anything that fails while several downloads run at once is tried again on its own, because
    some servers refuse parallel requests. require_all: raise Incomplete if a file still fails,
    so a source is never used with pieces missing.
    """
    import threading
    quiet = len(files) > 1 if quiet is None else quiet
    target = source / subdir if subdir else source
    done, lock = [0], threading.Lock()

    def one(item):
        name, url, size = item
        dest = _safe_under(target, name)
        if dest is None:
            return None
        path = _fetch(url, dest, budget, size, quiet=quiet)
        with lock:
            done[0] += 1
            if quiet and len(files) > 1 and (done[0] % 10 == 0 or done[0] == len(files)):
                ui.info("{} of {} files".format(done[0], len(files)))
        if path is not None and name.lower().endswith((".tar", ".tar.gz", ".tgz", ".zip")):
            try:
                n = _safe_extract(path, path.parent / re.sub(r"(\.tar\.gz|\.tgz|\.tar|\.zip)$", "", path.name))
                ui.info("{} extracted ({} files)".format(name, n))
            except Exception as err:
                # Keep the file: the AI can still be told about it, and one bad
                # archive must not end the build with a traceback.
                ui.warn("{} could not be unpacked ({}); leaving it as it is".format(name, err))
        return path

    if len(files) > 1:
        ui.info("downloading {} files, {} at a time".format(len(files), min(WORKERS, len(files))))
    results = _parallel(files, one)
    failed = [item for item, _path, err in results if err is not None]
    if failed and len(files) > 1 and WORKERS > 1:
        ui.info("{} file(s) failed while downloading in parallel; trying them one at a time".format(len(failed)))
        time.sleep(3)
        again = {item: (path, err) for item, path, err in _parallel(failed, one, workers=1)}
        results = [(item, again[item][0], again[item][1]) if item in again else (item, path, err)
                   for item, path, err in results]
    got, errors = [], []
    for (name, _url, _size), path, err in results:
        if err is not None:
            ui.warn("could not download {}: {}".format(name, err))
            errors.append("{}: {}".format(name, err))
        elif path is not None:
            got.append(str(Path(path).resolve().relative_to(Path(source).resolve())))
    if errors_out is not None:
        errors_out.extend(errors)
    if errors and require_all:
        raise Incomplete("{} of {} file(s) could not be downloaded ({})".format(
            len(errors), len(files), "; ".join(errors[:3])))
    return got


# ------------------------------------------------------------------ catalog ----
#
# Every study's list of the files it offers. The fixed plans below download the cheapest set
# that works, fast and with no AI. Anything else in the catalog, or any public URL, can be
# asked for by the AI (FETCH), so a study with unusual files, a failed download, or a site this
# agent has never seen is handled by the same loop instead of by site-specific code.

class CatalogEntry(object):
    __slots__ = ("name", "url", "size", "note")

    def __init__(self, name, url, size=None, note=""):
        self.name, self.url, self.size, self.note = name, url, size, note


def catalog_add(ds, name, url, size=None, note=""):
    name = str(name).replace("\\", "/").lstrip("/")
    if name and url and name not in ds.catalog:
        ds.catalog[name] = CatalogEntry(name, url, size, note)


def _public_url_problem(url):
    """None for an http(s) URL on a public host; otherwise why it may not be fetched."""
    import ipaddress
    import socket
    from urllib.parse import urlparse
    parts = urlparse(url)
    if parts.scheme not in ("http", "https") or not parts.hostname:
        return "only http and https links can be fetched"
    try:
        addresses = {info[4][0] for info in socket.getaddrinfo(parts.hostname, None)}
    except OSError as err:
        return "the host does not resolve ({})".format(err)
    for address in addresses:
        try:
            if not ipaddress.ip_address(address.split("%")[0]).is_global:
                return "it points to a private or local address"
        except ValueError:
            return "unreadable address {}".format(address)
    return None


def _name_for_url(url):
    from urllib.parse import parse_qs, unquote, urlparse
    parts = urlparse(url)
    for key in ("file", "filename", "name"):
        vals = parse_qs(parts.query).get(key)
        if vals and vals[0].strip():
            return re.sub(r"[^A-Za-z0-9._-]", "_", unquote(vals[0].strip().rsplit("/", 1)[-1]))
    base = unquote(parts.path.rstrip("/").rsplit("/", 1)[-1]) or parts.hostname or "download"
    return re.sub(r"[^A-Za-z0-9._-]", "_", base)


def _visible_text(page):
    text = re.sub(r"(?is)<(script|style|noscript)[^>]*>.*?</\1>", " ", page)
    text = re.sub(r"(?s)<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def page_catalog(ds, url, page, limit=400):
    """A web page that is not itself a data file: its links become the catalog, its text a note."""
    title = re.search(r"(?is)<title[^>]*>(.*?)</title>", page)
    added = 0
    for m in re.finditer(r"""(?is)<a\s[^>]*?href\s*=\s*["']([^"'#]+)["'][^>]*>(.*?)</a>""", page):
        href, label = html.unescape(m.group(1)).strip(), _visible_text(m.group(2))[:80]
        if href.lower().startswith(("mailto:", "javascript:", "tel:")):
            continue
        link = _https(urljoin(url, href))
        if not link.lower().startswith("http"):
            continue
        name = _name_for_url(link)
        if name in ds.catalog:
            if ds.catalog[name].url == link:
                continue
            name = "{}-{}".format(added + 1, name)
        catalog_add(ds, name, link, None, "link on the study page" + (": " + label if label else ""))
        added += 1
        if added >= limit:
            break
    about = "Study page {}{}: {}".format(url, " ({})".format(_visible_text(title.group(1))[:120]) if title else "",
                                         _visible_text(page)[:2500])
    return added, about


def fetch_requested(ds, source, cfg, items, assume_yes=False):
    """Download what the AI asked for: catalog names, name patterns (* ?), or public URLs.
    Returns (report lines for the AI, names downloaded)."""
    import fnmatch
    source = Path(source)
    budget = Budget(cfg, assume_yes)
    report, chosen = [], OrderedDict()
    for raw in items:
        item = raw.strip().strip("`'\"").strip()
        if not item:
            continue
        if item in ds.catalog:
            matches = [ds.catalog[item]]
        elif any(ch in item for ch in "*?["):
            matches = [e for n, e in ds.catalog.items()
                       if fnmatch.fnmatch(n, item) or fnmatch.fnmatch(n.rsplit("/", 1)[-1], item)]
        else:
            matches = [e for n, e in ds.catalog.items() if n.rsplit("/", 1)[-1] == item]
            if not matches and re.match(r"^(https?|ftp)://", item, re.I):
                url = _https(item)
                problem = _public_url_problem(url)
                if problem:
                    report.append("{}: not fetched, {}".format(item, problem))
                    continue
                matches = [CatalogEntry(_name_for_url(url), url, None, "asked for by URL")]
        if not matches:
            report.append("{}: not in the list of files this study offers (ask by exact name, pattern or URL)"
                          .format(item))
            continue
        if len(matches) > 1000:
            report.append("{}: matches {} files; ask for fewer".format(item, len(matches)))
            continue
        for entry in matches:
            chosen.setdefault(entry.name, entry)
    files = []
    for name, entry in chosen.items():
        if (source / name).is_file() and (source / name).stat().st_size > 0:
            report.append("{}: already downloaded".format(name))
        else:
            files.append((name, entry.url, int(entry.size) if entry.size else None))
    errors, got = [], []
    if files:
        ui.info("downloading {} file(s) the AI asked for".format(len(files)))
        got = _download_all(files, source, budget, errors_out=errors)
    for name in got:
        path = source / name
        report.append("{}: downloaded ({:,} bytes)".format(name, path.stat().st_size if path.is_file() else 0))
    report += ["{}: FAILED ({})".format(*e.split(": ", 1)) if ": " in e else e for e in errors]
    if budget.skipped:
        report.append("not downloaded because of the size limits: " + ", ".join(budget.skipped))
    ds.failures += errors
    return report, got



# ------------------------------------------------------------------ choosing files ----
#
# Some studies offer the same data several times (raw counts, normalized counts, TPM, per
# tissue, per cohort). Downloading all of it can take minutes. When several files add up to
# a lot, their first lines are read with one small request each, and the AI picks the ones
# the conversion needs. The rest stay in the catalog in case the conversion asks for them.

PICK_ABOVE = float(os.environ.get("HEGEMON_PICK_ABOVE_MB", "100")) * 1e6


def _peek(url, limit=16384):
    """The first lines of a remote file, from one small byte-range request (no download)."""
    try:
        _, _, body = net.request(url, headers={"Range": "bytes=0-{}".format(limit - 1)}, timeout=30,
                                 retries=1, max_bytes=limit)
    except net.NetError:
        return ""
    if body[:2] == b"\x1f\x8b":
        try:
            body = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(body)
        except zlib.error:
            return "(gzip data)"
    elif body[:2] == b"PK" or b"\x00" in body[:1024]:
        return "(binary file)"
    lines = body.decode("utf-8", "replace").splitlines()[:4]
    return "\n".join(line[:400] for line in lines)


def _pick(ds, candidates, picker, what, minimum=None):
    """candidates: [(name, url, size)]. Returns the ones to download now."""
    minimum = PICK_ABOVE if minimum is None else minimum
    if picker is None or len(candidates) < 2 or len(candidates) > 40:
        return candidates
    if any(size is None for _, _, size in candidates):   # sizes decide whether choosing is worth it
        sized = _parallel(candidates, lambda it: it[2] if it[2] is not None else net.remote_size(it[1]))
        candidates = [(item[0], item[1], size if err is None else None) for item, size, err in sized]
    total = sum(size or 0 for _, _, size in candidates)
    if total < minimum:
        return candidates
    ui.info("{} files ({:.0f} MB) could be downloaded; reading their first lines so the AI picks what is needed"
            .format(len(candidates), total / 1e6))
    heads = {item[0]: text for item, text, err in _parallel(candidates, lambda it: _peek(it[1])) if err is None}
    try:
        chosen, why = picker(ds, [(name, size, heads.get(name, "")) for name, _url, size in candidates], what)
    except Exception as err:   # choosing is an optimization; if it fails, download them all as before
        ui.warn("could not choose between the files ({}); downloading all of them".format(err))
        return candidates
    keep = [c for c in candidates if c[0] in set(chosen)]
    if not keep:
        ui.warn("the AI chose none of the files; downloading all of them")
        return candidates
    kept = sum(size or 0 for _, _, size in keep)
    ui.info("the AI chose {} of {} files ({:.0f} of {:.0f} MB): {}{}".format(
        len(keep), len(candidates), kept / 1e6, total / 1e6, ", ".join(c[0] for c in keep[:6]),
        " ({})".format(why) if why else ""))
    for name, url, size in candidates:
        if (name, url, size) not in keep:
            catalog_add(ds, name, url, size, "not downloaded: the AI chose other files first")
    return keep


def _pick_kinds(ds, per_sample, picker):
    """Per-sample files of several kinds (counts, FPKM, ...): the AI picks the kinds. per_sample: [(name, url, size)]."""
    if picker is None:
        return per_sample
    kinds = OrderedDict()
    for item in per_sample:
        kinds.setdefault(re.sub(r"^GSM\d+", "GSM#", item[0].rsplit("/", 1)[-1]), []).append(item)
    if len(kinds) < 2 or min(len(v) for v in kinds.values()) < 2:
        return per_sample
    examples = [("{} ({} files)".format(kind, len(items)), None, items[0][1]) for kind, items in kinds.items()]
    heads = {e[0]: text for e, text, err in _parallel(examples, lambda e: _peek(e[2])) if err is None}
    try:
        chosen, why = picker(ds, [(label, None, heads.get(label, "")) for label, _s, _u in examples],
                             "kinds of per-sample files (one file of each kind per sample)")
    except Exception as err:
        ui.warn("could not choose between the per-sample file kinds ({}); downloading all of them".format(err))
        return per_sample
    keep_kinds = [kind for kind in kinds if "{} ({} files)".format(kind, len(kinds[kind])) in set(chosen)]
    if not keep_kinds:
        return per_sample
    ui.info("the AI chose these per-sample files: {}{}".format(", ".join(keep_kinds), " ({})".format(why) if why else ""))
    return [item for kind in keep_kinds for item in kinds[kind]]



def fetch_geo(ds, source, cfg, assume_yes=False, platform=None, min_samples=0, max_samples=0, picker=None):
    """Download the smallest set of files that reproduces the notebook's inputs.

    Order of preference, cheapest first, matching what an operator picks by hand:
      1. GSM data tables inside the SOFT file          (notebook microarray path)
      2. NCBI-computed counts for a sequencing series  (notebook "single file" path)
      3. series-level processed files that are not archives
      4. each GSM's own supplementary files, into GSE###_RAW/   (notebook "separate files" path)
      5. an archive bundle such as GSE###_RAW.tar      (last resort; same content as 4)
    Everything the series offers also goes into ds.catalog, and anything that failed into
    ds.failures, so the AI can ask for other files instead of the build stopping here.
    Returns (soft_path, gse, notes, processed_sources); gse has every sample's metadata but only the
    platform tables and the first two sample tables per platform.
    """
    source = Path(source)
    source.mkdir(parents=True, exist_ok=True)
    budget = Budget(cfg, assume_yes)
    base = geo_series_url(ds.accession)
    soft = _fetch(base + "soft/{}_family.soft.gz".format(ds.accession),
                  source / "{}_family.soft.gz".format(ds.accession), budget)
    if soft is None:
        raise net.NetError("The family SOFT file is required (it holds the sample metadata).")
    gse = read_geo_soft(soft, sample_tables=2)   # all metadata + what the evidence shows; one quick pass
    if min_samples:
        per_platform = Counter((g.metadata.get("platform_id") or ["?"])[0] for g in gse.gsms.values())
        largest = max(per_platform.values()) if per_platform else 0
        if largest < min_samples:
            raise TooFewSamples(largest, min_samples)  # decided from the SOFT file, before any data download
    if max_samples:
        per_platform = Counter((g.metadata.get("platform_id") or ["?"])[0] for g in gse.gsms.values())
        largest = max(per_platform.values()) if per_platform else 0
        if largest > max_samples:
            raise TooManySamples(largest, max_samples)
    types = " ".join(gse.metadata.get("type", [])).lower()
    if looks_single_cell(" ".join(gse.metadata.get("title", [])), " ".join(gse.metadata.get("overall_design", [])),
                         [name for _g, name, _u in geo_sample_files(gse)]):
        raise SingleCell("{} is a single-cell study (from its title, design or per-cell files)".format(ds.accession))
    notes, processed, too_large, failures = [], [], [], []

    # what the series offers (the AI can ask for any of it later)
    per_sample = geo_sample_files(gse, platform)
    for gsm, name, url in per_sample:
        catalog_add(ds, "{}_RAW/{}".format(ds.accession, name), url, None,
                    "raw per-sample file" if _is_raw(name) else "per-sample file of " + gsm)
    try:
        listing = [(n, u, sz) for n, u, sz in _listing(base + "suppl/") if n != "filelist.txt"]
    except net.NetError as err:
        listing = []
        if "HTTP 404" in str(err):
            notes.append("The series has no supplementary folder on GEO.")
        else:
            failures.append("could not list {}suppl/: {}".format(base, err))
    for name, url, size in listing:
        low = name.lower()
        catalog_add(ds, name, url, size, "raw data" if _is_raw(name) else
                    "archive" if low.endswith((".tar", ".tar.gz", ".tgz", ".zip")) else "series file")

    # 1. tables already in the SOFT file
    rows = [int((g.metadata.get("data_row_count") or ["0"])[0] or 0) for g in gse.gsms.values()]
    with_tables = sum(1 for r in rows if r > 0)
    if with_tables == len(rows) and with_tables and "sequencing" not in types:
        notes.append("Sample tables are in the SOFT file (the notebook's microarray path); nothing else downloaded.")
        ds.failures += failures
        return soft, gse, notes, ["GSM data tables in the SOFT file ({} samples)".format(with_tables)]

    # 2. NCBI-computed counts: a few MB, and the notebook's single-file path
    if "sequencing" in types:
        links = _ncbi_counts_links(ds.accession)
        if links:
            files = []
            for link in links:
                match = re.search(r"file=([^&]+)", link)
                name = match.group(1) if match else link.split("?", 1)[0].rsplit("/", 1)[-1]
                files.append((name, link, None))
                catalog_add(ds, name, link, None, "GEO's computed counts")
            _show_plan("GEO's own computed counts", files, _sum_sizes([(f[1],) for f in files]))
            try:
                got = _download_all(files, source, budget, require_all=True)
            except Incomplete as err:
                # counts without their annotation file (or the reverse) cannot be converted: use the next source
                got = []
                for name, _url, _size in files:
                    leftover = source / name
                    if leftover.is_file():
                        leftover.unlink()
                failures.append("GEO's computed counts: {}".format(err))
                notes.append("GEO's computed counts could not be downloaded ({}); using the series' own files."
                             .format(err))
                ui.warn("GEO's computed counts could not be downloaded; trying the series' own files")
            if got:
                processed += got
                notes.append("Used GEO's NCBI-computed counts ({}); the supplementary archives were not needed."
                             .format(", ".join(got)))
                return soft, gse, notes, processed
        else:
            notes.append("GEO has no NCBI-computed counts for this series.")

    def small_enough(size):
        return size is None or size <= cfg.max_file_gb * 1e9

    # 3. series-level processed files that are not archives
    archives, plain = [], []
    for name, url, size in listing:
        if _is_raw(name):
            continue
        (archives if name.lower().endswith((".tar", ".tar.gz", ".tgz", ".zip")) else plain).append((name, url, size))
    plain = _pick(ds, plain, picker, "series-level processed files")
    if plain:
        _show_plan("the series' processed files", plain, sum((s or 0) for _, _, s in plain) or None)
        errors = []
        got = _download_all(plain, source, budget, errors_out=errors)
        processed += got
        if got:
            notes.append("Series-level processed files: {}".format(", ".join(got)))
        if errors:
            failures.append("series files: " + "; ".join(errors[:5]))

    # 4. per-sample files, fetched individually into GSE###_RAW/ (no bundle)
    need_bundle = False
    if not processed:
        keep = [(name, url, None) for _gsm, name, url in per_sample if not _is_raw(name)]
        skipped_raw = len(per_sample) - len(keep)
        keep = _pick_kinds(ds, keep, picker)
        if keep:
            total = _sum_sizes([(k[1],) for k in keep])
            _show_plan("each sample's own supplementary file", keep, total)
            ui.info("these are the same files GSE{}_RAW.tar contains, so the bundle is skipped"
                    .format(ds.accession[3:]))
            if total is not None and total > cfg.max_total_gb * 1e9 and not ui.yes_no(
                    "That is {:.1f} GB in total. Download it?".format(total / 1e9), False,
                    True if assume_yes else None):
                too_large.append("{} per-sample files totalling {:.1f} GB".format(len(keep), total / 1e9))
            else:
                errors = []
                got = _download_all(keep, source, budget, subdir="{}_RAW".format(ds.accession), errors_out=errors)
                processed += got
                if got:
                    notes.append("Fetched {} per-sample file(s) into {}_RAW/ (the layout the notebook untars to); "
                                 "the {}_RAW.tar bundle was not downloaded."
                                 .format(len(got), ds.accession, ds.accession))
                if errors:
                    # the bundle holds the same files; what did arrive is kept either way
                    need_bundle = True
                    failures.append("per-sample files: {} of {} could not be downloaded ({})".format(
                        len(errors), len(keep), "; ".join(errors[:3])))
                    notes.append("{} of {} per-sample files could not be downloaded; trying the {}_RAW bundle for "
                                 "them.".format(len(errors), len(keep), ds.accession))
        if skipped_raw:
            notes.append("Skipped {} raw per-sample file(s) (sequence/alignment/array formats).".format(skipped_raw))

    # 5. archive bundles, only if nothing above produced data (or per-sample files are missing)
    if archives and (not processed or need_bundle):
        wanted = [a for a in archives if "_raw" in a[0].lower()] if processed else archives
        for name, url, size in wanted:
            if not small_enough(size):
                too_large.append("{} ({:.1f} GB)".format(name, (size or 0) / 1e9))
                continue
            _show_plan("the archive " + name, [(name, url, size)], size)
            errors = []
            got = _download_all([(name, url, size)], source, budget, errors_out=errors)
            processed += [g for g in got if g not in processed]
            if got:
                notes.append("{} was {}.".format(name, "used for the missing per-sample files" if need_bundle
                                                 else "the only processed source available"))
            if errors:
                failures.append(errors[0])
    if failures:
        ds.failures += failures
        notes.append("Download problems (other files can be asked for): " + " | ".join(failures)[:1500])
    if too_large or budget.skipped:
        notes.append("NOT downloaded because of size: " + "; ".join(too_large + budget.skipped) +
                     ". Raise HEGEMON_MAX_FILE_GB / HEGEMON_MAX_TOTAL_GB, or pass --yes, to include them.")
    if with_tables:
        processed.append("GSM data tables in the SOFT file ({} samples)".format(with_tables))
    return soft, gse, notes, processed


def _is_page(url):
    """True when the URL answers with a web page rather than a file."""
    try:
        _, headers, body = net.request(url, timeout=60, retries=1, max_bytes=4096)
    except net.NetError:
        return False
    kind = str({k.lower(): v for k, v in headers.items()}.get("content-type", "")).lower()
    return "text/html" in kind or net._looks_like_page(body[:512])


def fetch_urls(ds, source, cfg, urls, assume_yes=False):
    """Any dataset from any site: each URL is either a file (downloaded and unpacked) or a web page
    (its links become the study's catalog, its text a note), so the AI can pick the files.

    The conversion itself is repository-agnostic - the notebook blueprint plus the AI
    handle whatever files are here - so this is all that a new repository needs.
    """
    source = Path(source)
    source.mkdir(parents=True, exist_ok=True)
    budget = Budget(cfg, assume_yes)
    files, notes, pages = [], [], 0
    for url in urls:
        url = _https(url)
        name = _name_for_url(url)
        if not name.lower().endswith((".gz", ".zip", ".tar", ".tgz", ".bz2", ".xz", ".txt", ".tsv", ".csv",
                                      ".xlsx", ".xls", ".h5", ".h5ad", ".rds", ".soft", ".json")) and _is_page(url):
            page = net.get_text(url, timeout=60, max_bytes=5 * 1024 * 1024)
            added, about = page_catalog(ds, url, page)
            notes.append(about)
            ui.info("{} is a web page: {} link(s) on it are offered to the AI".format(url, added))
            pages += 1
            continue
        catalog_add(ds, name, url, None, "URL you gave")
        files.append((name, url, net.remote_size(url)))
    got = []
    if files:
        _show_plan("the URLs you gave", files, sum((s or 0) for _, _, s in files) or None)
        errors = []
        got = _download_all(files, source, budget, quiet=False, errors_out=errors)
        ds.failures += errors
    if not got and not pages:
        raise net.NetError("None of the {} URL(s) could be downloaded.".format(len(urls)))
    if got:
        notes.append("Downloaded {} file(s) from the URLs given.".format(len(got)))
    return notes, got


def _biostudies_files(acc):
    info = net.get_json("{}/studies/{}/info".format(BIOSTUDIES, quote(acc)), timeout=30)
    base = str(info.get("httpLink") or info.get("ftpLink") or "").rstrip("/")
    if base.startswith("ftp://"):
        base = "https://" + base[len("ftp://"):]
    items = []
    for url in ("{}/files/{}?pageSize=1000&page=0".format(BIOSTUDIES, quote(acc)),
                "{}/studies/{}/files?limit=1000&offset=0".format(BIOSTUDIES, quote(acc))):
        try:
            payload = net.get_json(url, timeout=60, retries=1)
        except net.NetError:
            continue
        rows = payload if isinstance(payload, list) else next(
            (payload.get(k) for k in ("files", "items", "list", "data") if isinstance(payload.get(k), list)), [])
        items = [r for r in rows if isinstance(r, dict)]
        if items:
            break
    if not items:
        study = net.get_json("{}/studies/{}".format(BIOSTUDIES, quote(acc)), timeout=60)
        items = list(_walk_files(study))
    out = []
    for item in items:
        path = item.get("path") or item.get("name") or item.get("fileName")
        if not path or str(item.get("type", "")).lower() == "directory":
            continue
        url = "{}/Files/{}".format(base, quote(path)) if base else \
            "https://www.ebi.ac.uk/biostudies/files/{}/{}".format(acc, quote(path))
        out.append((path, url, item.get("size")))
    return out


def _walk_files(node):
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "files":
                for f in value:
                    for leaf in (f if isinstance(f, list) else [f]):
                        if isinstance(leaf, dict):
                            yield leaf
            else:
                for leaf in _walk_files(value):
                    yield leaf
    elif isinstance(node, list):
        for value in node:
            for leaf in _walk_files(value):
                yield leaf


def ae_files_base(acc):
    """The study's Files/ folder, e.g. https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/604/E-MTAB-7604/Files/"""
    try:
        info = net.get_json("{}/studies/{}/info".format(BIOSTUDIES, quote(acc)), timeout=30, retries=1)
        base = str(info.get("httpLink") or info.get("ftpLink") or "").rstrip("/")
        if base.startswith("ftp://"):
            base = "https://" + base[len("ftp://"):]
        if base:
            return base + "/Files/"
    except (net.NetError, ValueError):
        pass
    m = re.match(r"^(E-[A-Z]+-)(\d+)$", acc)
    if not m:
        raise net.NetError("Not an ArrayExpress accession: " + acc)
    return "https://ftp.ebi.ac.uk/biostudies/fire/{}/{}/{}/Files/".format(m.group(1), m.group(2)[-3:].zfill(3), acc)


def sdrf_data_files(path):
    """The processed per-sample file names the SDRF lists.

    MAGE-TAB keeps raw files in 'Array Data File' and processed ones in
    'Derived Array Data File' / 'Derived Array Data Matrix File' / 'Processed Data File',
    so only the derived/processed columns are used here.
    """
    import pandas as pd
    table = pd.read_csv(str(path), sep="\t", dtype=str, keep_default_na=False)
    names = []
    for col in table.columns:
        low = col.lower()
        if "data" not in low or "file" not in low:
            continue
        if "derived" in low or "processed" in low or "normalized" in low or "normalised" in low:
            names += [v.strip() for v in table[col] if v.strip()]
    return list(dict.fromkeys(names))


def fetch_arrayexpress(ds, source, cfg, assume_yes=False, min_samples=0, max_samples=0, picker=None):
    """Returns (notes, processed_files). The SDRF is required, as is every file it names."""
    acc = ds.accession
    source = Path(source)
    source.mkdir(parents=True, exist_ok=True)
    budget = Budget(cfg, assume_yes)
    notes, tried = [], []
    base = ae_files_base(acc)
    ui.info("study folder: " + base)
    sdrf = None
    for url in (base + acc + ".sdrf.txt",
                "https://www.ebi.ac.uk/biostudies/files/{0}/{0}.sdrf.txt".format(acc)):
        try:
            sdrf = _fetch(url, source / (acc + ".sdrf.txt"), budget)
            break
        except net.NetError as err:
            tried.append("{} ({})".format(url, err))
    if sdrf is None:
        raise net.NetError("Could not download {}.sdrf.txt, which holds ArrayExpress's sample metadata. Not "
                           "continuing without it.\n  tried: ".format(acc) + "\n  tried: ".join(tried))
    ui.ok("sample metadata: " + sdrf.name)
    n_samples = None
    if min_samples or max_samples or picker:
        import pandas as pd
        table = pd.read_csv(str(sdrf), sep="\t", dtype=str, keep_default_na=False)
        col = next((c for c in table.columns if c.lower().startswith("source name")), table.columns[0])
        n_samples = count = table[col].str.split(".").str[0].nunique()
        if min_samples and count < min_samples:
            raise TooFewSamples(count, min_samples)
        if max_samples and count > max_samples:
            raise TooManySamples(count, max_samples)
    try:
        _fetch(base + acc + ".idf.txt", source / (acc + ".idf.txt"), budget)
    except net.NetError:
        notes.append("No IDF file (optional).")

    idf_text = ""
    idf = source / (acc + ".idf.txt")
    if idf.is_file():
        idf_text = idf.read_text(errors="replace")
    sdrf_text = sdrf.read_text(errors="replace")
    wanted = sdrf_data_files(sdrf)

    # What the study folder holds: the BioStudies file API (the only one that sees inside
    # subfolders such as counts/) merged with the FTP directory listing.
    catalog = OrderedDict()
    try:
        for path, url, size in _biostudies_files(acc):
            catalog.setdefault(path, (url, size))
    except (net.NetError, ValueError, KeyError) as err:
        notes.append("BioStudies file list unavailable ({}); using the FTP listing".format(err))
    try:
        for name, url, size in _listing(base):
            if not name.endswith("/"):
                catalog.setdefault(name, (url, size))
    except net.NetError as err:
        notes.append("Could not list {} ({})".format(base, err))
    catalog.pop(acc + ".sdrf.txt", None)
    catalog.pop(acc + ".idf.txt", None)
    for name, (url, size) in catalog.items():
        low = name.lower()
        catalog_add(ds, name, url, int(size) if size else None,
                    "raw data" if _is_raw(name) or ".raw." in low else
                    "archive" if low.endswith((".zip", ".tar", ".tar.gz", ".tgz")) else "study file")

    sc_columns = [c for c in sdrf_text.split("\n", 1)[0].split("\t")
                  if "single cell" in c.lower() or "library construction" in c.lower()]
    sc_values = " ".join(v for v in re.findall(r"10x|drop-?seq|chromium", sdrf_text, re.I)[:3])
    if looks_single_cell(ds.title, idf_text[:20000], list(catalog) + wanted,
                         extra=" ".join(c for c in sc_columns if "single cell" in c.lower()) or
                         ("10x genomics" if sc_values else "")):
        raise SingleCell("{} is a single-cell study (from its IDF, SDRF or per-cell files)".format(acc))

    got, too_large, skipped_raw = [], [], 0
    import threading
    fetch_lock = threading.Lock()

    def fetch(name, url, size, quiet=False):
        dest = _safe_under(source, name)
        if dest is None:
            notes.append("Refused an unsafe file path from the repository: " + name)
            return None
        path = _fetch(url, dest, budget, int(size) if size else None, quiet=quiet)
        with fetch_lock:
            if path is None:
                too_large.append(name)
                return None
            got.append(name)
        if name.lower().endswith((".zip", ".tar", ".tar.gz", ".tgz")):
            try:
                n = _safe_extract(path, path.parent / re.sub(r"(\.tar\.gz|\.tgz|\.tar|\.zip)$", "", path.name))
                notes.append("{} extracted ({} files).".format(name, n))
            except Exception as err:
                ui.warn("{} could not be unpacked ({}); leaving it as it is".format(name, err))
        return path

    def present():
        files = [q for q in source.rglob("*") if q.is_file()]
        return {str(q.relative_to(source)) for q in files} | {q.name for q in files}

    def missing_now():
        have = present()
        return [n for n in wanted if n not in have and n.rsplit("/", 1)[-1] not in have]

    if wanted:
        # The SDRF says which files hold the processed data. Fetch those and nothing else,
        # the way an operator reads the SDRF and downloads just that.
        ui.info("the SDRF names {} processed file(s); fetching only those".format(len(wanted)))
        loose, not_loose = [], []
        for name in wanted:
            if name in catalog:
                loose.append((name,) + catalog[name])
                continue
            url = base + quote(name)
            exists, size = _probe(url)
            (loose if exists else not_loose).append((name, url, size))
        if n_samples and len(wanted) < n_samples / 2.0:
            # a few big matrices, not one file per sample: they may be alternatives of each other
            chosen = _pick(ds, loose, picker, "processed data files named in the SDRF")
            if len(chosen) < len(loose):
                dropped = {c[0] for c in loose} - {c[0] for c in chosen}
                wanted = [w for w in wanted if w not in dropped]
                loose = chosen
        if len(loose) > 1:
            ui.info("downloading {} files, {} at a time".format(len(loose), min(WORKERS, len(loose))))

        def gone(err):
            return isinstance(err, net.NetError) and ("HTTP 404" in str(err) or "HTTP 410" in str(err))
        results = _parallel(loose, lambda it: fetch(it[0], it[1], it[2], quiet=len(loose) > 1))
        not_loose += [item for item, _path, err in results if err is not None and gone(err)]
        retry = [(item, err) for item, _path, err in results if err is not None and not gone(err)]
        if retry and len(loose) > 1 and WORKERS > 1:
            # some servers refuse parallel requests: what failed that way is tried on its own
            ui.info("{} file(s) failed while downloading in parallel; trying them one at a time".format(len(retry)))
            time.sleep(3)
            again = _parallel([item for item, _err in retry], lambda it: fetch(it[0], it[1], it[2]), workers=1)
            not_loose += [item for item, _path, err in again if err is not None and gone(err)]
            retry = [(item, err) for item, _path, err in again if err is not None and not gone(err)]
        if retry:
            raise retry[0][1]
        if missing_now():
            # Older studies keep their processed files inside the study's archives.
            archives = [(n, u, s) for n, (u, s) in catalog.items()
                        if n.lower().endswith((".zip", ".tar", ".tar.gz", ".tgz")) and not _is_raw(n)
                        and not re.search(r"(^|[._-])raw([._-]|$)", n.lower().rsplit("/", 1)[-1])]
            archives.sort(key=lambda a: (0 if re.search(r"process|derived|normali[sz]", a[0], re.I) else 1,
                                         a[2] or 0))
            for name, url, size in archives:
                if not missing_now():
                    break
                ui.info("looking for {} inside {}".format(", ".join(missing_now()[:3]), name))
                fetch(name, url, size)
        still = missing_now()
        if still:
            if too_large:
                raise FilesTooLarge("{} needs {} but the size limits refused: {}".format(
                    acc, ", ".join(still[:5]), ", ".join(too_large)))
            problem = ("The SDRF names {} as the processed data, but ArrayExpress does not serve {} (checked the "
                       "study folder{}).".format(", ".join(still[:5]), "it" if len(still) == 1 else "them",
                                                " and its archives" if catalog else ""))
            ds.failures.append(problem)
            notes.append(problem + " Other files the study offers can be asked for.")
            ui.warn(problem)
        else:
            notes.append("Fetched the {} processed file(s) the SDRF names; nothing else was downloaded."
                         .format(len(wanted)))
    else:
        # No processed files named: take the study's non-raw files, smallest first.
        rest = [(n, u, s) for n, (u, s) in catalog.items()]
        keep = []
        for name, url, size in rest:
            if _is_raw(name) or ".raw." in name.lower():
                skipped_raw += 1
            else:
                keep.append((name, url, size))
        keep.sort(key=lambda x: x[2] if x[2] is not None else float("inf"))
        ui.info("the SDRF names no processed file; taking the study's {} non-raw file(s), smallest first"
                .format(len(keep)))
        for name, url, size in keep:
            try:
                fetch(name, url, size, quiet=len(keep) > 12)
            except net.NetError as err:
                notes.append("Could not download {} ({})".format(name, err))
        if not got and too_large:
            raise FilesTooLarge("{}'s files were all refused by the size limits: {}".format(acc, ", ".join(too_large)))
        notes.append("The SDRF names no per-sample data files; used the study's own files.")
    if too_large:
        notes.append("NOT downloaded because of size: " + ", ".join(too_large))
    if skipped_raw:
        notes.append("Skipped {} raw files.".format(skipped_raw))
    return notes, got
