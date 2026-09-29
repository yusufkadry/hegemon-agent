"""Publish a finished build to your Hegemon site and prove it works through the real PHP pages.

Adds one [section] to explore.conf (the same keys hegemon.php reads) and one
explore.php?key= link to index.html, then checks, over HTTP:
  1. explore.php?key=KEY lists the dataset with the right sample count
     (Hegemon silently leaves out datasets whose expr file it cannot open)
  2. go=getdatajson returns a gene's row through idx pointer + expr file
  3. go=getthrjson returns that probe's thresholds from the thr file
On any failure both files are restored byte-for-byte.
"""
import fcntl
import json
import os
import re
import time
import urllib.parse
from pathlib import Path

from . import net, ui

PREFERRED_GENES = ("GAPDH", "ACTB", "Gapdh", "Actb", "B2M", "RPLP0")


class PublishError(Exception):
    pass


def _sections(text):
    sections, current = {}, None
    for line in text.splitlines():
        m = re.match(r"^\s*\[(.+)\]\s*$", line)
        if m:
            current = m.group(1)
            sections[current] = {}
        elif current is not None and "=" in line:
            k, v = line.split("=", 1)
            sections[current][k.strip()] = v.strip()
    return sections


def _atomic_write(path, data):
    path = Path(path)
    tmp = path.with_name(".{}.agent-tmp".format(path.name))
    tmp.write_bytes(data)
    try:
        os.chmod(str(tmp), os.stat(str(path)).st_mode & 0o7777)
    except OSError:
        pass
    os.replace(str(tmp), str(path))


def _unique(base, used):
    name, n = base, 2
    while name in used:
        name = "{}-{}".format(base, n)
        n += 1
    return name


def _json_from(body):
    body = body.strip()
    try:
        return json.loads(body)
    except ValueError:
        start, end = body.find("["), body.rfind("]")
        if start >= 0 and end > start:
            return json.loads(body[start:end + 1])
        raise PublishError("Hegemon did not return JSON: " + body[:300])


def pick_gene(run_dir, prefix, idx_path=None, expr_path=None):
    """A gene with at least two numeric values, preferring housekeeping genes."""
    best = None
    idx_path = idx_path or Path(run_dir) / "{}-idx.txt".format(prefix)
    expr_path = expr_path or Path(run_dir) / "{}-expr.txt".format(prefix)
    with open(str(idx_path), "rb") as fi, open(str(expr_path), "rb") as fe:
        fi.readline()
        for line in fi:
            probe, ptr, name, _ = (line.decode("utf-8", "replace").rstrip("\r\n").split("\t") + ["", "", "", ""])[:4]
            gene = name.split(" /// ")[0].strip()
            if not gene or gene == "---" or re.search(r"\s", gene):
                continue
            fe.seek(int(ptr))
            vals = fe.readline().decode("utf-8", "replace").rstrip("\r\n").split("\t")[2:]
            if sum(1 for v in vals if v.strip()) < 2:
                continue
            if gene in PREFERRED_GENES:
                return gene
            best = best or gene
    if not best:
        raise PublishError("No gene with values found in idx for the live check")
    return best


def live_check(cfg, section, key, run_dir, prefix, n_samples):
    base = cfg.site_url + "/explore.php"
    try:
        page = net.get_text(base + "?key=" + urllib.parse.quote(key), timeout=60, retries=1)
    except net.NetError as err:
        raise PublishError("explore.php crashed or was unreachable with the new dataset added ({}). Usually PHP "
                           "could not open one of the dataset files: check that the web server can read the build "
                           "folder and every folder above it.".format(err))
    option = re.search(r'<option value="{}"[^>]*>(.*?)</option>'.format(re.escape(section)), page, re.S)
    if not option:
        raise PublishError("explore.php does not list the dataset. Hegemon leaves out datasets whose expr file the "
                           "web server cannot open: check that the build folder and every folder above it are "
                           "readable by the web server.")
    shown = re.search(r"\(n = (-?\d+)\)", option.group(1))
    if not shown or int(shown.group(1)) != n_samples:
        raise PublishError("explore.php shows n = {} but the dataset has {} samples".format(
            shown.group(1) if shown else "?", n_samples))
    ui.ok("explore.php lists {} with n = {}".format(section, n_samples))
    if os.environ.get("HEGEMON_LIVE_CHECK", "full") == "list":
        ui.warn("HEGEMON_LIVE_CHECK=list: skipping the gene and threshold lookups")
        return

    gene = pick_gene(run_dir, prefix)
    q = "&A={0}&B={0}&id={1}&key={2}".format(urllib.parse.quote(gene), urllib.parse.quote(section),
                                              urllib.parse.quote(key))
    try:
        rows = _json_from(net.get_text(base + "?go=getdatajson" + q, timeout=120, retries=1))
    except net.NetError as err:
        raise PublishError("Hegemon's getdatajson (the endpoint the lab's HegemonUtil.py uses) failed for {}: {}. "
                           "If that endpoint is broken for every dataset on this server, set "
                           "HEGEMON_LIVE_CHECK=list.".format(gene, err))
    if len(rows) < 2 or len(rows[0]) != n_samples + 2 or any(len(r) != n_samples + 2 for r in rows[1:]):
        raise PublishError("getdatajson for {} returned {} rows of the wrong shape".format(gene, len(rows)))
    probes = {r[0] for r in rows[1:]}
    ui.ok("Hegemon returned {} ({} probe row(s) x {} samples)".format(gene, len(rows) - 1, n_samples))

    try:
        thr = _json_from(net.get_text(base + "?go=getthrjson" + q, timeout=120, retries=1))
    except net.NetError as err:
        raise PublishError("Hegemon's getthrjson failed for {}: {}".format(gene, err))
    if not any(isinstance(r, list) and r and r[0] in probes for r in thr):
        raise PublishError("getthrjson did not read {}'s thresholds from the thr file".format(gene))
    ui.ok("Hegemon read {}'s thresholds from the thr file".format(gene))
    check_metadata_live(cfg, section, key, n_samples)


def check_metadata_live(cfg, section, key, n_samples):
    """Prove Hegemon can read the sample metadata AND group samples by it.

    Two failures look like "the dataset opens but nothing can be plotted":
    - survival's ArrayIds do not match the expr columns, so Hegemon returns every value blank;
    - the only metadata is titles/IDs (all different) or organism (all the same), which group nothing.
    Both are refused here, through the same endpoints the lab's HegemonUtil.py uses.
    """
    from .checks import comparison_field, groups_text, why_not_grouping
    base = cfg.site_url + "/explore.php"
    q = "&id={}&key={}".format(urllib.parse.quote(section), urllib.parse.quote(key))
    try:
        headers = _json_from(net.get_text(base + "?go=getpatientinfojson" + q, timeout=60, retries=1))
    except net.NetError as err:
        raise PublishError("Hegemon's getpatientinfojson failed: {}. The survival file is not readable by the web "
                           "server.".format(err))
    if not isinstance(headers, list) or len(headers) < 4:
        raise PublishError("Hegemon lists no sample metadata fields for this dataset (it returned {}). Without them "
                           "the dataset opens but nothing can be grouped or plotted.".format(headers))
    fields = [(i, h) for i, h in enumerate(headers) if str(h).startswith(("c ", "n "))]
    if not fields:
        raise PublishError("Hegemon lists survival columns {} but none is a 'c ' or 'n ' metadata field, so the "
                           "dataset cannot be grouped or plotted.".format(headers[1:]))
    groupable, rejected, any_value = [], [], False
    for index, name in fields[:80]:
        try:
            payload = _json_from(net.get_text(
                base + "?go=getpatientdatajson&clinical={}".format(index) + q, timeout=120, retries=1))
        except net.NetError as err:
            raise PublishError("Hegemon's getpatientdatajson failed for '{}': {}".format(name, err))
        if not (isinstance(payload, list) and len(payload) == 2):
            raise PublishError("getpatientdatajson returned an unexpected shape for '{}'".format(name))
        expr_header, values = payload
        cells = values[2:]
        if len(expr_header) - 2 != n_samples or len(cells) != n_samples:
            raise PublishError("Hegemon matched {} of {} samples for '{}'".format(len(cells), n_samples, name))
        if any(str(v).strip() for v in cells):
            any_value = True
        counts = comparison_field(name, cells, n_samples)
        if counts:
            groupable.append((name, counts))
            if len(groupable) >= 3:
                break
        else:
            rejected.append("{} ({})".format(name, why_not_grouping(name, cells, n_samples)))
    if not any_value:
        raise PublishError(
            "Hegemon read the metadata fields {} but every sample's value came back blank. That is the failure where "
            "the dataset opens and no plot can be made: the ArrayIds in survival.txt do not match the sample columns "
            "in expr.txt as Hegemon joins them.".format([n for _, n in fields][:5]))
    if not groupable:
        raise PublishError(
            "Hegemon reads the metadata, but no field splits the samples into groups it can compare, so the page "
            "would load with nothing to plot by. Fields: {}".format("; ".join(rejected[:8])))
    for name, counts in groupable:
        ui.ok("Hegemon can group samples by '{}': {}".format(name, groups_text(counts)))


# ------------------------------------------------------------------ the page's opening genes
#
# Hegemon's groups.php draws the plot and the "Select Patient Information" dropdown only
# when BOTH Gene A and Gene B are found in the selected dataset; otherwise the page shows
# nothing. It defaults to TYROBP / FCER1G, which many datasets lack. So every published
# link names two genes that exist in every dataset under its key, with cmd=explore so the
# page draws the plot on load.

HOUSEKEEPING = ("ACTB", "GAPDH", "B2M", "RPLP0", "PPIA", "TUBB", "HPRT1", "TBP", "PGK1", "YWHAZ")


def idx_names(idx_path):
    """Names as Hegemon's readIndexFile looks them up: Name split on ' /// ', upper-cased."""
    names = set()
    with open(str(idx_path), "rb") as fh:
        fh.readline()
        for line in fh:
            parts = line.decode("utf-8", "replace").rstrip("\r\n").split("\t")
            if len(parts) != 4:
                continue
            for token in parts[2].split(" /// "):
                token = token.strip().upper()
                if token and token != "---":
                    names.add(token)
    return names


def choose_genes(idx_paths, requested=None):
    """(A, B, how) present in every dataset, or None."""
    sets = []
    for path in idx_paths:
        try:
            sets.append(idx_names(path))
        except OSError:
            return None
    if not sets:
        return None
    common = set.intersection(*sets)
    wanted = [g.upper() for g in (requested or []) if g]
    picked = [g for g in wanted if g in common][:2]
    how = "for the topic" if picked else ""
    for g in HOUSEKEEPING:
        if len(picked) == 2:
            break
        if g in common and g not in picked:
            picked.append(g)
            how = how + (" + " if how else "") + "housekeeping"
    if len(picked) < 2:
        plain = sorted(g for g in common if re.match(r"^[A-Z][A-Z0-9-]{1,14}$", g) and not g.startswith(("LOC", "ENS")))
        for g in plain:
            if len(picked) == 2:
                break
            if g not in picked:
                picked.append(g)
                how = how or "shared by every dataset"
    return (picked[0], picked[1], how) if len(picked) == 2 else None


def key_link(key, pair):
    return "explore.php?key={}".format(key) + ("&A={}&B={}&cmd=explore".format(pair[0], pair[1]) if pair else "")


def _set_homepage_href(text, key, href):
    pattern = re.compile(r'(href=["\'])explore\.php\?key={}(?:&[^"\']*)?(["\'])'.format(re.escape(key)))
    return pattern.sub(lambda m: m.group(1) + href.replace("&", "&amp;") + m.group(2), text)


def check_page_opens(cfg, section, key, pair):
    """The page must draw the plot and the metadata dropdown with these genes, as a visitor sees it."""
    q = "?go=explore&A={}&B={}&id={}&key={}".format(urllib.parse.quote(pair[0]), urllib.parse.quote(pair[1]),
                                                     urllib.parse.quote(section), urllib.parse.quote(key))
    try:
        html = net.get_text(cfg.site_url + "/explore.php" + q, timeout=120, retries=1)
    except net.NetError as err:
        raise PublishError("Hegemon could not draw {} vs {}: {}".format(pair[0], pair[1], err))
    if 'id="clinical0"' not in html or "groupPlot" not in html:
        raise PublishError("Hegemon did not draw the plot and metadata dropdown for {} vs {}".format(*pair))
    options = re.findall(r'<option value="\d+">([^<]+)</option>', html)
    meta = [o for o in options if o.startswith(("c ", "n "))]
    ui.ok("the page opens on {} vs {} with the plot and the Select Patient Information dropdown ({} fields)"
          .format(pair[0], pair[1], len(meta)))



def register(cfg, receipt, label=None, check=True, key=None, link_label=None, genes=None):
    site = cfg.site_root
    conf, index = site / "explore.conf", site / "index.html"
    if not conf.is_file() or not index.is_file():
        raise PublishError("No explore.conf/index.html in {} (set HEGEMON_SITE_ROOT)".format(site))
    run_dir, prefix = Path(receipt["run_dir"]), receipt["prefix"]
    files = {k: str(run_dir / "{}-{}.txt".format(prefix, k)) for k in ("expr", "idx", "survival", "ih", "info")}
    label = re.sub(r"[<>\r\n\t]", " ", label or "{} {}".format(receipt["accession"], receipt.get("title", "")))[:160].strip()

    with open(str(site / ".hegemon-agent.lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        old_conf, old_index = conf.read_bytes(), index.read_bytes()
        conf_text, index_text = old_conf.decode("utf-8", "replace"), old_index.decode("utf-8", "replace")
        sections = _sections(conf_text)
        used_keys = {v.get("key", "") for v in sections.values()} | set(re.findall(r"explore\.php\?key=([^\"'&<\s]+)", index_text))
        section = _unique(prefix, set(sections))
        if key:
            # A shared key: several datasets carry it, and Hegemon's explore.php?key=<key>
            # lists every one of them in its dropdown.
            key = str(key).strip()
            if not re.match(r"^[A-Za-z0-9._+-]+$", key):
                raise PublishError("A key may only contain letters, digits, '.', '_', '+' or '-' (Hegemon splits the "
                                   "key field on ':'); got {!r}".format(key))
            for name, values in sections.items():
                if values.get("source", "").strip() == receipt["accession"] and key in values.get("key", "").split(":"):
                    raise PublishError("{} is already published under key '{}' as section [{}].".format(
                        receipt["accession"], key, name))
        else:
            key = _unique(re.sub(r"[^a-z0-9]", "", prefix.lower()) or "dataset", used_keys)
        block = ("\n[{}]\nname= {}\nexpr= {}\nindex= {}\nsurvival= {}\nindexHeader= {}\ninfo= {}\nkey= {}\nsource= {}\n"
                 .format(section, label, files["expr"], files["idx"], files["survival"], files["ih"], files["info"],
                         key, receipt["accession"]))
        pos = index_text.lower().rfind("</table>")
        if pos < 0:
            raise PublishError("index.html has no </table>; not guessing where the link goes")
        if re.search(r'explore\.php\?key={}(?=["\'&<\s])'.format(re.escape(key)), index_text):
            new_index = index_text  # this key already has a homepage link; the dataset joins it
        else:
            shown = re.sub(r"[<>\r\n\t]", " ", str(link_label or label))[:160].strip()
            row = '<tr>\n<td>\n<a href="explore.php?key={}"> {} </a>\n</td>\n</tr>\n'.format(key, shown)
            new_index = index_text[:pos] + row + index_text[pos:]

        backup = site / ".hegemon-agent-backups" / time.strftime("%Y%m%d-%H%M%S")
        backup.mkdir(parents=True, exist_ok=True)
        (backup / "explore.conf").write_bytes(old_conf)
        (backup / "index.html").write_bytes(old_index)
        try:
            _atomic_write(conf, old_conf.rstrip(b"\n") + b"\n" + block.encode("utf-8"))
            _atomic_write(index, new_index.encode("utf-8"))
            if check:
                live_check(cfg, section, key, run_dir, prefix, int(receipt["samples"]))
            # Two genes every dataset under this key has, so the page opens on a plot.
            idx_paths = [v.get("index") for v in _sections(conf.read_text(errors="replace")).values()
                         if key in [k.strip() for k in v.get("key", "").split(":")] and v.get("index")]
            pair = choose_genes(idx_paths, genes)
            if pair:
                if check:
                    check_page_opens(cfg, section, key, pair)
                _atomic_write(index, _set_homepage_href(index.read_text(errors="replace"), key,
                                                        key_link(key, pair)).encode("utf-8"))
                ui.info("genes chosen {}: {} vs {}".format(pair[2], pair[0], pair[1]))
            else:
                ui.warn("no two genes are shared by every dataset under '{}'; the link opens without a plot".format(key))
        except Exception as err:
            _atomic_write(conf, old_conf)
            _atomic_write(index, old_index)
            restored = conf.read_bytes() == old_conf and index.read_bytes() == old_index
            raise PublishError("{}\n   explore.conf and index.html {}.".format(
                err, "were restored exactly" if restored else "COULD NOT BE RESTORED - copies are in " + str(backup)))
    link = "{}/{}".format(cfg.site_url, key_link(key, pair))
    return section, key, link + ("&id={}".format(urllib.parse.quote(section)) if pair else "")


def diagnose(cfg, only_key=None, min_samples=0):
    """Check the datasets already registered on this Hegemon site, one row each.

    This answers "the page opens but I cannot plot anything" without rebuilding:
    it uses the same live checks as publishing, against what is already in explore.conf.
    """
    conf = cfg.site_root / "explore.conf"
    if not conf.is_file():
        raise PublishError("No explore.conf in {} (set HEGEMON_SITE_ROOT)".format(cfg.site_root))
    rows = []
    for section, values in _sections(conf.read_text(errors="replace")).items():
        keys = [k.strip() for k in values.get("key", "").split(":") if k.strip()]
        if only_key and only_key not in keys:
            continue
        row = {"section": section, "source": values.get("source", ""), "name": values.get("name", ""),
               "keys": keys, "problems": [], "notes": []}
        rows.append(row)
        expr = Path(values.get("expr", ""))
        idx = Path(values.get("index", ""))
        surv = Path(values.get("survival", ""))
        for label, path in (("expr", expr), ("index", idx), ("survival", surv),
                            ("indexHeader", Path(values.get("indexHeader", "")))):
            if not str(path) or not path.is_file():
                row["problems"].append("{} file is missing: {}".format(label, path))
            elif path.stat().st_size == 0:
                row["problems"].append("{} file is empty".format(label))
        if row["problems"]:
            continue
        try:
            with open(str(expr), "r", encoding="utf-8", errors="replace") as fh:
                header = fh.readline().rstrip("\r\n").split("\t")
            n_samples = len(header) - 2
            row["samples"] = n_samples
            with open(str(surv), "r", encoding="utf-8", errors="replace") as fh:
                surv_header = fh.readline().rstrip("\r\n").split("\t")
            meta = [c for c in surv_header[3:] if c.startswith(("c ", "n "))]
            row["metadata_columns"] = len(meta)
            if not meta:
                row["problems"].append("survival.txt has no 'c ' or 'n ' metadata column, so nothing can be "
                                       "grouped or plotted")
            ids = set()
            with open(str(surv), "r", encoding="utf-8", errors="replace") as fh:
                fh.readline()
                for line in fh:
                    first = line.split("\t", 1)[0].strip()
                    if first:
                        ids.add(first)
            matched = sum(1 for s in header[2:] if s in ids)
            row["matched_samples"] = matched
            if matched == 0:
                row["problems"].append("none of the {} sample columns in expr.txt appear in survival.txt's ArrayId "
                                       "column, so every metadata value is blank in Hegemon".format(n_samples))
            elif matched < n_samples:
                row["notes"].append("{} of {} samples have metadata".format(matched, n_samples))
            if matched and meta:
                import csv
                import pandas as pd
                from .checks import comparison_field, groups_text
                table = pd.read_csv(str(surv), sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)
                groups = [(c, comparison_field(c, table[c].tolist(), len(table))) for c in meta if c in table.columns]
                groups = [(c, g) for c, g in groups if g]
                row["groups"] = [(c, groups_text(g)) for c, g in groups[:3]]
                if not groups:
                    row["problems"].append("no metadata column splits the samples into groups (only IDs, titles or "
                                           "constant fields), so there is nothing to plot by")
            if min_samples and n_samples < min_samples:
                row["problems"].append("only {} samples (fewer than {}), too few points for Hegemon's "
                                       "comparisons".format(n_samples, min_samples))
        except OSError as err:
            row["problems"].append("could not read the files: {}".format(err))
            continue
        if not keys:
            row["notes"].append("no key, so it is listed for every visitor")
        key = keys[0] if keys else ""
        try:
            page = net.get_text("{}/explore.php?key={}".format(cfg.site_url, urllib.parse.quote(key)),
                                timeout=60, retries=1)
            if not re.search(r'<option value="{}"'.format(re.escape(section)), page):
                row["problems"].append("Hegemon does not list it at explore.php?key={} (the web server probably "
                                       "cannot read the data files)".format(key))
            else:
                gene = pick_gene(None, None, idx, expr)
                row["gene"] = gene
                q = "&A={0}&B={0}&id={1}&key={2}".format(urllib.parse.quote(gene),
                                                         urllib.parse.quote(section), urllib.parse.quote(key))
                data = _json_from(net.get_text(cfg.site_url + "/explore.php?go=getdatajson" + q, timeout=120,
                                               retries=1))
                if len(data) < 2:
                    row["problems"].append("gene lookup returned nothing for " + gene)
                if not any("splits the samples" in pr or "ArrayId" in pr for pr in row["problems"]):
                    try:
                        check_metadata_live(cfg, section, key, row.get("samples", 0))
                    except PublishError as err:
                        row["problems"].append(str(err).split(". ")[0])
        except (net.NetError, PublishError, OSError, ValueError) as err:
            row["problems"].append("live check failed: {}".format(err))
    return rows


# ------------------------------------------------------------------ unpublish

def _section_spans(text):
    heads = list(re.finditer(r"(?m)^[ \t]*\[([^\]\r\n]+)\][ \t]*\r?$", text))
    spans = []
    for i, m in enumerate(heads):
        end = heads[i + 1].start() if i + 1 < len(heads) else len(text)
        spans.append((m.group(1).strip(), m.start(), end))
    return spans


def _remove_homepage_link(text, key):
    """Remove the <tr> rows whose link is explore.php?key=<key>. Returns (text, rows_removed)."""
    removed = 0
    pattern = re.compile(r'explore\.php\?key={}(?=["\'&<\s])'.format(re.escape(key)))
    for m in reversed(list(pattern.finditer(text))):
        start = text.lower().rfind("<tr", 0, m.start())
        end = text.lower().find("</tr>", m.end())
        if start < 0 or end < 0 or text.lower().find("</tr>", start, m.start()) >= 0:
            continue  # not inside a table row we can identify; leave it rather than guess
        end += len("</tr>")
        if text[end:end + 1] == "\n":
            end += 1
        text = text[:start] + text[end:]
        removed += 1
    return text, removed


def plan_unpublish(cfg, sections=None, key=None):
    conf = cfg.site_root / "explore.conf"
    if not conf.is_file():
        raise PublishError("No explore.conf in {} (set HEGEMON_SITE_ROOT)".format(cfg.site_root))
    values = _sections(conf.read_text(errors="replace"))
    unknown = [s for s in (sections or []) if s not in values]
    if unknown:
        raise PublishError("Not in explore.conf: {}".format(", ".join(unknown)))
    remove, strip = [], []
    for name, v in values.items():
        keys = [k.strip() for k in v.get("key", "").split(":") if k.strip()]
        if sections and name in sections:
            remove.append(name)
        elif key and key in keys:
            (remove if len(keys) == 1 else strip).append(name)
    return remove, strip, values


def unpublish(cfg, sections=None, key=None, check=True):
    """Take datasets off the site: their explore.conf sections and any homepage link left pointing
    at nothing. The data files stay where they are. Both site files are backed up first and
    restored byte-for-byte if the page stops loading."""
    site = cfg.site_root
    conf, index = site / "explore.conf", site / "index.html"
    with open(str(site / ".hegemon-agent.lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        remove, strip, values = plan_unpublish(cfg, sections, key)
        if not remove and not strip:
            raise PublishError("Nothing in explore.conf matches" + (" key '{}'".format(key) if key else ""))
        old_conf, old_index = conf.read_bytes(), index.read_bytes() if index.is_file() else b""
        text = old_conf.decode("utf-8", "replace")
        out, pos = [], 0
        for name, start, end in _section_spans(text):
            out.append(text[pos:start])
            pos = end
            block = text[start:end]
            if name in remove:
                continue
            if name in strip:
                block = re.sub(r"(?m)^([ \t]*key[ \t]*=[ \t]*)(.*)$", lambda m: m.group(1) + ":".join(
                    k for k in m.group(2).split(":") if k.strip() and k.strip() != key), block, count=1)
            out.append(block)
        out.append(text[pos:])
        new_conf = re.sub(r"\n{3,}", "\n\n", "".join(out))

        still_used = set()
        for v in _sections(new_conf).values():
            still_used |= {k.strip() for k in v.get("key", "").split(":") if k.strip()}
        affected = set([key] if key else [])
        for name in remove + strip:
            affected |= {k.strip() for k in values[name].get("key", "").split(":") if k.strip()}
        dead = sorted(k for k in affected if k not in still_used)
        new_index, links = old_index.decode("utf-8", "replace"), 0
        for k in dead:
            new_index, n = _remove_homepage_link(new_index, k)
            links += n

        backup = site / ".hegemon-agent-backups" / ("unpublish-" + time.strftime("%Y%m%d-%H%M%S"))
        backup.mkdir(parents=True, exist_ok=True)
        (backup / "explore.conf").write_bytes(old_conf)
        (backup / "index.html").write_bytes(old_index)
        try:
            _atomic_write(conf, new_conf.encode("utf-8"))
            if index.is_file():
                _atomic_write(index, new_index.encode("utf-8"))
            if check:
                page = net.get_text(cfg.site_url + "/explore.php", timeout=60, retries=1)
                low = page.casefold()
                if "select dataset" not in low or "fatal error" in low:
                    raise PublishError("explore.php no longer renders after the change")
        except Exception as err:
            _atomic_write(conf, old_conf)
            if index.is_file():
                _atomic_write(index, old_index)
            raise PublishError("{} - explore.conf and index.html were restored from {}".format(err, backup))
    return remove, strip, dead, links, backup


def relink(cfg, only_key=None, genes=None, check=True):
    """Point each key's homepage link at two genes every dataset under it has (with cmd=explore),
    after proving the page draws the plot and the metadata dropdown with them. Returns rows."""
    site = cfg.site_root
    conf, index = site / "explore.conf", site / "index.html"
    if not conf.is_file() or not index.is_file():
        raise PublishError("No explore.conf/index.html in {} (set HEGEMON_SITE_ROOT)".format(site))
    rows = []
    with open(str(site / ".hegemon-agent.lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        old_index = index.read_bytes()
        text = old_index.decode("utf-8", "replace")
        by_key = {}
        for section, v in _sections(conf.read_text(errors="replace")).items():
            for k in [k.strip() for k in v.get("key", "").split(":") if k.strip()]:
                by_key.setdefault(k, []).append((section, v.get("index")))
        for key, members in sorted(by_key.items()):
            if only_key and key != only_key:
                continue
            row = {"key": key, "datasets": len(members), "pair": None, "linked": False, "problem": ""}
            rows.append(row)
            pair = choose_genes([i for _s, i in members if i], genes)
            if not pair:
                row["problem"] = "no two genes are shared by every dataset under this key"
                continue
            row["pair"] = pair
            if check:
                try:
                    check_page_opens(cfg, members[0][0], key, pair)
                except PublishError as err:
                    row["problem"] = str(err)
                    continue
            new = _set_homepage_href(text, key, key_link(key, pair))
            row["linked"] = new != text or "key={}&".format(key) in text
            if new == text and not re.search(r'explore\.php\?key={}(?=["\'&<\s])'.format(re.escape(key)), text):
                row["problem"] = "no homepage link for this key (use its URL directly)"
            text = new
        if text.encode("utf-8") != old_index:
            backup = site / ".hegemon-agent-backups" / ("relink-" + time.strftime("%Y%m%d-%H%M%S"))
            backup.mkdir(parents=True, exist_ok=True)
            (backup / "index.html").write_bytes(old_index)
            _atomic_write(index, text.encode("utf-8"))
    return rows
