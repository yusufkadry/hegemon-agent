"""Describe the downloaded files so the AI can write the dataset-specific code.

Shows what a lab member would look at in the notebook before editing it: the GEO
sample metadata, platform columns, the first and last lines of tables (count
files hide htseq summary rows at the bottom), spreadsheet sheets, HDF5 keys.
"""
from collections import Counter, OrderedDict, defaultdict
import bz2
import gzip
import os
import re
from pathlib import Path

MAX_TOTAL = 70000
MAX_FILE = 3000


def _open_text(path):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    if p.endswith(".bz2"):
        return bz2.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def _cut(text, n):
    text = str(text)
    return text if len(text) <= n else text[:n] + "...[+{} chars]".format(len(text) - n)


def _head(path, n=12, width=260):
    lines = []
    try:
        with _open_text(path) as fh:
            for line in fh:
                lines.append(_cut(line.rstrip("\r\n"), width))
                if len(lines) >= n:
                    break
    except (OSError, EOFError, UnicodeError) as err:
        lines.append("<unreadable: {}>".format(err))
    return lines


def _tail(path, n=6, width=200, limit=150e6):
    if os.path.getsize(str(path)) > limit:
        return None
    from collections import deque
    buf = deque(maxlen=n)
    total = 0
    try:
        with _open_text(path) as fh:
            for line in fh:
                buf.append(_cut(line.rstrip("\r\n"), width))
                total += 1
    except (OSError, EOFError, UnicodeError):
        return None
    return list(buf), total


def _is_texty(name):
    low = name.lower()
    for ext in (".gz", ".bz2"):
        if low.endswith(ext):
            low = low[:-len(ext)]
    return low.endswith((".txt", ".tsv", ".csv", ".tab", ".count", ".counts", ".soft", ".sdrf", ".idf",
                         ".mtx", ".gct", ".res", ".dat", ".xls")) or "." not in low.rsplit("/", 1)[-1]


def _describe_file(path, rel):
    size = path.stat().st_size
    out = ["--- {} ({:,} bytes)".format(rel, size)]
    low = path.name.lower()
    try:
        if low.endswith((".xlsx", ".xlsm")):
            import openpyxl
            wb = openpyxl.load_workbook(str(path), read_only=True, data_only=True)
            for ws in wb.worksheets[:8]:
                out.append("  sheet '{}':".format(ws.title))
                for i, row in enumerate(ws.iter_rows(values_only=True)):
                    out.append("    " + _cut("\t".join("" if v is None else str(v) for v in row), 240))
                    if i >= 5:
                        break
            wb.close()
        elif low.endswith((".h5", ".h5ad", ".hdf5", ".loom")):
            import h5py
            with h5py.File(str(path), "r") as h5:
                items = []
                h5.visititems(lambda k, v: items.append("  {} {}".format(
                    k, getattr(v, "shape", "") if hasattr(v, "shape") else "(group)")))
                out.extend(items[:60])
        elif low.endswith((".rds", ".rda", ".rdata")):
            out.append("  R data file (readable with the 'rdata' package if installed)")
        elif _is_texty(low):
            out.extend("  " + line for line in _head(path))
            tail = _tail(path)
            if tail:
                lines, total = tail
                out.append("  ... last lines (file has {:,} lines):".format(total))
                out.extend("  " + line for line in lines)
        else:
            with open(str(path), "rb") as fh:
                out.append("  binary, starts with {}".format(fh.read(16)))
    except Exception as err:  # evidence must never stop the build
        out.append("  <could not inspect: {}>".format(err))
    return _cut("\n".join(out), MAX_FILE)


def _soft_summary(soft_path, platform, gse=None):
    """gse: the series as already read (metadata, platform tables, first sample tables), so the big
    SOFT file is not read again here."""
    if gse is None or not _has_tables(gse, platform):
        from .geo_soft import read_geo_soft
        gse = read_geo_soft(soft_path, sample_tables=2)
    m = gse.metadata
    out = ["GEO FAMILY SOFT ({}) parsed with read_geo_soft():".format(Path(soft_path).name)]
    for key in ("title", "type", "platform_id", "summary", "overall_design", "relation"):
        for v in m.get(key, [])[:6]:
            out.append("  series {}: {}".format(key, _cut(v, 900)))
    samples = [n for n, g in gse.gsms.items() if (g.metadata.get("platform_id") or [""])[0] == platform]
    out.append("  samples on {}: {} (series total {})".format(platform, len(samples), len(gse.gsms)))
    gpl = gse.gpls.get(platform)
    if gpl is not None:
        out.append("PLATFORM {} table: {} rows; columns:".format(platform, len(gpl.table)))
        for col, desc in gpl.columns["description"].items():
            out.append("  {}: {}".format(col, _cut(desc, 120)))
        if len(gpl.table):
            out.append("  first rows:")
            out.append("  " + _cut(gpl.table.head(3).to_csv(sep="\t", index=False), 1500).replace("\n", "\n  "))
    keys = OrderedDict()
    for name in samples:
        for v in gse.gsms[name].metadata.get("characteristics_ch1", []):
            k = v.split(":", 1)[0].strip() if ":" in v else "(no key)"
            keys.setdefault(k, Counter())[v.split(":", 1)[-1].strip()] += 1
    if keys:
        out.append("characteristics_ch1 keys across these samples (count, example values):")
        for k, values in keys.items():
            out.append("  {} ({}): {}".format(k, sum(values.values()), _cut(", ".join(v for v, _ in values.most_common(4)), 200)))
    out.append("GSM -> title (first 15):")
    for name in samples[:15]:
        out.append("  {} -> {}".format(name, _cut((gse.gsms[name].metadata.get("title") or [""])[0], 120)))
    for name in samples[:2]:
        g = gse.gsms[name]
        out.append("SAMPLE {} metadata:".format(name))
        for k, vals in g.metadata.items():
            if "protocol" in k or k in ("contact_name", "contact_email", "contact_institute"):
                continue
            out.append("  {}: {}".format(k, _cut(" | ".join(vals), 300)))
        if len(g.table):
            out.append("  table: {} rows, columns {}".format(len(g.table), list(g.table.columns)))
            for col, desc in g.columns["description"].items():
                out.append("    #{}: {}".format(col, _cut(desc, 200)))
            out.append("  " + g.table.head(3).to_csv(sep="\t", index=False).replace("\n", "\n  "))
        else:
            out.append("  (no data table in the SOFT file for this sample)")
    other = Counter((g.metadata.get("platform_id") or ["?"])[0] for g in gse.gsms.values())
    if len(other) > 1:
        out.append("samples per platform in this series: {}".format(dict(other)))
    return gse, samples, "\n".join(out)


def _has_tables(gse, platform):
    """False when the series was read without the tables shown here (then it is read again)."""
    entries = [gse.gpls.get(platform)] + [
        g for g in gse.gsms.values() if (g.metadata.get("platform_id") or [""])[0] == platform][:2]
    for entry in entries:
        if entry is None:
            continue
        count = str((entry.metadata.get("data_row_count") or ["0"])[0]).strip()
        if count.isdigit() and int(count) > 0 and len(entry.table) == 0:
            return False
    return True


def _size(n):
    if not n:
        return "size unknown"
    n = float(n)
    for unit in ("bytes", "KB", "MB", "GB"):
        if n < 1000 or unit == "GB":
            return "{:.0f} {}".format(n, unit) if unit == "bytes" else "{:.1f} {}".format(n, unit)
        n /= 1000.0


def catalog_text(ds, source_dir, max_lines=120):
    """What the study offers that is not downloaded, and what failed: the AI can ask for any of it."""
    source_dir = Path(source_dir)
    catalog = getattr(ds, "catalog", None) or {}
    waiting = [e for name, e in catalog.items() if not (source_dir / name).is_file()]
    out = []
    if waiting:
        out.append("\nFILES THIS STUDY OFFERS THAT ARE NOT DOWNLOADED ({} files). To get any of them, reply with "
                   "lines 'FETCH: <name, name pattern with *, or URL>' and no code:".format(len(waiting)))
        groups = OrderedDict()
        for e in waiting:
            groups.setdefault((re.sub(r"\d+", "#", e.name), e.note.split(" of ")[0]), []).append(e)
        shown = 0
        for (pattern, note), members in groups.items():
            if len(out) > max_lines:
                out.append("  ... and {} more files".format(len(waiting) - shown))
                break
            if len(members) > 4:
                total = sum(e.size or 0 for e in members)
                out.append("  {} files like {} [{}]{} e.g. {}".format(
                    len(members), pattern, note, " ({} in total)".format(_size(total)) if total else "",
                    ", ".join(e.name for e in members[:3])))
            else:
                for e in members:
                    out.append("  {} ({}) [{}]".format(e.name, _size(e.size), _cut(e.note, 120)))
            shown += len(members)
    failures = getattr(ds, "failures", None) or []
    if failures:
        out.append("\nDOWNLOAD PROBLEMS (these files did not arrive; ask for another file, or retry one if it may "
                   "have been temporary):")
        out.extend("  " + _cut(f, 300) for f in failures[:20])
    return "\n".join(out)


def build(ds, source_dir, platform=None, soft_path=None, notes=(), gse=None):
    source_dir = Path(source_dir)
    parts = ["DATASET {} ({})".format(ds.accession, ds.repository),
             "title: " + (ds.title or ""),
             "type: " + (ds.kind or ""),
             "organism: " + (ds.organism or "")]
    if ds.summary:
        parts.append("summary: " + _cut(ds.summary, 1500))
    for note in notes:
        parts.append("download note: " + note)
    samples = None
    if soft_path:
        gse, samples, text = _soft_summary(soft_path, platform, gse)
        parts.append(text)

    files = sorted(p for p in source_dir.rglob("*") if p.is_file() and not p.name.endswith(".part"))
    parts.append("\nFILES under SOURCE_DIR ({} files):".format(len(files)))
    groups = defaultdict(list)
    for p in files:
        rel = str(p.relative_to(source_dir))
        pattern = re.sub(r"\d+", "#", rel)
        groups[(str(p.parent), pattern)].append((p, rel))
    for (_, pattern), members in sorted(groups.items(), key=lambda kv: kv[1][0][1]):
        if len(members) > 4:
            parts.append("{} files like {} (e.g. {})".format(len(members), pattern,
                         ", ".join(rel for _, rel in members[:25]) + (" ..." if len(members) > 25 else "")))
        else:
            for p, rel in members:
                parts.append("{} ({:,} bytes)".format(rel, p.stat().st_size))
    parts.append(catalog_text(ds, source_dir))
    parts.append("\nFILE CONTENTS:")
    for (_, pattern), members in sorted(groups.items(), key=lambda kv: kv[1][0][1]):
        for p, rel in members[:2]:
            if p.name.endswith("_family.soft.gz"):
                continue
            parts.append(_describe_file(p, rel))
        if sum(len(x) for x in parts) > MAX_TOTAL:
            parts.append("[evidence truncated at {} characters]".format(MAX_TOTAL))
            break
    return gse, samples, "\n".join(parts)
