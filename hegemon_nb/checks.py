"""Checks every build must pass before the lab's idx code and script run.

These are the output rules from the notebooks plus the sample-identity rules:
exact one-to-one sample matching, nothing invented, nothing silently dropped.
Messages are written so they can be handed straight back to the AI.
"""
import csv
import re
from collections import Counter, OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

ENSEMBL = re.compile(r"^ENS[A-Z]*[GT]\d{6,}")


class Result(object):
    def __init__(self):
        self.errors, self.warnings, self.stats = [], [], {}

    @property
    def ok(self):
        return not self.errors

    def text(self):
        lines = ["ERROR: " + e for e in self.errors] + ["WARNING: " + w for w in self.warnings]
        return "\n".join(lines)


def header_of(path):
    with open(str(path), "rb") as fh:
        return fh.readline().decode("utf-8", "replace").rstrip("\r\n").split("\t")


def _field_counts(path, expected, label, res):
    bad, first = 0, None
    with open(str(path), "rb") as fh:
        fh.readline()
        for n, line in enumerate(fh, start=2):
            count = line.count(b"\t") + 1
            if count != expected:
                bad += 1
                first = first or (n, count)
    if bad:
        res.errors.append("{}: {} rows have the wrong number of tab-separated fields (line {} has {}, "
                          "header has {}). A value contains a tab/newline or rows are misaligned."
                          .format(label, bad, first[0], first[1], expected))


def _table(path):
    return pd.read_csv(str(path), sep="\t", dtype=str, keep_default_na=False, na_filter=False,
                       quoting=csv.QUOTE_NONE)


_BOOLEAN_CELLS = (b"\tTrue", b"\ttrue", b"\tTRUE", b"\tFalse", b"\tfalse", b"\tFALSE")


def _expr_numbers(path, ncols):
    """The whole expr file with its values read by pandas' C number parser (fast), or None.

    None means some value is not a plain number, so the slower cell-by-cell check below runs
    and names it. That parser rejects everything pd.to_numeric rejects except true/false
    spellings, so a file containing one of those always takes the slow check.
    """
    import mmap
    with open(str(path), "rb") as fh:
        try:
            view = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        except (ValueError, OSError):
            return None
        try:
            if any(view.find(cell) >= 0 for cell in _BOOLEAN_CELLS):
                return None
        finally:
            view.close()
    kinds = dict([(0, str), (1, str)] + [(i, "float64") for i in range(2, ncols)])
    try:
        return pd.read_csv(str(path), sep="\t", header=0, names=list(range(ncols)), dtype=kinds,
                           keep_default_na=False, na_values={i: [""] for i in range(2, ncols)},
                           quoting=csv.QUOTE_NONE)
    except (ValueError, TypeError, pd.errors.ParserError):
        return None


def _examples(values, n=5):
    values = list(values)
    return ", ".join(values[:n]) + (" ..." if len(values) > n else "")


def check_core(run_dir, prefix, repository, series_samples=None, single_cell=False, genome_ok=False):
    res = Result()
    run_dir = Path(run_dir)
    paths = {k: run_dir / "{}-{}.txt".format(prefix, k) for k in ("expr", "survival", "ih")}
    missing = [p.name for p in paths.values() if not p.is_file() or p.stat().st_size == 0]
    if missing:
        res.errors.append("Missing or empty output files: " + ", ".join(missing))
        return res

    # ---- expr
    header = header_of(paths["expr"])
    samples = header[2:]
    if header[:2] != ["ProbeID", "Name"]:
        res.errors.append("expr header must start with 'ProbeID\\tName'; it starts with {}".format(header[:3]))
    if len(samples) < 2:
        res.errors.append("expr has {} sample columns; need at least 2".format(len(samples)))
    dups = sorted({s for s in samples if samples.count(s) > 1}) if len(samples) < 5000 else \
        sorted(pd.Series(samples)[pd.Series(samples).duplicated()].unique())
    if dups:
        res.errors.append("Duplicate sample columns in expr: " + _examples(dups))
    if any(not s.strip() or s != s.strip() for s in samples):
        res.errors.append("Blank or space-padded sample column names in expr")
    before = len(res.errors)
    _field_counts(paths["expr"], len(header), "expr", res)

    ids, n_rows, n_vals, n_blank, n_bad, n_inf, n_neg = [], 0, 0, 0, 0, 0, 0
    names_blank = names_eq_id = ens_ids = 0
    vmax, vmin, integers = -np.inf, np.inf, True
    bad_examples = []
    fast = _expr_numbers(paths["expr"], len(header)) if len(res.errors) == before and len(header) > 2 else None
    reader = [fast] if fast is not None else pd.read_csv(
        str(paths["expr"]), sep="\t", dtype=str, keep_default_na=False, na_filter=False,
        quoting=csv.QUOTE_NONE, chunksize=20000)
    for chunk in reader:
        n_rows += len(chunk)
        col_id, col_name = chunk.iloc[:, 0], chunk.iloc[:, 1]
        ids.extend(col_id.tolist())
        names_blank += int((col_name.str.strip() == "").sum())
        names_eq_id += int((col_name == col_id).sum())
        ens_ids += int(col_id.str.match(ENSEMBL).sum())
        if fast is not None:   # every value parsed as a number; the only NaNs are blank cells
            nums = chunk.iloc[:, 2:].to_numpy(dtype="float64")
            blank = np.isnan(nums)
        else:
            raw = chunk.iloc[:, 2:]
            blank = (raw == "").to_numpy()
            nums = raw.apply(pd.to_numeric, errors="coerce").to_numpy(dtype="float64")
            bad = ~blank & np.isnan(nums)
            if bad.any():
                n_bad += int(bad.sum())
                if len(bad_examples) < 5:
                    bad_examples.extend(np.unique(raw.to_numpy()[bad])[:5].tolist())
        finite = np.isfinite(nums)
        n_inf += int(np.isinf(nums).sum())
        n_vals += int(finite.sum())
        n_blank += int(blank.sum())
        if finite.any():
            vals = nums[finite]
            vmax, vmin = max(vmax, vals.max()), min(vmin, vals.min())
            n_neg += int((vals < 0).sum())
            if integers and not np.all(vals == np.floor(vals)):
                integers = False
    res.stats.update(samples=len(samples), probes=n_rows, values=n_vals, blank_values=n_blank,
                     min=None if n_vals == 0 else round(float(vmin), 4),
                     max=None if n_vals == 0 else round(float(vmax), 4))
    id_series = pd.Series(ids)
    if (id_series.str.strip() == "").any():
        res.errors.append("{} rows have a blank ProbeID".format(int((id_series.str.strip() == "").sum())))
    dup_ids = id_series[id_series.duplicated()].unique()
    if len(dup_ids):
        res.errors.append("{} ProbeIDs repeat (e.g. {}). ProbeID must be unique; Name may repeat. Keep every row "
                          "and make the IDs unique (e.g. SYMBOL, SYMBOL_2) instead of dropping rows."
                          .format(len(dup_ids), _examples(dup_ids)))
    if n_bad:
        res.errors.append("{} expression values are not numbers (e.g. {}). Missing values must be blank."
                          .format(n_bad, _examples(bad_examples)))
    if n_inf:
        res.errors.append("{} expression values are infinite".format(n_inf))
    if n_vals == 0:
        res.errors.append("expr has no numeric values")
    elif integers and vmax > 100:
        res.errors.append("Values are whole numbers up to {:.0f}: these look like raw counts. The notebooks apply "
                          "CPM then log2(x+1) (cpm_log2) to RNA-seq.".format(vmax))
    elif vmax > 100:
        res.errors.append("Values go up to {:.4g}: not log scale. The notebooks take log2(x+1) of unlogged data "
                          "(and CPM first for RNA-seq).".format(vmax))
    if n_neg:
        res.warnings.append("{} negative values (fine for log-ratios/centered data; unexpected after log2(x+1))"
                            .format(n_neg))
    if n_rows and names_blank == n_rows:
        res.errors.append("Every Name is blank; Name must hold the gene symbol (Symbol:Definition for arrays)")
    elif n_rows and names_blank / float(n_rows) > 0.2:
        res.warnings.append("{:.0%} of rows have a blank Name (not searchable by gene)".format(names_blank / float(n_rows)))
    if n_rows and genome_ok and ens_ids / float(n_rows) > 0.5 and names_eq_id / float(n_rows) > 0.5:
        res.errors.append("Names are Ensembl IDs. Map them to gene names with lab_gene_names() like the notebook.")

    # ---- survival
    sh = header_of(paths["survival"])
    if len(sh) < 3 or sh[0].lower() != "arrayid" or sh[1] != "time" or sh[2] != "status":
        res.errors.append("survival header must start with ArrayId, time, status; got {}".format(sh[:4]))
    extra = [c for c in sh[3:] if not (c.startswith("c ") or c.startswith("n "))]
    if extra:
        res.errors.append("survival metadata columns must start with 'c ' (or 'n '): " + _examples(extra))
    _field_counts(paths["survival"], len(sh), "survival", res)
    surv = _table(paths["survival"])
    s_ids = surv.iloc[:, 0].tolist() if len(surv.columns) else []
    if len(set(s_ids)) != len(s_ids):
        res.errors.append("survival has repeated ArrayIds: " + _examples(sorted({x for x in s_ids if s_ids.count(x) > 1})))
    meta_cols = [c for c in surv.columns[3:]]
    filled = [c for c in meta_cols if (surv[c].str.strip() != "").any()]
    if not filled:
        res.errors.append("survival has no sample metadata values (every metadata column is empty or absent)")
    res.stats["metadata_columns"] = len(filled)
    for col in ("time", "status"):
        if col in surv.columns:
            vals = surv[col][surv[col].str.strip() != ""]
            if len(vals) and pd.to_numeric(vals, errors="coerce").isna().any():
                res.errors.append("survival '{}' must be numeric or blank".format(col))

    # ---- ih
    ih_header = header_of(paths["ih"])
    if ih_header != ["ArrayID", "ArrayHeader", "ClinicalHeader"]:
        res.errors.append("ih header must be ArrayID, ArrayHeader, ClinicalHeader; got {}".format(ih_header))
    _field_counts(paths["ih"], len(ih_header), "ih", res)
    ih = _table(paths["ih"])
    if len(ih.columns) >= 3:
        if ih.iloc[:, 0].duplicated().any():
            res.errors.append("ih has repeated ArrayIDs")
        headers = set(ih.iloc[:, 1])
        no_ih = [s for s in samples if s not in headers]
        no_expr = [h for h in ih.iloc[:, 1] if h not in set(samples)]
        if no_ih:
            res.errors.append("{} expr columns have no ih row (e.g. {})".format(len(no_ih), _examples(no_ih)))
        if no_expr:
            res.errors.append("{} ih rows point to no expr column (e.g. {})".format(len(no_expr), _examples(no_expr)))
        if set(ih.iloc[:, 0]) != set(s_ids):
            res.errors.append("ih ArrayIDs and survival ArrayIds are not the same set of samples")
        blank = int((ih.iloc[:, 2].str.strip() == "").sum())
        if blank:
            res.errors.append("{} ih rows have a blank ClinicalHeader".format(blank))

    # ---- GEO sample identity
    if repository == "GEO" and not single_cell and samples:
        non_gsm = [s for s in samples if not re.match(r"^GSM\d+$", s)]
        if non_gsm:
            res.errors.append("GEO expr columns must be GSM accessions like the notebook ({} are not, e.g. {}). Map "
                              "each column to its GSM by exact title/characteristic/file name, and leave out columns "
                              "that match no GSM (they have no metadata).".format(len(non_gsm), _examples(non_gsm)))
        elif series_samples:
            outside = [s for s in samples if s not in set(series_samples)]
            if outside:
                res.errors.append("Columns are GSMs from another platform/series: " + _examples(outside))
            absent = [s for s in series_samples if s not in set(samples)]
            if absent:
                res.warnings.append("{} of {} samples on this platform have no expression data in the chosen "
                                    "source (e.g. {})".format(len(absent), len(series_samples), _examples(absent)))
    return res


# ------------------------------------------------------------------ metadata that groups samples
#
# Hegemon is used to compare groups of samples. A metadata column is only useful for that
# if it splits the samples into real groups. `title` (different for every sample) and
# `organism` (the same for every sample) are always filled but group nothing, so a check
# that merely asks "is there metadata?" passes datasets that cannot be plotted by anything.

MISSING_VALUES = {"", "na", "n/a", "nan", "none", "null", "-", "--", "?", "unknown", "not available",
                  "not applicable", "not collected", "missing", "not provided", "not reported", "not specified"}
TECHNICAL_FIELDS = {
    "title", "geo_accession", "organism", "taxid", "platform_id", "molecule", "label",
    "submission_date", "last_update_date", "channel_count", "data_row_count", "sample_id", "sample_name",
    "sample_title", "barcode", "arrayid", "array_id", "assay_name", "scan_name", "extract_name",
    "labeled_extract_name", "hybridization_name", "technology_type", "material_type", "protocol_ref",
    "performer", "term_source_ref", "term_accession_number", "array_data_file", "derived_array_data_file",
    "array_design_ref", "instrument_model", "library_strategy", "library_source", "library_selection",
}


def field_stem(name):
    """'c disease_ch1' -> 'disease'; 'c Characteristics[organism]' -> 'organism'; Comment[...] -> 'comment:...'."""
    s = str(name).strip()
    if s[:2] in ("c ", "n "):
        s = s[2:]
    s = s.strip().lower()
    m = re.match(r"^(characteristics|factor ?value|comment|parameter value|unit)\s*\[(.*)\](\.\d+)?$", s)
    if m:
        if m.group(1) in ("comment", "unit"):
            return "comment:" + m.group(2)
        s = m.group(2)
    s = re.sub(r"_ch\d+$", "", s)
    return re.sub(r"[^a-z0-9]+", "_", s).strip("_")


def comparison_field(name, values, n_samples):
    """Counter of groups if this field splits samples into comparable groups, else None.

    Rules: not an identifier/technical field; filled for at least half the samples; at least
    two groups; groups average two or more samples (so not an ID column); and at least two
    groups have two or more samples each (so there is something to compare).
    """
    stem = field_stem(name)
    if (stem in TECHNICAL_FIELDS or stem.startswith("comment:") or stem.startswith("contact")
            or "protocol" in stem):
        return None
    filled = [str(v).strip() for v in values]
    filled = [v for v in filled if v.lower() not in MISSING_VALUES]
    if len(filled) < max(2, 0.5 * n_samples):
        return None
    counts = Counter(filled)
    if len(counts) < 2 or len(counts) > len(filled) / 2.0:
        return None
    if sum(1 for c in counts.values() if c >= 2) < 2:
        return None
    return counts


def why_not_grouping(name, values, n_samples):
    stem = field_stem(name)
    if (stem in TECHNICAL_FIELDS or stem.startswith("comment:") or stem.startswith("contact")
            or "protocol" in stem):
        return "identifier/technical field"
    filled = [str(v).strip() for v in values]
    filled = [v for v in filled if v.lower() not in MISSING_VALUES]
    counts = Counter(filled)
    if len(filled) < max(2, 0.5 * n_samples):
        return "filled for only {} of {} samples".format(len(filled), n_samples)
    if len(counts) < 2:
        return "the same value for every sample ({})".format(next(iter(counts)) if counts else "")
    if len(counts) > len(filled) / 2.0:
        return "{} different values for {} samples (an identifier, not a group)".format(len(counts), len(filled))
    return "no two groups have 2+ samples"


def groups_text(counts, top=4):
    items = counts.most_common()
    text = " / ".join("{} {}".format(v[:40], c) for v, c in items[:top])
    return text + (" / +{} more".format(len(items) - top) if len(items) > top else "")


def grouping_fields(run_dir, prefix):
    """[(column, Counter)] for every survival column that groups the samples."""
    surv = _table(Path(run_dir) / "{}-survival.txt".format(prefix))
    n = len(surv)
    out = []
    for col in surv.columns[3:]:
        if str(col).startswith(("c ", "n ")):
            counts = comparison_field(col, surv[col].tolist(), n)
            if counts:
                out.append((col, counts))
    return out


def describe_fields(run_dir, prefix, limit=15):
    surv = _table(Path(run_dir) / "{}-survival.txt".format(prefix))
    n = len(surv)
    lines = []
    for col in list(surv.columns[3:])[:limit]:
        lines.append("  {}: {}".format(col, why_not_grouping(col, surv[col].tolist(), n)))
    if len(surv.columns) - 3 > limit:
        lines.append("  ... and {} more columns".format(len(surv.columns) - 3 - limit))
    return "\n".join(lines) or "  (no metadata columns at all)"


def expected_metadata_fields(meta_source, samples):
    """Metadata columns the source holds for these samples, named the way the notebooks name them.

    GEO: 'c <key>_chN' for each characteristics key and 'c source_name_chN' (as geo_survival_ih).
    ArrayExpress: 'c ' + every Characteristics[...] / Factor Value[...] SDRF column (as the E-MTAB notebook).
    """
    if not meta_source:
        return []
    expected = OrderedDict()
    kind = meta_source.get("kind")
    if kind == "GEO" and meta_source.get("gse") is not None:
        gse = meta_source["gse"]
        for s in samples:
            gsm = gse.gsms.get(s)
            if gsm is None:
                continue
            for col, values in gsm.metadata.items():
                if "_ch" not in col or "protocol" in col or "label" in col or "channel_count" in col:
                    continue
                values = [str(v) for v in values if str(v).strip()]
                if not values:
                    continue
                if col.startswith("characteristics"):
                    suffix = "_ch" + col.split("_ch")[-1]
                    if ":" in values[0]:
                        for v in values:
                            if ":" in v:
                                key, value = [p.strip() for p in v.split(":", 1)]
                                if value:
                                    expected["c " + key + suffix] = None
                            else:
                                expected["c " + col] = None
                    else:
                        expected["c " + col] = None
                elif col.startswith("source_name"):
                    expected["c " + col] = None
    elif kind == "ArrayExpress" and meta_source.get("sdrf") and Path(str(meta_source["sdrf"])).is_file():
        table = pd.read_csv(str(meta_source["sdrf"]), sep="\t", dtype=str, keep_default_na=False)
        for col in table.columns:
            low = col.lower()
            if low.startswith(("characteristics", "factor value", "factorvalue")):
                if (table[col].astype(str).str.strip() != "").any():
                    expected["c " + col] = None
    return list(expected)


def missing_metadata_fields(run_dir, prefix, meta_source):
    if not meta_source:
        return []
    samples = header_of(Path(run_dir) / "{}-expr.txt".format(prefix))[2:]
    have = set(header_of(Path(run_dir) / "{}-survival.txt".format(prefix)))
    return [f for f in expected_metadata_fields(meta_source, samples) if f not in have]
