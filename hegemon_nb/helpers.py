"""The lab notebook's fixed steps, as functions for the per-dataset script.

The dataset-specific part (which files, which columns, how sample IDs are
spelled) is written by the AI for each dataset. The steps the notebook does the
same way every time live here so they are done identically every time:

  read_geo_soft      GEOparse.get_GEO(...)                    (GSE notebook cell 5)
  geo_survival_ih    survival + ih from GSM metadata           (cells 9-13)
  geo_sample_matrix  GSM VALUE tables merged on ID_REF         (cell 16)
  geo_platform_names Symbol/Definition from the GPL table      (cell 18)
  lab_gene_names     Ensembl ID -> gene name, lab genome files (cell 24)
  cpm_log2 / cpm / log2p1   normalize_total(1e6) + log1p(base=2)   (cell 29, E-MTAB cell 8)
  write_native       to_csv(sep='\\t', index=False) of the three files

Deliberate fixes relative to the notebook are marked FIX below; each one keeps
data the notebook would drop or duplicate.
"""
from collections import OrderedDict
import os
import re

import numpy as np
import pandas as pd

from .geo_soft import read_geo_soft  # noqa: F401  (re-exported for scripts)

GENOME_DIR = os.environ.get("HEGEMON_GENOME_DIR", "/booleanfs2/sahoo/Data/SeqData/genome")
HUMAN_GENOME = "Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.t.txt"
MOUSE_GENOME = "Mus_musculus.GRCm38.94.chr_patch_hapl_scaff.t.txt"


def gsm_ids(gse, platform=None):
    """GSM names in SOFT order, optionally only those on one GPL."""
    out = []
    for name, gsm in gse.gsms.items():
        if platform is None or (gsm.metadata.get("platform_id") or [""])[0] == platform:
            out.append(name)
    return out


def geo_survival_ih(gse, gsm_list=None):
    """Notebook cell 9: one survival row and one ih row per GSM.

    Keeps 'title' and every *_ch* field except protocol/label/channel_count.
    Characteristics like 'age: 65' become column 'c age_ch1'.
    FIX: splits on the first ':' only (the notebook turned 'dose: 5: mg' or 'age:65'
    into 'N/A'), keeps the real channel suffix, never splits titles, and joins
    repeated values with '; ' instead of silently creating extra rows for a sample.
    """
    names = list(gsm_list) if gsm_list is not None else list(gse.gsms)
    rows, ih = [], OrderedDict([("ArrayID", []), ("ArrayHeader", []), ("ClinicalHeader", [])])
    order = OrderedDict([("ArrayId", None), ("time", None), ("status", None)])
    for name in names:
        gsm = gse.gsms[name]
        row = {"ArrayId": name, "time": "", "status": ""}
        ih["ArrayID"].append(name)
        ih["ArrayHeader"].append(name)
        ih["ClinicalHeader"].append((gsm.metadata.get("title") or [name])[0])
        for col, values in gsm.metadata.items():
            if not (col == "title" or "_ch" in col):
                continue
            if "protocol" in col or "label" in col or "channel_count" in col:
                continue
            values = [str(v) for v in values]
            # FIX: the notebook also split titles like 'Liver: patient 3' into a bogus
            # column; only characteristics fields hold 'key: value' pairs.
            if col.startswith("characteristics") and values and ":" in values[0]:
                suffix = "_ch" + col.split("_ch")[-1]
                for v in values:
                    if ":" in v:
                        key, value = [p.strip() for p in v.split(":", 1)]
                        key = "c " + key + suffix
                    else:
                        key, value = "c " + col, v.strip()
                    row[key] = value if key not in row or not row[key] else row[key] + "; " + value
                    order[key] = None
            else:
                key = "c " + col
                row[key] = "; ".join(v.strip() for v in values)
                order[key] = None
        rows.append(row)
    survival = pd.DataFrame(rows, columns=list(order)).fillna("")
    return survival, pd.DataFrame(ih)


def geo_sample_matrix(gse, platform, value_col="VALUE", id_col="ID_REF", fill_missing=0):
    """Notebook cell 16: outer-merge each GSM table (ID_REF, VALUE) into one matrix.

    Returns columns ID, GSM1, GSM2, ... Missing values become 0 like the notebook
    (gse_df.replace(np.NaN, 0)); pass fill_missing=None to keep them blank.
    """
    pieces = []
    for name in gsm_ids(gse, platform):
        table = gse.gsms[name].table
        if table is None or table.empty or id_col not in table.columns or value_col not in table.columns:
            continue
        part = table[[id_col, value_col]].drop_duplicates(subset=[id_col], keep="first")
        pieces.append(part.set_index(id_col)[value_col].rename(name))
    if not pieces:
        raise ValueError("No GSM on {} has a table with columns {} and {}".format(platform, id_col, value_col))
    first = pieces[0].index
    if first.is_unique and all(p.index.dtype == first.dtype and p.index.equals(first) for p in pieces[1:]) \
            and len({p.name for p in pieces}) == len(pieces):
        # Every sample lists the same probes in the same order, so the outer merge is just the
        # columns side by side (what pd.concat returns in that case, without aligning 400+ indexes).
        matrix = pd.DataFrame(OrderedDict((p.name, p.to_numpy()) for p in pieces), index=first)
    else:
        matrix = pd.concat(pieces, axis=1, join="outer", sort=False)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    if fill_missing is not None:
        matrix = matrix.fillna(fill_missing)
    matrix.index.name = "ID"
    return matrix.reset_index()


_SYMBOL_RENAMES = {
    "GeneSymbol": "Symbol", "Gene Symbol": "Symbol", "gene symbol": "Symbol",
    "GENE_SYMBOL": "Symbol", "GENE_SYM": "Symbol", "Gene_Name": "Symbol",
    "geneName": "Definition", "Gene Name": "Definition", "Gene Title": "Definition",
    "Gene Description": "Definition", "Description": "Definition", "DESCRIPTION": "Definition",
    "description": "Definition", "Gene_Desc": "Definition",
}


def geo_platform_names(gpl, symbol_col=None, definition_col=None):
    """Notebook cell 18: return DataFrame [ID, Symbol, Definition] from the GPL table.

    Uses gene_assignment when present, else the notebook's column renames. Pass
    symbol_col/definition_col when the platform uses other names (the notebook
    asks the operator to add them to the rename list).
    """
    table = gpl.table.copy()
    if symbol_col or definition_col:
        if symbol_col:
            table["Symbol"] = table[symbol_col]
        if definition_col:
            table["Definition"] = table[definition_col]
    elif "gene_assignment" in table.columns:
        symbols, titles = [], []
        for value in table["gene_assignment"]:
            value = str(value)
            if value != "---" and value != "nan":
                parts = [y.split(" // ") for y in value.split(" /// ")]
                titles.append(" /// ".join(p[2] if len(p) > 2 else "ERR" for p in parts))
                symbols.append(" /// ".join(p[1] if len(p) > 2 else "ERR" for p in parts))
            else:
                symbols.append("---")
                titles.append("---")
        table["Symbol"], table["Definition"] = symbols, titles
    else:
        table = table.rename(columns=_SYMBOL_RENAMES)
        table = table.loc[:, ~table.columns.duplicated()]
    if "Definition" not in table.columns and "Symbol" in table.columns:
        table["Definition"] = table["Symbol"]
    if "Symbol" not in table.columns:
        raise ValueError(
            "No gene symbol column in the GPL table. Columns are: {}. "
            "Pass symbol_col=... (and definition_col=...).".format(list(gpl.table.columns)))
    return table[["ID", "Symbol", "Definition"]]


def name_from_symbol(symbol, definition):
    """Notebook: Name = Symbol + ':' + Definition (blank when either is missing)."""
    return symbol.astype("object") + ":" + definition.astype("object")


def lab_gene_names(probe_ids):
    """Notebook cell 24: map Ensembl IDs to gene names with the lab's genome tables.

    Returns a list aligned with probe_ids (None when unmapped).
    FIX: version suffixes are stripped before lookup and the table is de-duplicated,
    so versioned IDs map and genes are not repeated once per transcript.
    """
    ids = pd.Series([str(x) for x in probe_ids])
    if ids.str.contains("ENSG").any():
        fname, cols = HUMAN_GENOME, [1, 2]
    elif ids.str.contains("ENST").any():
        fname, cols = HUMAN_GENOME, [0, 2]
    elif ids.str.contains("ENSMUST").any():
        fname, cols = MOUSE_GENOME, [0, 2]
    elif ids.str.contains("ENSMUSG").any():
        fname, cols = MOUSE_GENOME, [1, 2]
    else:
        return [None] * len(ids)
    path = os.path.join(GENOME_DIR, fname)
    if not os.path.isfile(path):
        raise FileNotFoundError("Lab genome annotation not found: " + path)
    table = pd.read_csv(path, sep="\t", header=None, usecols=cols, dtype=str)
    table.columns = ["ProbeID", "Name"]
    table["ProbeID"] = table["ProbeID"].str.split(".").str[0]
    lookup = table.drop_duplicates("ProbeID").set_index("ProbeID")["Name"]
    stripped = ids.str.split(".").str[0]
    return [lookup.get(x) for x in stripped]


def cpm(values):
    """sc.pp.normalize_total(adata, target_sum=1e6) with samples as columns."""
    values = values.apply(pd.to_numeric, errors="coerce").astype("float64")
    totals = values.sum(axis=0)
    totals = totals.where(totals > 0, 1.0)
    return values.divide(totals, axis=1) * 1e6


def log2p1(values):
    """sc.pp.log1p(adata, base=2) and the microarray np.log2(x+1)."""
    values = values.apply(pd.to_numeric, errors="coerce").astype("float64")
    return np.log2(values + 1.0)


def cpm_log2(values):
    """Notebook RNA-seq default: CPM, then log2(x+1)."""
    return log2p1(cpm(values))


_BAD = re.compile(r"[\t\r\n]+")


def _clean(value):
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return _BAD.sub(" ", str(value)).strip()


def write_native(prefix, expr, survival, ih, directory="."):
    """Write PREFIX-expr.txt, PREFIX-survival.txt and PREFIX-ih.txt like the notebooks do."""
    expr = expr.copy()
    if list(expr.columns[:2]) != ["ProbeID", "Name"]:
        raise ValueError("expr must start with columns ProbeID, Name; got {}".format(list(expr.columns[:3])))
    expr["ProbeID"] = [_clean(x) for x in expr["ProbeID"]]
    expr["Name"] = [_clean(x) for x in expr["Name"]]
    expr.columns = ["ProbeID", "Name"] + [_clean(c) for c in expr.columns[2:]]
    survival = survival.copy()
    survival.columns = [_clean(c) for c in survival.columns]
    for col in survival.columns:
        survival[col] = [_clean(x) for x in survival[col]]
    ih = ih.copy()
    if list(ih.columns) != ["ArrayID", "ArrayHeader", "ClinicalHeader"]:
        raise ValueError("ih columns must be ArrayID, ArrayHeader, ClinicalHeader")
    for col in ih.columns:
        ih[col] = [_clean(x) for x in ih[col]]
    paths = {}
    for kind, frame in (("expr", expr), ("survival", survival), ("ih", ih)):
        path = os.path.join(directory, "{}-{}.txt".format(prefix, kind))
        _to_csv(frame, path)
        paths[kind] = path
    return paths


_SHARED = []


def _csv_rows(bounds):
    start, end = bounds
    return _SHARED[0].iloc[start:end].to_csv(sep="\t", index=False, header=False)


def _to_csv(frame, path):
    """frame.to_csv(path, sep='\t', index=False), byte for byte, with big tables written in parallel.

    Turning ~20 million numbers into text is most of the time a big dataset spends here; rows are
    formatted independently, so pieces of rows are formatted by a few processes and joined in order.
    """
    rows, cols = frame.shape
    workers = min(4, os.cpu_count() or 1)
    threshold = int(os.environ.get("HEGEMON_PARALLEL_WRITE_CELLS", "4000000"))
    if rows * cols >= threshold and workers > 1 and rows > 1 and hasattr(os, "fork"):
        try:
            import multiprocessing
            step = -(-rows // (workers * 2))
            bounds = [(a, min(rows, a + step)) for a in range(0, rows, step)]
            _SHARED[:] = [frame]
            try:
                with multiprocessing.get_context("fork").Pool(workers) as pool:
                    parts = pool.map(_csv_rows, bounds)
            finally:
                _SHARED[:] = []
            with open(path, "w", encoding="utf-8", newline="") as fh:
                fh.write(frame.iloc[:0].to_csv(sep="\t", index=False))
                for part in parts:
                    fh.write(part)
            return
        except Exception:
            pass   # fall back to one process
    frame.to_csv(path, sep="\t", index=False)
