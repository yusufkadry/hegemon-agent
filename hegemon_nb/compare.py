"""Compare an agent build with a folder the lab made by hand with the notebook."""
import csv
from pathlib import Path

import numpy as np
import pandas as pd


def _find(folder, kind):
    hits = sorted(Path(folder).glob("*-{}.txt".format(kind)))
    if len(hits) != 1:
        raise SystemExit("Expected one *-{}.txt in {}, found {}".format(kind, folder, len(hits)))
    return hits[0]


def _read(path, **kw):
    return pd.read_csv(str(path), sep="\t", dtype=str, keep_default_na=False, na_filter=False,
                       quoting=csv.QUOTE_NONE, **kw)


def _show(values, n=5):
    values = sorted(values)
    return ", ".join(values[:n]) + (" ..." if len(values) > n else "")


def compare(agent_dir, lab_dir, tol=1e-3):
    lines, same = [], True
    a_expr, l_expr = _read(_find(agent_dir, "expr")), _read(_find(lab_dir, "expr"))
    a_s, l_s = list(a_expr.columns[2:]), list(l_expr.columns[2:])
    common_s = [s for s in l_s if s in set(a_s)]
    lines.append("samples: agent {}, lab {}, shared {}".format(len(a_s), len(l_s), len(common_s)))
    if set(a_s) != set(l_s):
        same = False
        if set(l_s) - set(a_s):
            lines.append("  only in lab:   " + _show(set(l_s) - set(a_s)))
        if set(a_s) - set(l_s):
            lines.append("  only in agent: " + _show(set(a_s) - set(l_s)))

    a_expr = a_expr.drop_duplicates(a_expr.columns[0]).set_index(a_expr.columns[0])
    l_expr = l_expr.drop_duplicates(l_expr.columns[0]).set_index(l_expr.columns[0])
    common_p = l_expr.index.intersection(a_expr.index)
    lines.append("probes: agent {}, lab {}, shared {}".format(len(a_expr), len(l_expr), len(common_p)))
    if len(common_p) != len(a_expr) or len(common_p) != len(l_expr):
        same = False
        lines.append("  only in lab:   " + _show(set(l_expr.index) - set(a_expr.index)))
        lines.append("  only in agent: " + _show(set(a_expr.index) - set(l_expr.index)))
    if len(common_p):
        names_equal = (a_expr.loc[common_p, "Name"].str.strip() == l_expr.loc[common_p, "Name"].str.strip()).mean()
        lines.append("Name identical for {:.1%} of shared probes".format(names_equal))
        same &= names_equal > 0.999
    if len(common_p) and common_s:
        a = a_expr.loc[common_p, common_s].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        b = l_expr.loc[common_p, common_s].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        both = ~np.isnan(a) & ~np.isnan(b)
        diff = np.abs(a - b)[both]
        blank_mismatch = int((np.isnan(a) != np.isnan(b)).sum())
        within = float((diff <= tol).mean()) if diff.size else 1.0
        lines.append("values: {:.2%} within {} (max difference {:.3g}); blank in one file only: {}".format(
            within, tol, float(diff.max()) if diff.size else 0.0, blank_mismatch))
        same &= within > 0.9999 and blank_mismatch == 0

    a_sv, l_sv = _read(_find(agent_dir, "survival")), _read(_find(lab_dir, "survival"))
    a_cols, l_cols = set(a_sv.columns[3:]), set(l_sv.columns[3:])
    lines.append("survival columns: agent {}, lab {}, shared {}".format(len(a_cols), len(l_cols), len(a_cols & l_cols)))
    if a_cols != l_cols:
        same = False
        if l_cols - a_cols:
            lines.append("  only in lab:   " + _show(l_cols - a_cols))
        if a_cols - l_cols:
            lines.append("  only in agent: " + _show(a_cols - l_cols))
    a_sv = a_sv.drop_duplicates(a_sv.columns[0]).set_index(a_sv.columns[0])
    l_sv = l_sv.drop_duplicates(l_sv.columns[0]).set_index(l_sv.columns[0])
    ids = l_sv.index.intersection(a_sv.index)
    for col in sorted(a_cols & l_cols):
        eq = (a_sv.loc[ids, col].str.strip() == l_sv.loc[ids, col].str.strip()).mean() if len(ids) else 1.0
        if eq < 1:
            same = False
            lines.append("  '{}' differs for {:.0%} of samples".format(col, 1 - eq))
    a_ih, l_ih = _read(_find(agent_dir, "ih")).set_index("ArrayID"), _read(_find(lab_dir, "ih")).set_index("ArrayID")
    ids = l_ih.index.intersection(a_ih.index)
    if len(ids):
        eq = (a_ih.loc[ids, "ClinicalHeader"] == l_ih.loc[ids, "ClinicalHeader"]).mean()
        lines.append("ih ClinicalHeader identical for {:.1%} of shared samples".format(eq))
        same &= eq == 1
    lines.append("RESULT: " + ("MATCH - same samples, probes, values and metadata as the lab's files" if same
                               else "DIFFERENCES - see above"))
    return same, "\n".join(lines)
