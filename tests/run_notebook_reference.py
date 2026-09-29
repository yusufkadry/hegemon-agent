"""Run the lab's own notebook cells (from blueprints/) on the local fixtures -> reference files.

Only changes to the notebook code: dataset names and paths, GEOparse reads the local
SOFT file instead of downloading, Jupyter-only lines (%xmode, !cmd, the IPython
import) are dropped, display() is a no-op, a DataFrame.append shim is added for
pandas >= 2 (the lab's Jupyter has an older pandas), and for the microarray run the
scanpy/anndata imports are dropped (only the RNA-seq cells use them).

usage: python run_notebook_reference.py geo|ae FIXTURES_DIR OUT_DIR
"""
import os
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

BLUEPRINTS = Path(__file__).resolve().parents[1] / "blueprints"


def cells(path):
    parts = re.split(r"===== CELL (\d+) \((CODE|MARKDOWN)\) =====", Path(path).read_text())
    return {int(parts[i]): parts[i + 2] for i in range(1, len(parts), 3) if parts[i + 1] == "CODE"}


def clean(code):
    return "\n".join(l for l in code.splitlines() if not l.strip().startswith(("%", "!")))


if not hasattr(pd.DataFrame, "append"):
    def _append(self, other, sort=False, **kw):
        return pd.concat([self, other], sort=sort)
    pd.DataFrame.append = _append

which, fixtures, out = sys.argv[1], Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve()
ns = {"display": lambda *a, **k: None, "__name__": "__notebook__"}

if which == "geo":
    c = cells(BLUEPRINTS / "GSE51984-jupyter-process.txt")
    work = out / "GSE999001-GPL999"
    work.mkdir(parents=True, exist_ok=True)
    os.chdir(str(work))
    soft = str(fixtures / "geo" / "GSE999001_family.soft.gz")
    subs = [("'GSE51984'", "'GSE999001'"), ("experimentType = 'rna'", "experimentType = 'm'"),
            ("norm = 'cpm'", "norm = 'None'"), ("takeLog = 'True'", "takeLog = 'False'"),
            ("GEOparse.get_GEO(geo=accessionID)", "GEOparse.get_GEO(filepath={!r}, silent=True)".format(soft))]
    for n in (3, 5, 7, 9, 11, 13, 16, 18, 20, 31):
        code = "\n".join(l for l in clean(c[n]).splitlines()
                         if not re.match(r"\s*(from IPython|import scanpy|from anndata)", l))
        for a, b in subs:
            code = code.replace(a, b)
        exec(compile(code, "notebook cell {}".format(n), "exec"), ns)
elif which == "ae":
    c = cells(BLUEPRINTS / "E-MTAB-7604-make-Hegemon-files.txt")
    work = out / "E-MTAB-9999"
    if work.exists():
        shutil.rmtree(str(work))
    shutil.copytree(str(fixtures / "ae"), str(work))
    os.chdir(str(work))
    for n in (1, 6, 7, 8, 9, 10, 11, 12):
        exec(compile(clean(c[n]).replace("E-MTAB-7604", "E-MTAB-9999"), "notebook cell {}".format(n), "exec"), ns)
    shutil.rmtree("counts")
    os.remove("E-MTAB-9999.sdrf.txt")
print("notebook reference written to", work, sorted(os.listdir(str(work))))
