# Stand-in for an AI-written script (local tests only): the E-MTAB notebook adapted to E-MTAB-9999.
import glob
import os
import pandas as pd
from hegemon_nb.helpers import cpm_log2, write_native

HTSEQ = ["no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"]
df = None
for f in sorted(glob.glob(os.path.join(SOURCE_DIR, "counts", "*.count"))):
    arr = os.path.basename(f).split(".")[0]          # notebook: f.split('/')[1].split('.')[0]
    part = pd.read_csv(f, sep="\t", header=0, names=["ProbeID", arr], dtype={"ProbeID": str})
    part = part[~part["ProbeID"].str.lstrip("_").isin(HTSEQ)].set_index("ProbeID")
    if df is None:
        df = part
    elif arr in df.columns:
        df[arr] = df[arr].add(part[arr], fill_value=0)  # notebook sums lane files of one sample
    else:
        df = df.join(part, how="outer")
df = df.fillna(0)
expr = cpm_log2(df)                                  # notebook: normalize_total(1e6) + log1p(base=2)
expr.insert(0, "Name", df.index)                      # notebook: Name = ProbeID
expr = expr.reset_index()
surv = pd.read_csv(os.path.join(SOURCE_DIR, "E-MTAB-9999.sdrf.txt"), sep="\t", dtype=str, keep_default_na=False)
surv["Source Name"] = surv["Source Name"].str.split(".").str[0]
surv = surv.rename(columns={"Source Name": "ArrayId"})
surv.columns = ["c " + x if x != "ArrayId" else x for x in surv.columns]
surv.insert(1, "time", "")
surv.insert(2, "status", "")
ih = pd.DataFrame({"ArrayID": surv["ArrayId"], "ArrayHeader": surv["ArrayId"], "ClinicalHeader": surv["ArrayId"]})
write_native(PREFIX, expr, surv, ih)
print("source: counts/*.count; {} genes x {} samples; CPM + log2(x+1); metadata: {}".format(
    len(expr), expr.shape[1] - 2, ", ".join(surv.columns[3:])))
