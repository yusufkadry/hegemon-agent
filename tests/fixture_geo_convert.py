# Stand-in for an AI-written script (local tests only): the GEO notebook's microarray path.
import os
from hegemon_nb.helpers import (read_geo_soft, geo_sample_matrix, geo_platform_names, name_from_symbol,
                                geo_survival_ih, write_native)

gse = read_geo_soft(os.path.join(SOURCE_DIR, PREFIX.split("-")[0] + "_family.soft.gz"))
platform = PREFIX.split("-")[-1]
expr = geo_sample_matrix(gse, platform)            # VALUE is 'RMA log2 signal' -> takeLog False
ann = geo_platform_names(gse.gpls[platform])       # 'Gene Symbol' / 'Gene Title' via the notebook's renames
data = ann.merge(expr, how="right", on="ID")
data.insert(1, "Name", name_from_symbol(data["Symbol"], data["Definition"]))
data = data.drop(columns=["Symbol", "Definition"]).rename(columns={"ID": "ProbeID"})
survival, ih = geo_survival_ih(gse, list(data.columns[2:]))
write_native(PREFIX, data, survival, ih)
print("source: GSM VALUE tables; {} probes x {} samples; no normalization (RMA log2); metadata: {}".format(
    len(data), data.shape[1] - 2, ", ".join(survival.columns[3:])))
