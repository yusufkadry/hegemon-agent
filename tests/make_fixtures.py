"""Create small synthetic datasets that exercise the notebook paths (local tests only).

geo/GSE999001_family.soft.gz  microarray series, 6 samples, Affy-style GPL table
    (Gene Symbol / Gene Title columns, a multi-gene probe, an unannotated probe,
    one missing value, a title containing ':' and a 'treatment: drug: 5 mg' field)
ae/counts/*.count + ae/E-MTAB-9999.sdrf.txt  htseq counts with the 5 summary rows
    at the bottom and one sample split across two lane files
genome/Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.t.txt  tiny lab-style table
"""
import gzip
import sys
from pathlib import Path

import numpy as np

base = Path(sys.argv[1])
rng = np.random.RandomState(7)

# ---------------------------------------------------------------- GEO ----
genes = [("1007_s_at", "DDR1", "discoidin domain receptor tyrosine kinase 1"),
         ("200001_at", "GAPDH", "glyceraldehyde-3-phosphate dehydrogenase"),
         ("200002_at", "ACTB", "actin beta"),
         ("200003_s_at", "RPL5 /// SNORA66", "ribosomal protein L5 /// small nucleolar RNA, H/ACA box 66"),
         ("200004_at", "", ""),
         ("AFFX-BioB-5_at", "", "")]
for i in range(24):
    genes.append(("21{:04d}_at".format(i), "GENE{}".format(i + 1), "made-up gene {}".format(i + 1)))
samples = [("GSM900001", "Tumor 1", "tissue: tumor", "age: 61"),
           ("GSM900002", "Tumor 2", "tissue: tumor", "age: 57"),
           ("GSM900003", "Tumor: relapse 3", "tissue: tumor", "age: 70"),
           ("GSM900004", "Normal 1", "tissue: normal", "age: 49"),
           ("GSM900005", "Normal 2", "tissue: normal", ""),
           ("GSM900006", "Normal 3", "tissue: normal", "age: 66")]
lines = ["^DATABASE = GeoMiame", "!Database_name = Gene Expression Omnibus (GEO)",
         "^SERIES = GSE999001", "!Series_title = Synthetic tumor vs normal test series",
         "!Series_summary = Local test fixture.", "!Series_type = Expression profiling by array",
         "!Series_platform_id = GPL999"] + ["!Series_sample_id = " + s[0] for s in samples]
lines += ["^PLATFORM = GPL999", "!Platform_title = Synthetic Affymetrix-style array", "!Platform_organism = Homo sapiens",
          "#ID = Affymetrix Probe Set ID", "#Gene Title = Entrez Gene name", "#Gene Symbol = Entrez Gene symbol",
          "!platform_table_begin", "ID\tGene Title\tGene Symbol"]
lines += ["{}\t{}\t{}".format(pid, title, sym) for pid, sym, title in genes] + ["!platform_table_end"]
for n, (gsm, title, tissue, age) in enumerate(samples):
    lines += ["^SAMPLE = " + gsm, "!Sample_title = " + title, "!Sample_type = RNA",
              "!Sample_source_name_ch1 = human " + tissue.split(": ")[1], "!Sample_organism_ch1 = Homo sapiens",
              "!Sample_characteristics_ch1 = " + tissue]
    if age:
        lines.append("!Sample_characteristics_ch1 = " + age)
    if n == 2:
        lines.append("!Sample_characteristics_ch1 = treatment: drug: 5 mg")
    lines += ["!Sample_molecule_ch1 = total RNA", "!Sample_label_ch1 = biotin",
              "!Sample_extract_protocol_ch1 = standard", "!Sample_platform_id = GPL999",
              "!Sample_data_row_count = {}".format(len(genes)),
              "#ID_REF = ", "#VALUE = RMA log2 signal", "!sample_table_begin", "ID_REF\tVALUE"]
    for j, (pid, _, _) in enumerate(genes):
        value = "" if (n == 4 and j == 7) else "{:.4f}".format(4 + 8 * rng.rand() + (2 if n < 3 and j < 10 else 0))
        lines.append("{}\t{}".format(pid, value))
    lines.append("!sample_table_end")
(base / "geo").mkdir(parents=True, exist_ok=True)
with gzip.open(str(base / "geo" / "GSE999001_family.soft.gz"), "wt") as fh:
    fh.write("\n".join(lines) + "\n")

# ----------------------------------------------------------------- AE ----
symbols = ["GAPDH", "ACTB", "TP53", "MYC", "EGFR", "CD4", "CD8A", "FOXP3"] + ["GENE{}".format(i) for i in range(1, 33)]
summary = ["__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"]
counts_dir = base / "ae" / "counts"
counts_dir.mkdir(parents=True, exist_ok=True)
sdrf = ["Source Name\tCharacteristics[organism]\tCharacteristics[disease]\tCharacteristics[sex]\tComment[ENA_RUN]"]
files = {"S01": ["S01"], "S02": ["S02"], "S03": ["S03.L1", "S03.L2"], "S04": ["S04"], "S05": ["S05"]}
for k, (sample, parts) in enumerate(files.items()):
    for part in parts:
        with open(str(counts_dir / (part + ".count")), "w") as fh:
            fh.write("gene_id\tcount\n")
            for g in symbols:
                fh.write("{}\t{}\n".format(g, int(rng.poisson(50 + 500 * (g in ("GAPDH", "ACTB"))))))
            for s in summary:
                fh.write("{}\t{}\n".format(s, int(rng.poisson(1000))))
    sdrf.append("{}.1\tHomo sapiens\t{}\t{}\tERR00{}".format(
        sample, "psoriasis" if k < 3 else "normal", "female" if k % 2 else "male", k))
(base / "ae" / "E-MTAB-9999.sdrf.txt").write_text("\n".join(sdrf) + "\n")

# ------------------------------------------------------------- genome ----
(base / "genome").mkdir(parents=True, exist_ok=True)
with open(str(base / "genome" / "Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.t.txt"), "w") as fh:
    for tid, gid, name in (("ENST0000001", "ENSG00000111640", "GAPDH"), ("ENST0000002", "ENSG00000111640", "GAPDH"),
                           ("ENST0000003", "ENSG00000075624", "ACTB")):
        fh.write("{}\t{}\t{}\n".format(tid, gid, name))
print("fixtures in", base)
