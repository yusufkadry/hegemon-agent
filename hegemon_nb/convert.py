"""The AI writes the dataset-specific cells of the lab notebook; fixed code judges the result.

One AI call per attempt, no separate reviewer. Every attempt's reply, script and
log are kept in the build folder. Failures are fed back as concrete errors.
"""
import ast
import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

from . import checks, ui
from .config import AGENT_DIR, BLUEPRINT_DIR

REASONS = ("NO_PROCESSED_DATA", "NO_SAMPLE_MATCH", "NOT_EXPRESSION", "MISSING_FILES", "NO_GROUPING_METADATA",
           "NOT_HOSTED")
MAX_FETCH_ROUNDS = 3


class Unconvertible(Exception):
    def __init__(self, reason, detail):
        Exception.__init__(self, "{}: {}".format(reason, detail))
        self.reason, self.detail = reason, detail


class GaveUp(Exception):
    pass


INSTRUCTIONS = """You are processing one public gene-expression dataset into native Hegemon files for Professor Debashis Sahoo's lab (BooleanLab, UCSD).

The lab does this by hand in a Jupyter notebook. The real lab notebooks are included below as BLUEPRINTS. Write the code a lab member writes when adapting that notebook to THIS dataset: same procedure and output format, changing only what depends on this dataset's files (which file holds the expression values, how sample IDs are spelled, which metadata table to use).

Write ONE Python 3.8 script. Two variables exist before your code runs:
  SOURCE_DIR  absolute path of the downloaded dataset files (read only)
  PREFIX      the file prefix, e.g. GSE51984-GPL11154 or E-MTAB-7604
The script runs inside the output folder and must write exactly three files there:
  PREFIX-expr.txt, PREFIX-survival.txt, PREFIX-ih.txt
Do not make idx, thr, info, vinfo or bv files: the notebooks' last two cells run unchanged after your script.

OUTPUT FORMAT (from the notebooks)
- expr: tab-separated; ProbeID, Name, then one column per sample. GEO sample columns are GSM accessions. For arrays Name is "Symbol:Definition" (name_from_symbol). For RNA-seq Name is the gene symbol.
- survival: ArrayId, time, status, then one column per metadata field named "c <field>". time and status stay blank unless the study has real survival time/event data.
- ih: ArrayID, ArrayHeader, ClinicalHeader. ArrayID = ArrayHeader = the expr column name. ClinicalHeader = the GEO sample title (ArrayExpress: the sample name).

WHERE THE FILES ARE
SOURCE_DIR holds what the notebook expects the operator to have fetched:
- <GSE>_family.soft.gz            the family SOFT file (GEO metadata and, for microarrays, the sample tables)
- <GSE>_RAW/                      one file per sample, the same layout the notebook gets from
                                  `tar -xf GSE###_RAW.tar`, so glob it the way the notebook does
- other processed files sit directly in SOURCE_DIR (a counts matrix, a series matrix, an SDRF, an extracted folder)
Read what is actually listed in the evidence below; do not assume a file that is not there.

GETTING MORE FILES
The agent downloaded the cheapest set of files that usually works. The evidence also lists every other file the
study offers ("FILES THIS STUDY OFFERS THAT ARE NOT DOWNLOADED") and any download that failed. The study may come
from any website. If the data or the sample metadata you need is not in SOURCE_DIR, do not guess and do not give
up: reply with one line per file and no code,
  FETCH: <exact name from that list, a name pattern such as GSE123_RAW/*_counts.txt.gz, or a public https URL>
The files are downloaded (checked, archives unpacked) and you get the evidence again. Prefer the smallest files
that hold processed values plus the per-sample metadata.

HELPERS: from hegemon_nb.helpers import ... They are the notebook's fixed cells; use them instead of rewriting those cells.
  read_geo_soft(path) -> gse   same object model GEOparse gives the notebook: gse.metadata, gse.gsms[GSM].metadata/.table/.columns, gse.gpls[GPL].table. The family SOFT file is SOURCE_DIR/<GSE>_family.soft.gz
  gsm_ids(gse, platform) -> GSMs on that platform in SOFT order
  geo_survival_ih(gse, gsm_list) -> (survival, ih)   notebook cells 9-13; pass the GSMs in expr column order
  geo_sample_matrix(gse, platform, value_col="VALUE") -> DataFrame [ID, GSM...]   notebook cell 16 (missing -> 0 like the notebook)
  geo_platform_names(gpl, symbol_col=None, definition_col=None) -> DataFrame [ID, Symbol, Definition]   notebook cell 18; pass column names when the platform uses others
  name_from_symbol(symbol, definition) -> Series "Symbol:Definition"
  lab_gene_names(ids) -> gene names for Ensembl gene/transcript IDs from the lab genome tables (notebook cell 24)
  cpm(values), log2p1(values), cpm_log2(values)   notebook normalization (normalize_total(target_sum=1e6), log1p(base=2)); samples are columns
  write_native(PREFIX, expr, survival, ih)   writes the three files like the notebooks' to_csv calls; use it

RULES
1. Keep every gene/probe row of the expression source you choose, and every sample that has both expression values and metadata. Do not average rows or samples.
2. Match samples by exact identifiers only (GSM accession, sample title, characteristic value, file name, SDRF Source Name, a sample sheet's ID column). Never assign data to a sample by position or by guessing. Every expr column needs exactly one survival row and one ih row, and vice versa. Use what matches: data columns that match no sample in the metadata are left out (they cannot be plotted), and samples that have metadata but no data are left out. Print how many samples were matched and how many of each kind were left out, with a few examples. Some columns not matching is never a reason to give up.
3. Never invent metadata. Keep every observed field: for GEO use geo_survival_ih (it keeps every characteristic);
   for ArrayExpress prefix every SDRF column with 'c ' exactly as the E-MTAB notebook does. If SOURCE_DIR also has a
   study-supplied per-sample table (clinical data, sample sheet, xlsx/csv/tsv with one row per sample), join it by
   exact sample IDs and add its columns as 'c <column>'. Hegemon is used to COMPARE GROUPS of samples, so the
   annotations that split samples (disease, treatment, tissue, genotype, time, response) are the point of the file.
   A file whose only metadata is titles and organism is useless. Unknown values stay blank.
4. Read identifier and metadata columns as text (dtype=str, keep_default_na=False) so leading zeros and literal NA survive.
5. Normalize like the notebooks: RNA-seq counts, and also RPKM/FPKM/TPM values, get cpm_log2 (the GSE51984 notebook used norm='cpm', takeLog=True on RPKM files). Values already on a log scale stay as they are. Microarray VALUE columns stay as they are when already log2 (RMA, 'log2' in the VALUE description, values mostly under 20); otherwise log2p1. Missing expression values become 0, as in the notebooks.
6. RNA-seq count files: drop htseq summary rows (no_feature, ambiguous, too_low_aQual, not_aligned, alignment_not_unique, with or without leading underscores) and sum files that belong to the same sample, as the notebooks do.
7. Ensembl IDs: strip version suffixes from ProbeID and fill Name with lab_gene_names(), like the notebook; leave Name blank where unmapped.
8. ProbeIDs must be unique. If a table uses gene symbols as IDs and some repeat, keep every row and make the IDs unique (SYMBOL, SYMBOL_2, ...) with Name = the symbol.
9. Python 3.8, pandas 2.0.3 (no DataFrame.append: use pd.concat), numpy 1.24. Installed: <<PACKAGES>>. scanpy and GEOparse are not installed; the helpers replace them.
10. Read only from SOURCE_DIR and <<GENOME>>. Write only in the current folder. No network, no subprocess, no os.system, no deleting or moving files.
11. End by printing: the expression source, samples x genes, the normalization applied, and the metadata columns.

The blueprints contain lab quirks (DataFrame.append, uninitialized variables, merges that duplicate rows). Follow their procedure, not their bugs; the helpers already fix those cells.

SAMPLE METADATA
- GEO: the family SOFT file holds every sample's metadata. read_geo_soft + geo_survival_ih always produce it.
- ArrayExpress: PREFIX.sdrf.txt is the sample table and is always downloaded.
- Any other site: the metadata is in a sample sheet, a series matrix, a clinical table or the file names. If it is
  not in SOURCE_DIR, FETCH it from the study's file list.

Reply in ONE of three ways:
(a) FETCH lines (see GETTING MORE FILES) when a file you need is not in SOURCE_DIR.
(b) 2-4 short lines stating your decisions (source file, how samples are matched and how many, normalization),
    then the complete script in one ```python block.
(c) Only if the dataset truly cannot be converted, one line and no code:
    UNCONVERTIBLE: <REASON>: <one sentence citing the evidence>
    where REASON is one of NO_PROCESSED_DATA (only raw reads/images, no processed values in SOURCE_DIR or in the
    study's file list), NO_SAMPLE_MATCH (fewer than two data columns can be matched to samples), NOT_EXPRESSION
    (not a gene/probe-by-sample table), NOT_HOSTED (the file that holds the data is listed but cannot be
    downloaded, and no other file holds it), NO_GROUPING_METADATA (every source was checked and the samples carry
    no annotation that splits them into groups). A missing file is never a reason: FETCH it.
"""

FORBIDDEN_IMPORTS = {"subprocess", "socket", "urllib", "urllib2", "urllib3", "requests", "http", "httpx",
                     "ftplib", "smtplib", "telnetlib", "paramiko", "ctypes", "multiprocessing", "importlib",
                     "shutil", "asyncio", "webbrowser", "pty", "signal", "pickle"}
OS_FORBIDDEN = {"system", "popen", "spawnl", "spawnle", "spawnlp", "spawnv", "spawnve", "spawnvp", "execv",
                "execve", "execl", "execlp", "execvp", "fork", "kill", "remove", "unlink", "rmdir", "removedirs",
                "rename", "renames", "replace", "chmod", "chown", "symlink", "link", "truncate", "chdir"}
ANY_FORBIDDEN = {"unlink", "rmtree", "rmdir", "removedirs", "symlink_to", "hardlink_to", "system", "popen"}
BUILTINS = {"eval", "exec", "compile", "__import__", "breakpoint", "input"}


def guard(code, allowed_roots):
    """Screen AI code before it runs on the shared server. Best effort, not a sandbox."""
    try:
        tree = ast.parse(code)
    except SyntaxError as err:
        return ["Python syntax error: {} (line {})".format(err.msg, err.lineno)]
    problems = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".")[0] in FORBIDDEN_IMPORTS:
                    problems.add("import {} is not allowed".format(alias.name))
        elif isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] in FORBIDDEN_IMPORTS:
                problems.add("import from {} is not allowed".format(node.module))
        elif isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Name) and func.id in BUILTINS:
                problems.add("{}() is not allowed".format(func.id))
            elif isinstance(func, ast.Attribute):
                if isinstance(func.value, ast.Name) and func.value.id == "os" and func.attr in OS_FORBIDDEN:
                    problems.add("os.{}() is not allowed".format(func.attr))
                elif func.attr in ANY_FORBIDDEN:
                    problems.add(".{}() is not allowed".format(func.attr))
        elif isinstance(node, ast.Constant) and isinstance(node.value, str):
            s = node.value
            if len(s) > 1 and s.startswith("/") and "\n" not in s and not any(s.startswith(r) for r in allowed_roots):
                problems.add("absolute path outside the dataset folders: {}".format(s[:80]))
            if "../" in s:
                problems.add("relative path leaving the folder: {}".format(s[:80]))
    return sorted(problems)


def installed_packages():
    names = ["pandas", "numpy", "scipy", "openpyxl", "xlrd", "h5py", "anndata", "pyarrow", "rdata", "scanpy",
             "GEOparse", "statsmodels"]
    found = []
    for name in names:
        if importlib.util.find_spec(name) is not None:
            try:
                from importlib.metadata import version
                found.append("{} {}".format(name, version(name)))
            except Exception:
                found.append(name)
    return ", ".join(found)


def blueprints():
    parts = []
    for path in sorted(BLUEPRINT_DIR.glob("*.txt")):
        parts.append("\n=== BLUEPRINT: {} ===\n{}".format(path.name, path.read_text(errors="replace")))
    return "\n".join(parts)


def parse_fetch(reply):
    """The files the AI asked for: one per 'FETCH: ...' line."""
    items = []
    for line in reply.splitlines():
        m = re.match(r"^\s*[-*]?\s*FETCH\s*:\s*(.+?)\s*$", line, re.I)
        if m:
            items.append(m.group(1).strip().strip("`'\""))
    return [i for i in items if i]


def extract(reply):
    """Returns (code, decision_lines, unconvertible_or_None)."""
    blocks = re.findall(r"```(?:python|py)?[ \t]*\n(.*?)```", reply, re.S)
    if blocks:
        code = max(blocks, key=len)
        before = reply.split("```", 1)[0].strip()
        return code, [l.strip() for l in before.splitlines() if l.strip()][:6], None
    m = re.search(r"UNCONVERTIBLE:\s*([A-Z_]+)\s*:\s*(.+)", reply)
    if m and m.group(1) in REASONS:
        return None, [], (m.group(1), m.group(2).strip())
    return None, [], None


PREAMBLE = ("# Written by hegemon-agent's AI for one dataset (the notebook cells a lab member edits).\n"
            "SOURCE_DIR = {source!r}\n"
            "PREFIX = {prefix!r}\n"
            "import sys as _sys; _sys.path.insert(0, {agent!r}); import warnings as _w; "
            "_w.filterwarnings('ignore', message=r'Columns \\(.*\\) have mixed types')\n")
PREAMBLE_LINES = 4


def _fix_line_numbers(text, script_name):
    def shift(m):
        return "{}, line {}".format(m.group(1), max(1, int(m.group(2)) - PREAMBLE_LINES))
    return re.sub(r'(File "[^"]*{}")?, line (\d+)'.format(re.escape(script_name)),
                  lambda m: shift(m) if m.group(1) else m.group(0), text)


def run_script(code, script_path, run_dir, source_dir, prefix, cfg, log_path):
    script_path.write_text(PREAMBLE.format(source=str(source_dir), prefix=prefix, agent=str(AGENT_DIR)) + code)
    env = {
        "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
        "HOME": os.environ.get("HOME", ""),
        "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8",
        "PYTHONNOUSERSITE": "1", "PYTHONUNBUFFERED": "1",
        "HEGEMON_GENOME_DIR": str(cfg.genome_dir),
        "OMP_NUM_THREADS": "4", "OPENBLAS_NUM_THREADS": "4", "MKL_NUM_THREADS": "4",
    }  # deliberately no API keys in the script's environment
    for name in ("HEGEMON_PARALLEL_WRITE_CELLS",):
        if name in os.environ:
            env[name] = os.environ[name]
    with open(str(log_path), "w") as log, ui.Heartbeat("conversion script running", every=30):
        try:
            rc = subprocess.run([sys.executable, str(script_path)], cwd=str(run_dir), env=env,
                                stdout=log, stderr=subprocess.STDOUT, timeout=cfg.script_timeout).returncode
        except subprocess.TimeoutExpired:
            rc = "timeout after {:.0f}s".format(cfg.script_timeout)
    with open(str(log_path), "r", errors="replace") as fh:
        text = fh.read()
    return rc, _fix_line_numbers(text[-8000:], script_path.name)


BUILTIN_GEO_ARRAY = """# Built into the agent: the GSE notebook's microarray cells (9-20), run with no AI.
import os
import numpy as np
import pandas as pd
from hegemon_nb.helpers import (read_geo_soft, geo_sample_matrix, geo_platform_names, name_from_symbol,
                                geo_survival_ih, write_native, log2p1)

accession, platform = PREFIX.split("-", 1)
gse = read_geo_soft(os.path.join(SOURCE_DIR, accession + "_family.soft.gz"))

# cell 16: every GSM's VALUE table merged on ID_REF (missing -> 0, as in the notebook)
expr = geo_sample_matrix(gse, platform)
values = expr.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
flat = values.to_numpy(dtype="float64")
flat = flat[np.isfinite(flat)]
top = float(np.percentile(flat, 99)) if flat.size else 0.0
low = float(flat.min()) if flat.size else 0.0
# takeLog: intensities get log2(x+1); log values and log-ratios stay as they are
if low >= 0 and top >= 100:
    values = log2p1(values)
    transform = "log2(x+1)"
else:
    transform = "none (already log scale)"
expr = pd.concat([expr[["ID"]], values], axis=1)

# cell 18: gene symbol and definition from the platform table
ann = geo_platform_names(gse.gpls[platform]).drop_duplicates("ID")
ann["ID"] = ann["ID"].astype(str)
expr["ID"] = expr["ID"].astype(str)

# cell 20: Name = Symbol:Definition, ProbeID first, every expr row kept
data = ann.merge(expr, how="right", on="ID")
data.insert(1, "Name", name_from_symbol(data["Symbol"], data["Definition"]))
data = data.drop(columns=["Symbol", "Definition"]).rename(columns={"ID": "ProbeID"})

# cells 9-13: survival and ih from the GSM metadata, in expr column order
survival, ih = geo_survival_ih(gse, list(data.columns[2:]))
write_native(PREFIX, data, survival, ih)
print("built-in notebook conversion: GSM VALUE tables from the SOFT file; {} probes x {} samples; "
      "normalization: {}".format(len(data), data.shape[1] - 2, transform))
"""

GROUPS_FEEDBACK = (
    "The files pass every format check, but no metadata column splits the samples into groups that Hegemon can "
    "compare. Each column is an identifier, constant, or mostly unique:\n{}\nLook in SOURCE_DIR, and in the "
    "study's list of files not downloaded (FETCH it), for a study-supplied per-sample table (clinical data, sample "
    "sheet, series matrix, xlsx/csv/tsv with one row per sample) and join its columns by exact sample IDs as "
    "'c <column>'. If no such table exists anywhere, reply exactly: UNCONVERTIBLE: NO_GROUPING_METADATA: <what "
    "annotation the samples do have>.")


def convert(ds, prefix, work, run_dir, source_dir, evidence_text, series_samples, llm, cfg, script_file=None,
            meta_source=None, builtin=None, fetcher=None):
    """Returns (code, result, output, attempts). attempts == 0 means the built-in script did it, no AI.

    fetcher(items) -> (report, new_evidence): downloads the files the AI asks for with FETCH lines and
    rebuilds the evidence, so missing or failed files are handled inside this loop."""
    work, run_dir = Path(work), Path(run_dir)
    logs = work / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    instructions = (INSTRUCTIONS.replace("<<PACKAGES>>", installed_packages())
                    .replace("<<GENOME>>", str(cfg.genome_dir)) + blueprints())
    base = "PREFIX = {}\nSOURCE_DIR = {}\n\n{}".format(prefix, source_dir, evidence_text)
    genome_ok = any((cfg.genome_dir / n).is_file() for n in (
        "Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.t.txt", "Mus_musculus.GRCm38.94.chr_patch_hapl_scaff.t.txt"))

    def judge(code, label):
        """Run one script and apply every check. Returns (verdict, info)."""
        problems = guard(code, [str(source_dir), str(cfg.genome_dir), str(run_dir), str(work)])
        if problems:
            text = "The script was not run, because: " + "; ".join(problems)
            ui.warn(text)
            return "rejected", text
        for kind in ("expr", "survival", "ih"):
            out = run_dir / "{}-{}.txt".format(prefix, kind)
            if out.exists():
                out.unlink()
        script = work / "convert-{}.py".format(label)
        rc, output = run_script(code, script, run_dir, source_dir, prefix, cfg, logs / "convert-{}.log".format(label))
        if rc != 0:
            last = [l for l in output.strip().splitlines() if l.strip()][-1:] or ["(no output)"]
            ui.warn("script stopped ({}): {}".format(rc, last[0][:200]))
            return "crashed", "The script stopped with exit status {}. Output:\n{}".format(rc, output)
        result = checks.check_core(run_dir, prefix, ds.repository, series_samples, genome_ok=genome_ok)
        if not result.ok:
            for err in result.errors:
                ui.warn(err)
            return "checks", "The script ran, but its files failed the checks:\n{}\n\nScript output (end):\n{}".format(
                result.text(), output[-3000:])
        # Every metadata field the source holds must reach survival.txt.
        missing = checks.missing_metadata_fields(run_dir, prefix, meta_source)
        if missing:
            ui.warn("survival.txt dropped {} metadata field(s) the source has: {}".format(
                len(missing), ", ".join(missing[:8]) + (" ..." if len(missing) > 8 else "")))
            return "fields", ("The files pass the format checks, but survival.txt is missing metadata fields that the "
                              "source holds for these samples: {}. Keep every field, named exactly like the notebooks "
                              "name them (for GEO call geo_survival_ih with the expr GSMs; for ArrayExpress prefix every "
                              "SDRF column with 'c ').".format(", ".join(missing)))
        # And at least one field has to split the samples into groups Hegemon can compare.
        groups = checks.grouping_fields(run_dir, prefix)
        if not groups and os.environ.get("HEGEMON_ALLOW_UNGROUPED") != "1":
            return "groups", checks.describe_fields(run_dir, prefix)
        result.stats["groups"] = [(col, checks.groups_text(counts)) for col, counts in groups]
        (work / "convert.py").write_text(script.read_text())
        return "ok", (result, output)

    feedback, previous, asked_for_groups, failed_before = None, None, False, False

    # The notebook's fixed cells first, when they can do the whole job: no AI time or cost.
    if builtin and not script_file and os.environ.get("HEGEMON_BUILTIN", "1") != "0":
        ui.step("Built-in conversion: the notebook's fixed microarray cells (no AI)")
        verdict, info = judge(builtin, "builtin")
        if verdict == "ok":
            ui.ok("converted by the notebook's own cells; no AI call needed")
            return builtin, info[0], info[1], 0
        if verdict == "groups":
            if llm is None:
                raise Unconvertible("NO_GROUPING_METADATA", "no metadata column splits the samples into groups:\n" + info)
            asked_for_groups, failed_before = True, True
            feedback, previous = GROUPS_FEEDBACK.format(info), builtin
            ui.warn("no metadata column splits the samples into groups; asking the AI to look for a sample table")
        else:
            ui.info("the built-in conversion did not pass ({}); the AI takes over".format(verdict))

    explicit_effort = (os.environ.get("HEGEMON_REASONING") or "").strip() or None
    attempt, fetch_rounds, hinted, why_again = 0, 0, False, ""
    while attempt < cfg.attempts:
        if script_file and attempt == 0:
            attempt += 1
            ui.step("Conversion attempt {} of {}".format(attempt, cfg.attempts))
            code = Path(script_file).read_text()
            ui.info("using the provided script " + str(script_file))
        else:
            if llm is None:
                raise GaveUp("The provided script failed and no AI is available to repair it:\n" + (feedback or ""))
            ui.step("Conversion attempt {} of {}{}".format(attempt + 1, cfg.attempts, why_again))
            why_asked, why_again = bool(why_again), ""
            text = base
            if feedback:
                if previous:
                    text += "\n\nPREVIOUS SCRIPT:\n```python\n{}\n```".format(previous)
                label = "WHAT WENT WRONG" if failed_before and not why_asked else "WHAT HAPPENED"
                text += "\n\n{}:\n{}\n\n{}".format(label, feedback[-9000:], (
                    "Fix the script for this same dataset and return the complete corrected script, or FETCH what "
                    "you need." if previous else "Write the conversion script now, or FETCH what you still need."))
            # A first try at medium reasoning is about twice as fast; repairs think harder.
            effort = explicit_effort or ("high" if failed_before else "medium")
            reply = llm.ask(instructions, text, effort=effort, label="AI writing the conversion")
            (work / "ai-reply-{}.txt".format(attempt + fetch_rounds + 1)).write_text(reply)
            code, decisions, unconvertible = extract(reply)
            wanted = [] if code else parse_fetch(reply)
            if unconvertible and unconvertible[0] == "MISSING_FILES" and fetcher is not None and not hinted:
                hinted = True
                feedback = ("You replied MISSING_FILES ({}). A missing file is not a reason to stop: the study's file "
                            "list is in the evidence. Ask for the file with FETCH lines, or convert with what is in "
                            "SOURCE_DIR.".format(unconvertible[1]))
                ui.warn("the AI said a file is missing; it was told to ask for it instead")
                why_again = " (asked again)"
                continue
            if unconvertible:
                raise Unconvertible(*unconvertible)
            if wanted:
                if fetcher is None or fetch_rounds >= MAX_FETCH_ROUNDS:
                    feedback = ("No more files can be downloaded for this build. Convert with what is in SOURCE_DIR, "
                                "or reply UNCONVERTIBLE.")
                    ui.warn("the AI asked for more files ({}); no more downloads allowed".format(", ".join(wanted[:3])))
                    fetch_rounds += 1
                    if fetch_rounds > MAX_FETCH_ROUNDS + 1:
                        raise GaveUp("The AI kept asking for files instead of converting: " + ", ".join(wanted[:5]))
                    continue
                fetch_rounds += 1
                ui.step("The AI asked for {} more file(s)".format(len(wanted)))
                for item in wanted[:8]:
                    ui.info("FETCH " + item)
                report, new_evidence = fetcher(wanted)
                base = "PREFIX = {}\nSOURCE_DIR = {}\n\n{}".format(prefix, source_dir, new_evidence)
                feedback = "Result of your FETCH request:\n" + report
                for line in report.splitlines()[:12]:
                    (ui.warn if "FAILED" in line or "not fetched" in line or "not in the list" in line else ui.info)(line)
                why_again = " (with the files it asked for)"
                continue
            attempt += 1
            if not code:
                feedback = "Your reply contained no ```python block, no FETCH line and no valid UNCONVERTIBLE line."
                failed_before = True
                ui.warn(feedback)
                continue
            for line in decisions:
                ui.info("AI: " + line)
        verdict, info = judge(code, "attempt{}".format(attempt))
        if verdict == "ok":
            return code, info[0], info[1], attempt
        failed_before = True
        if verdict == "groups":
            if llm is None or asked_for_groups:
                raise Unconvertible("NO_GROUPING_METADATA", "no metadata column splits the samples into groups:\n" + info)
            asked_for_groups = True
            ui.warn("no metadata column splits the samples into groups; asking the AI to look for a sample table")
            feedback = GROUPS_FEEDBACK.format(info)
        else:
            feedback = info
        previous = code
    raise GaveUp(feedback or "no attempts made")
