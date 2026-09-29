"""Local tests for hegemon-agent. Run: python tests/run_local_tests.py

What runs for real here: the SOFT reader (vs GEOparse), helpers, checks, the screening
of AI code, the notebooks' idx code, whole builds from local files, publishing to
Hegemon's real public PHP (github.com/sahoo00/Hegemon) served by `php -S`, and the
AI repair loop against a local mock of the OpenAI endpoint.
What does NOT run here: the lab's jupyter_gse_processing (a FAKE stand-in writes
thr/info/vinfo/bv), real OpenAI calls, and real GEO/ArrayExpress downloads.
"""
import json
from collections import OrderedDict
import os
import shutil
import re
import socket
import subprocess
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
AGENT = HERE.parent
sys.path.insert(0, str(AGENT))
T = Path(os.environ.get("HEGEMON_TEST_DIR", "/tmp/hgtest"))
PY = sys.executable
HEGEMON_REPO = Path(os.environ.get("HEGEMON_PUBLIC_REPO", "/tmp/sahoo_Hegemon"))
RESULTS = []


class Skip(Exception):
    pass


def run(name, fn):
    try:
        fn()
        RESULTS.append((name, "PASS", ""))
        print("PASS  " + name, flush=True)
    except Skip as s:
        RESULTS.append((name, "SKIP", str(s)))
        print("SKIP  {} ({})".format(name, s), flush=True)
    except Exception as e:
        RESULTS.append((name, "FAIL", repr(e)))
        print("FAIL  {}\n{}".format(name, traceback.format_exc()), flush=True)


def env(**extra):
    e = dict(os.environ)
    # the server's runs folder is reached through a symlink (/booleanfs2 -> /mnt/booleanfs2), so the tests are too
    (T / "runs").mkdir(parents=True, exist_ok=True)
    if not (T / "runs-link").is_symlink():
        (T / "runs-link").symlink_to(T / "runs")
    e.update(HEGEMON_RUNS=str(T / "runs-link"), HEGEMON_SITE_ROOT=str(T / "site"), HEGEMON_SITE_URL="http://127.0.0.1:8765",
             HEGEMON_BOOLEANLAB_SCRIPT=str(HERE / "FAKE_jupyter_gse_processing.sh"),
             HEGEMON_GENOME_DIR=str(T / "fixtures" / "genome"), HEGEMON_OFFLINE="1", HOME=str(T / "home"),
             OPENAI_API_KEY="sk-local-test", HEGEMON_OPENAI_URL="http://127.0.0.1:8799/v1/responses",
             HEGEMON_ATTEMPTS="3")
    e.update(extra)
    return e


def cli(*args, **extra):
    p = subprocess.run([PY, str(AGENT / "hegemon.py")] + list(args), env=env(**extra), stdin=subprocess.DEVNULL,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=900)
    return p.returncode, p.stdout


def wait_port(port, seconds=15):
    end = time.time() + seconds
    while time.time() < end:
        with socket.socket() as s:
            if s.connect_ex(("127.0.0.1", port)) == 0:
                return
        time.sleep(0.2)
    raise RuntimeError("port {} did not open".format(port))


def start(cmd, port, user=None):
    if user:
        cmd = ["runuser", "-u", user, "--"] + cmd
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    wait_port(port)
    return proc


def stop(proc):
    proc.terminate()
    proc.wait(10)


def php_site():
    if not shutil.which("php") or not (HEGEMON_REPO / "explore.php").is_file():
        raise Skip("needs php and a clone of github.com/sahoo00/Hegemon at " + str(HEGEMON_REPO))
    site = T / "site"
    if site.exists():
        shutil.rmtree(str(site))
    shutil.copytree(str(HEGEMON_REPO), str(site), ignore=shutil.ignore_patterns(".git", "kim-2010-*"))
    # Test-only shim: Hegemon's public util.php calls non-static methods statically (U::getX), which PHP 8
    # (this sandbox) rejects and older PHP accepts. The lab server runs Hegemon's code unmodified.
    util = site / "util.php"
    util.write_text(util.read_text().replace("\n  function ", "\n  static function "))
    (site / "explore.conf").write_text("title=Hierachical Exploration of Gene Expression Microarrays Online\n")
    (site / "index.html").write_text("<html><body><h1>Datasets</h1>\n<table>\n<tr><td>existing</td></tr>\n</table>\n</body></html>\n")
    subprocess.run(["chmod", "-R", "a+rX", str(site)], check=True)
    return site


def last_build(prefix_glob):
    builds = sorted((T / "runs" / "builds").glob(prefix_glob))
    return builds[-1]


# ----------------------------------------------------------------------------- tests

_TRICKY_SOFT = (
    b"^DATABASE = GeoMiame\r\n!Database_name = Gene Expression Omnibus (GEO)\r\n"
    b"^SERIES = GSE999900\n!Series_title = Tricky series\n"
    b"!Series_summary = a value with a lone \r carriage return inside it\n"
    b"!Series_type = Expression profiling by array\n!Series_platform_id = GPL9\n"
    b"^PLATFORM = GPL9\n!Platform_title = \xc2\xb5-array Stra\xc3\x9fe\n#ID = probe \n#Symbol = sym\n"
    b"!platform_table_begin\nID\tSymbol\np1\tGAPDH  \np2\t\n\np3\tACTB\r\np4\tbad\xff\xfebytes\n"
    b"p5\tKRT16\xc2\xa0\n!platform_table_end\n"
    b"^SAMPLE = GSM1\n!Sample_title = s1\n!Sample_platform_id = GPL9\n"
    b"!Sample_characteristics_ch1 = tissue: skin\n#ID_REF = \n#VALUE = v\n!sample_table_begin\n"
    b"ID_REF\tVALUE\np1\t1.5\np2\t2.5\n   \np3\t3.5\t\t\n"
    b"p4\tnull\np5\t7\n!sample_table_end\n"
    b"^SAMPLE = GSM2\n!Sample_title = s2\n!Sample_platform_id = GPL9\n#ID_REF = \n#VALUE = v\n"
    b"!sample_table_begin\nID_REF\tVALUE\np1\t1\t\np2\t2\t\np3\t3\rp4\t4\np5\t5\t\n!sample_table_end\n"
    b"^SAMPLE = GSM3\n!Sample_title = s3\n!Sample_platform_id = GPL9\n#ID_REF = \n#VALUE = v\n"
    b"!sample_table_begin\nID_REF\tVALUE\np1\t10\np2\t20\np3\t30\np4\t40\np5\t50\n!sample_table_end"
)


def t_soft_reader():
    try:
        import GEOparse
        import logging
        logging.disable(logging.CRITICAL)
    except ImportError:
        raise Skip("GEOparse not installed")
    import gzip
    from hegemon_nb import geo_soft
    from hegemon_nb.geo_soft import read_geo_soft
    tricky = T / "tricky_family.soft.gz"
    tricky.write_bytes(gzip.compress(_TRICKY_SOFT))
    files = [T / "fixtures" / "geo" / "GSE999001_family.soft.gz", tricky] + \
        sorted(Path("/tmp/geoparse_repo/tests").glob("*_family.soft"))
    real_block = geo_soft._BLOCK
    try:
        # tiny blocks put every line, table start and table end on a block boundary somewhere
        for block in (real_block, 5, 17, 64):
            geo_soft._BLOCK = block
            for f in files:
                ref = GEOparse.GEOparse.parse_GSE(str(f), open_kwargs={"encoding": "utf-8"})
                mine = read_geo_soft(f)
                assert list(ref.gsms) == list(mine.gsms) and list(ref.gpls) == list(mine.gpls), (f, block)
                assert ref.metadata == mine.metadata, (f, block, ref.metadata, mine.metadata)
                for kind in ("gsms", "gpls"):
                    for n, a in getattr(ref, kind).items():
                        b = getattr(mine, kind)[n]
                        assert a.metadata == b.metadata, (f, block, n, a.metadata, b.metadata)
                        assert a.columns.equals(b.columns), (f, block, n, a.columns, b.columns)
                        assert a.table.equals(b.table) and list(a.table.dtypes) == list(b.table.dtypes), \
                            (f, block, n, a.table, b.table)
                # the quick reads give the same metadata; sample_tables=N keeps only the first N tables
                light, first1 = read_geo_soft(f, tables=False), read_geo_soft(f, sample_tables=1)
                assert light.metadata == mine.metadata and first1.metadata == mine.metadata
                per_platform = {}
                for n, g in mine.gsms.items():
                    assert light.gsms[n].metadata == g.metadata and light.gsms[n].table.empty
                    plat = (g.metadata.get("platform_id") or ["?"])[0]
                    per_platform[plat] = per_platform.get(plat, 0) + 1
                    if per_platform[plat] == 1:
                        assert first1.gsms[n].table.equals(g.table), (f, n)
                    else:
                        assert first1.gsms[n].table.empty, (f, n)
                for n, g in mine.gpls.items():
                    assert first1.gpls[n].table.equals(g.table)
    finally:
        geo_soft._BLOCK = real_block
    assert len(read_geo_soft(tricky).gsms["GSM2"].table) == 5   # lone \r split a row, as GEOparse does


def t_helpers():
    import numpy as np
    import pandas as pd
    os.environ["HEGEMON_GENOME_DIR"] = str(T / "fixtures" / "genome")
    import importlib
    import hegemon_nb.helpers as h
    importlib.reload(h)
    gse = h.read_geo_soft(T / "fixtures" / "geo" / "GSE999001_family.soft.gz")
    surv, ih = h.geo_survival_ih(gse, h.gsm_ids(gse, "GPL999"))
    row = surv.set_index("ArrayId").loc["GSM900003"]
    assert row["c title"] == "Tumor: relapse 3", row["c title"]
    assert row["c treatment_ch1"] == "drug: 5 mg"
    assert surv.set_index("ArrayId").loc["GSM900005", "c age_ch1"] == ""
    assert "c Tumor_ch1" not in surv.columns and list(surv.columns[:3]) == ["ArrayId", "time", "status"]
    assert list(ih["ClinicalHeader"])[2] == "Tumor: relapse 3"
    names = h.geo_platform_names(gse.gpls["GPL999"]).set_index("ID")
    assert names.loc["200003_s_at", "Symbol"] == "RPL5 /// SNORA66"
    nm = h.name_from_symbol(names["Symbol"], names["Definition"])
    assert nm.loc["200002_at"] == "ACTB:actin beta" and pd.isna(nm.loc["200004_at"])
    assert h.lab_gene_names(["ENSG00000111640.5", "ENSG00000075624", "ENSG00000999999"]) == ["GAPDH", "ACTB", None]
    counts = pd.DataFrame({"a": [10, 30, 60], "b": [0, 5, 5]})
    expect = np.log2(counts / counts.sum() * 1e6 + 1)
    assert np.allclose(h.cpm_log2(counts).to_numpy(), expect.to_numpy())


def t_guard():
    from hegemon_nb.convert import guard
    roots = ["/data/source", "/lab/genome"]
    for script in ("fixture_geo_convert.py", "fixture_ae_convert.py"):
        assert guard((HERE / script).read_text(), roots) == [], script
    ok = "import pandas as pd\ndf = pd.DataFrame()\ndf = df.replace(0, 1).rename(columns={})\nx=[1];x.remove(1)\n"
    assert guard(ok, roots) == []
    bad = {"import subprocess": "subprocess", "import os\nos.remove('a')": "os.remove", "import shutil": "shutil",
           "eval('1')": "eval", "open('/etc/passwd')": "absolute path", "from pathlib import Path\nPath('a').unlink()": "unlink",
           "open('../x')": "leaving the folder", "import urllib.request": "urllib", "def f(:\n": "syntax"}
    for code, word in bad.items():
        found = " ".join(guard(code, roots))
        assert word in found, (code, found)


def _good_trio(d, prefix="GSE1-GPL1", samples=("GSM1", "GSM2", "GSM3")):
    import pandas as pd
    from hegemon_nb.helpers import write_native
    expr = pd.DataFrame({"ProbeID": ["p1", "p2", "p3"], "Name": ["GAPDH:g", "ACTB:a", "TP53:t"]})
    for i, s in enumerate(samples):
        expr[s] = [5.5 + i, 7.25, 9.0 - i]
    n = len(samples)
    surv = pd.DataFrame({"ArrayId": list(samples), "time": "", "status": "",
                         "c tissue_ch1": ["a" if i < n // 2 else "b" for i in range(n)]})
    ih = pd.DataFrame({"ArrayID": list(samples), "ArrayHeader": list(samples),
                       "ClinicalHeader": ["t{}".format(i + 1) for i in range(n)]})
    d.mkdir(parents=True, exist_ok=True)
    write_native(prefix, expr, surv, ih, str(d))
    return expr, surv, ih


def t_checks():
    from hegemon_nb.checks import check_core
    from hegemon_nb.helpers import write_native
    from hegemon_nb import native
    d = T / "checks"
    if d.exists():
        shutil.rmtree(str(d))
    expr, surv, ih = _good_trio(d)
    res = check_core(d, "GSE1-GPL1", "GEO", ["GSM1", "GSM2", "GSM3"])
    assert res.ok, res.text()

    def expect(label, mutate_expr=None, mutate_surv=None, mutate_ih=None, word=""):
        e, s, i = expr.copy(), surv.copy(), ih.copy()
        e = mutate_expr(e) if mutate_expr else e
        s = mutate_surv(s) if mutate_surv else s
        i = mutate_ih(i) if mutate_ih else i
        write_native("GSE1-GPL1", e, s, i, str(d))
        r = check_core(d, "GSE1-GPL1", "GEO", ["GSM1", "GSM2", "GSM3"])
        assert not r.ok and word in r.text(), (label, r.text())

    expect("duplicate ProbeID", lambda e: e.assign(ProbeID=["p1", "p1", "p3"]), word="repeat")
    expect("raw counts", lambda e: e.assign(GSM1=[5000, 12, 7], GSM2=[300, 10, 1], GSM3=[80, 40, 2]), word="raw counts")
    expect("unlogged", lambda e: e.assign(GSM1=[5000.5, 12.1, 7.2]), word="not log scale")
    expect("text in values", lambda e: e.assign(GSM2=["nan", "7", "8"]), word="not numbers")
    expect("non-GSM columns", lambda e: e.rename(columns={"GSM3": "Sample3"}), None,
           lambda i: i.assign(ArrayID=["GSM1", "GSM2", "Sample3"], ArrayHeader=["GSM1", "GSM2", "Sample3"]), word="GSM")
    expect("missing ih row", None, None, lambda i: i.iloc[:2], word="no ih row")
    expect("blank ClinicalHeader", None, None, lambda i: i.assign(ClinicalHeader=["t1", "", "t3"]), word="ClinicalHeader")
    expect("unprefixed metadata", None, lambda s: s.rename(columns={"c tissue_ch1": "tissue"}), word="must start with 'c '")
    expect("no metadata", None, lambda s: s.assign(**{"c tissue_ch1": ""}), word="no sample metadata")
    write_native("GSE1-GPL1", expr, surv, ih, str(d))
    p = d / "GSE1-GPL1-expr.txt"
    p.write_text(p.read_text().replace("GAPDH:g", "GAPDH:g\textra"))
    r = check_core(d, "GSE1-GPL1", "GEO", ["GSM1", "GSM2", "GSM3"])
    assert not r.ok and "wrong number of tab-separated fields" in r.text(), r.text()

    write_native("GSE1-GPL1", expr, surv, ih, str(d))
    native.make_idx(d, "GSE1-GPL1")
    assert native.check_idx(d, "GSE1-GPL1") is None
    run_dir = d / "GSE1-GPL1"
    if run_dir.exists():
        shutil.rmtree(str(run_dir))
    run_dir.mkdir()
    for f in d.glob("GSE1-GPL1-*.txt"):
        shutil.copy(str(f), str(run_dir / f.name))
    subprocess.run(["bash", str(HERE / "FAKE_jupyter_gse_processing.sh"), "GSE1-GPL1"], cwd=str(run_dir), check=True,
                   stdout=subprocess.DEVNULL)
    assert native.check_support(run_dir, "GSE1-GPL1", 3).ok
    thr = run_dir / "GSE1-GPL1-thr.txt"
    good = thr.read_text()
    thr.write_text(good + "".join("other{}\t1\t0\t0.5\t1.5\n".format(k) for k in range(5)))
    r = native.check_support(run_dir, "GSE1-GPL1", 3)
    assert not r.ok and "not in this expr file" in r.text(), r.text()
    thr.write_text("junk line one\njunk\njunk\njunk\n")
    assert not native.check_support(run_dir, "GSE1-GPL1", 3).ok
    thr.write_text(good)
    bv = run_dir / "GSE1-GPL1-bv.txt"
    bv.write_text(bv.read_text().replace("\tGAPDH:g\t", "\tGAPDH:g\t1"))
    r = native.check_support(run_dir, "GSE1-GPL1", 3)
    assert not r.ok and "bit vectors" in r.text(), r.text()
    old = time.time() - 3600
    os.utime(str(thr), (old, old))
    assert "older than the expr" in native.check_support(run_dir, "GSE1-GPL1", 3).text()


def t_geo_build_publish():
    site = php_site()
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    try:
        rc, out = cli("build", "GSE999001", "--source", str(T / "fixtures" / "geo"),
                      "--script", str(HERE / "fixture_geo_convert.py"), "--yes")
    finally:
        stop(php)
    assert rc == 0, out
    for word in ("files pass the checks: 6 samples x 30 probes", "every pointer lands on its expr row",
                 "lists GSE999001-GPL999 with n = 6", "Hegemon returned GAPDH", "read GAPDH's thresholds", "live: "):
        assert word in out, (word, out)
    conf = (site / "explore.conf").read_text()
    assert "[GSE999001-GPL999]" in conf and "source= GSE999001" in conf
    index = (site / "index.html").read_text()
    assert 'href="explore.php?key=gse999001gpl999&amp;A=ACTB&amp;B=GAPDH&amp;cmd=explore"' in index, index
    assert "the page opens on ACTB vs GAPDH with the plot and the Select Patient Information dropdown" in out, out
    receipt = json.loads((last_build("GSE999001-*") / "build.json").read_text())
    assert receipt["status"] == "built" and receipt["published"]["key"] == "gse999001gpl999"


def t_geo_vs_notebook():
    lab = T / "lab" / "GSE999001-GPL999"
    if not lab.is_dir():
        raise Skip("run run_notebook_reference.py geo first")
    agent = last_build("GSE999001-*") / "GSE999001-GPL999"
    for kind in ("expr", "idx", "ih"):
        a, b = (agent / "GSE999001-GPL999-{}.txt".format(kind)).read_bytes(), (lab / "GSE999001-GPL999-{}.txt".format(kind)).read_bytes()
        assert a == b, "{} differs from the notebook's".format(kind)
    from hegemon_nb.compare import compare
    same, report = compare(agent, lab)
    print("      " + report.replace("\n", "\n      "))
    assert "only in lab:   c Tumor_ch1, c treatmen_ch1" in report and "only in agent: c treatment_ch1" in report


def t_publish_rollback():
    site = php_site()
    if not shutil.which("runuser"):
        raise Skip("needs runuser")
    rc, out = cli("build", "GSE999001", "--source", str(T / "fixtures" / "geo"),
                  "--script", str(HERE / "fixture_geo_convert.py"), "--no-publish")
    assert rc == 0, out
    work = last_build("GSE999001-*")
    conf_before, index_before = (site / "explore.conf").read_bytes(), (site / "index.html").read_bytes()
    os.chmod(str(work), 0o700)  # the web server (nobody) can no longer reach the files
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765, user="nobody")
    try:
        rc, out = cli("publish", str(work))
    finally:
        stop(php)
        os.chmod(str(work), 0o755)
    assert rc == 1 and ("does not list the dataset" in out or "explore.php crashed" in out) and "were restored exactly" in out, out
    assert (site / "explore.conf").read_bytes() == conf_before and (site / "index.html").read_bytes() == index_before


def _mock(plan):
    (T / "plan.json").write_text(json.dumps(plan))
    (T / "ai-log.jsonl").write_text("")
    return start([PY, str(HERE / "mock_openai.py"), str(T / "plan.json"), str(T / "ai-log.jsonl"), "8799"], 8799)


def _reply(script):
    return "Source: GSM VALUE tables.\nSamples: GSM accessions.\nNormalization: none (RMA log2).\n```python\n{}\n```".format(script)


def t_ai_repair_loop():
    good = (HERE / "fixture_geo_convert.py").read_text()
    broken = good.replace('gse.gpls[platform]', 'gse.gpls["GPL000"]')
    ai = _mock([{"text": _reply(broken)}, {"text": _reply(good)}])
    try:
        rc, out = cli("build", "GSE999001", "--source", str(T / "fixtures" / "geo"), "--no-publish",
                      HEGEMON_BUILTIN="0")
    finally:
        stop(ai)
    assert rc == 0, out
    calls = [json.loads(l) for l in (T / "ai-log.jsonl").read_text().splitlines()]
    assert len(calls) == 2, len(calls)
    first, second = calls[0]["body"], calls[1]["body"]
    assert first["reasoning"]["effort"] == "medium" and second["reasoning"]["effort"] == "high", \
        (first.get("reasoning"), second.get("reasoning"))
    assert "BLUEPRINT: GSE51984-jupyter-process.txt" in first["instructions"]
    assert "BLUEPRINT: E-MTAB-7604-make-Hegemon-files.txt" in first["instructions"]
    assert "GSM900003 -> Tumor: relapse 3" in first["input"]
    assert "WHAT WENT WRONG" in second["input"] and "KeyError: 'GPL000'" in second["input"], second["input"][-800:]
    assert calls[0]["auth"].startswith("Bearer sk-")
    receipt = json.loads((last_build("GSE999001-*") / "build.json").read_text())
    assert receipt["attempts"] == 2 and receipt["ai_calls"] == 2


def t_ai_outcomes():
    cases = [
        ([{"status": 401, "error": {"message": "Incorrect API key", "code": "invalid_api_key"}}], "AI_PROVIDER_AUTH"),
        ([{"status": 429, "error": {"message": "You exceeded your quota", "code": "insufficient_quota"}}], "AI_PROVIDER_QUOTA"),
        ([{"status": 404, "error": {"message": "The model does not exist", "code": "model_not_found"}}], "AI_PROVIDER_MODEL"),
        ([{"text": "UNCONVERTIBLE: NO_SAMPLE_MATCH: columns are S1..S6 with no link to the GSMs"}], "NO_SAMPLE_MATCH"),
    ]
    for plan, word in cases:
        ai = _mock(plan)
        try:
            rc, out = cli("build", "GSE999001", "--source", str(T / "fixtures" / "geo"), "--no-publish",
                          HEGEMON_BUILTIN="0")
        finally:
            stop(ai)
        assert rc == 2 and word in out, (word, out[-1500:])
        if word.startswith("AI_PROVIDER"):
            assert "not the dataset" in out


def t_ae_build_vs_notebook():
    rc, out = cli("build", "E-MTAB-9999", "--source", str(T / "fixtures" / "ae"),
                  "--script", str(HERE / "fixture_ae_convert.py"), "--no-publish")
    assert rc == 0, out
    lab = T / "lab" / "E-MTAB-9999"
    if not lab.is_dir():
        raise Skip("built OK; run run_notebook_reference.py ae for the comparison")
    import numpy as np
    import pandas as pd
    from hegemon_nb.compare import compare
    agent_dir = last_build("E-MTAB-9999-*") / "E-MTAB-9999"
    same, report = compare(agent_dir, lab)
    print("      " + report.replace("\n", "\n      "))
    a = pd.read_csv(str(agent_dir / "E-MTAB-9999-expr.txt"), sep="\t").set_index("ProbeID")
    b = pd.read_csv(str(lab / "E-MTAB-9999-expr.txt"), sep="\t").set_index("ProbeID").loc[a.index]
    for s in ("S01", "S02", "S04", "S05"):  # one count file each: identical to the notebook
        assert (a[s] - b[s]).abs().max() < 1e-4, s
    # S03 is split across S03.L1/S03.L2. The notebook's `df[arr] += df1[arr]` adds them by row position after
    # its outer merge re-sorted the rows, i.e. it adds counts of different genes. The agent sums by gene.
    read = lambda f: pd.read_csv(str(f), sep="\t", header=0, names=["g", "c"]).iloc[:-5].set_index("g")["c"]
    s03 = read(T / "fixtures/ae/counts/S03.L1.count") + read(T / "fixtures/ae/counts/S03.L2.count")
    correct = np.log2(s03 / s03.sum() * 1e6 + 1)
    assert (a.loc[correct.index, "S03"] - correct).abs().max() < 1e-9
    assert (b.loc[correct.index, "S03"] - correct).abs().max() > 0.5, "expected the notebook's lane-sum bug"
    print("      S03 (split across two files): agent = correct by-gene sum; notebook adds different genes' counts")


def t_real_data_through_php():
    """Hegemon's own sample dataset (real Illumina data, 256 samples): notebook idx code + publish + PHP lookups."""
    if not (HEGEMON_REPO / "kim-2010-bladder-expr.txt").is_file():
        raise Skip("needs the public Hegemon repo")
    site = php_site()
    from hegemon_nb import native
    from hegemon_nb.config import Config
    from hegemon_nb import publish
    run_dir = T / "runs" / "kim" / "KIM2010"
    if run_dir.parent.exists():
        shutil.rmtree(str(run_dir.parent))
    run_dir.mkdir(parents=True)
    for src, kind in (("expr", "expr"), ("os", "survival"), ("ih", "ih"), ("thr", "thr")):
        shutil.copy(str(HEGEMON_REPO / "kim-2010-bladder-{}.txt".format(src)), str(run_dir / "KIM2010-{}.txt".format(kind)))
    n = native.make_idx(run_dir, "KIM2010")
    assert native.check_idx(run_dir, "KIM2010") is None
    import csv
    import pandas as pd
    mine = pd.read_csv(str(run_dir / "KIM2010-idx.txt"), sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)
    shipped = pd.read_csv(str(HEGEMON_REPO / "kim-2010-bladder-idx.txt"), sep="\t", dtype=str, keep_default_na=False,
                          quoting=csv.QUOTE_NONE)
    for col in ("ProbeID", "Ptr", "Name"):
        assert (mine[col] == shipped[col]).all(), col
    # The notebook splits Name:Description on ':' and keeps the space; Hegemon's shipped idx split on ': '.
    assert (mine["Description"].str.lstrip(" ") == shipped["Description"]).all()
    real_thr = (run_dir / "KIM2010-thr.txt").read_bytes()
    subprocess.run(["bash", str(HERE / "FAKE_jupyter_gse_processing.sh"), "KIM2010"], cwd=str(run_dir), check=True,
                   stdout=subprocess.DEVNULL)
    (run_dir / "KIM2010-thr.txt").write_bytes(real_thr)  # keep Hegemon's real thresholds
    samples = len(open(str(run_dir / "KIM2010-expr.txt")).readline().split("\t")) - 2
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    cfg = Config()
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    try:
        section, key, url = publish.register(cfg, {"run_dir": str(run_dir), "prefix": "KIM2010", "accession": "GSE13507",
                                                   "title": "Kim 2010 bladder", "samples": samples})
    finally:
        stop(php)
    assert section == "KIM2010" and n == 48702 and samples == 256


def t_repo_parsers():
    """Parsers for GEO pages (the downloads themselves need the real server)."""
    from hegemon_nb import repos
    from hegemon_nb.convert import extract
    assert repos.geo_series_url("GSE51984") == "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51984/"
    assert repos.geo_series_url("GSE123") == "https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE123/"
    assert repos.detect_repository("e-mtab-7604") == "ArrayExpress" and repos.detect_repository("GSE1") == "GEO"
    d = T / "pages"
    d.mkdir(parents=True, exist_ok=True)
    (d / "suppl.html").write_text(
        '<html><body><h1>Index of /geo/series/GSE51nnn/GSE51984/suppl</h1><pre>Name   Last modified   Size  <hr>'
        '<a href="/geo/series/GSE51nnn/GSE51984/">Parent Directory</a>          -   \n'
        '<a href="GSE51984_RAW.tar">GSE51984_RAW.tar</a>        2013-11-05 16:47   12M  \n'
        '<a href="filelist.txt">filelist.txt</a>            2013-11-05 16:47  1.2K  \n<hr></pre></body></html>')
    server = start([PY, "-m", "http.server", "8790", "--bind", "127.0.0.1", "--directory", str(d)], 8790)
    try:
        listing = repos._listing("http://127.0.0.1:8790/suppl.html")
    finally:
        stop(server)
    assert [(n, round(s or 0)) for n, _, s in listing] == [("GSE51984_RAW.tar", 12000000), ("filelist.txt", 1200)], listing
    filelist = ("#Archive/File\tName\tTime\tSize\tType\nArchive\tGSE51984_RAW.tar\t11/05/2013\t12615680\tTAR\n"
                "File\tGSM1256809_0_01_RPKM.txt.gz\t11/05/2013\t1011123\tTXT\n"
                "File\tGSM1256810_0_02_RPKM.txt.gz\t11/05/2013\t1010412\tTXT\n")
    assert dict(repos._raw_tar_contents(filelist)) == {"TXT": 2}
    code, notes, unc = extract("Source: x\nNorm: none\n```python\nprint(1)\n```")
    assert code.strip() == "print(1)" and notes == ["Source: x", "Norm: none"] and unc is None
    assert extract("UNCONVERTIBLE: NO_PROCESSED_DATA: only FASTQ in SRA")[2] == ("NO_PROCESSED_DATA", "only FASTQ in SRA")


def t_shared_key():
    """Several datasets under one key: Hegemon lists them all at one URL, with one homepage link."""
    site = php_site()
    from hegemon_nb import native, net, publish
    from hegemon_nb.config import Config
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    cfg = Config()
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    try:
        receipts = []
        for name in ("GSE7001-GPL1", "GSE7002-GPL1", "GSE7003-GPL1"):
            run_dir = T / "shared" / name
            if run_dir.exists():
                shutil.rmtree(str(run_dir))
            _good_trio(run_dir, name, tuple("GSM{}{}".format(name[6], k) for k in range(1, 5)))
            native.make_idx(run_dir, name)
            subprocess.run(["bash", str(HERE / "FAKE_jupyter_gse_processing.sh"), name], cwd=str(run_dir),
                           check=True, stdout=subprocess.DEVNULL)
            native.open_permissions(run_dir, T / "shared")
            receipt = {"run_dir": str(run_dir), "prefix": name, "accession": name.split("-")[0],
                       "title": name + " study", "samples": 4}
            section, key, url = publish.register(cfg, receipt, key="tregcd", link_label="Treg in Crohn's (batch)")
            assert key == "tregcd", key
            receipts.append((section, url))
        page = net.get_text(cfg.site_url + "/explore.php?key=tregcd", timeout=30, retries=0)
        listed = re.findall(r'<option value="([^"]+)"', page)
        assert set(listed) == {"GSE7001-GPL1", "GSE7002-GPL1", "GSE7003-GPL1"}, listed
        assert page.count("(n = 4)") == 3, page.count("(n = 4)")
        index = (site / "index.html").read_text()
        assert index.count("explore.php?key=tregcd") == 1, "expected exactly one homepage link per key"
        assert "Treg in Crohn's (batch)" in index
        assert len({u.split("&id=")[0] for _, u in receipts}) == 1, "all datasets should share one page"
        assert all("&id=" + s in u for s, u in receipts), receipts
        assert 'href="explore.php?key=tregcd&amp;A=ACTB&amp;B=GAPDH&amp;cmd=explore"' in index, index
        # the same accession cannot be added twice under one key
        try:
            publish.register(cfg, {"run_dir": str(T / "shared" / "GSE7001-GPL1"), "prefix": "GSE7001-GPL1",
                                   "accession": "GSE7001", "title": "dup", "samples": 4}, key="tregcd")
            raise AssertionError("expected a duplicate to be refused")
        except publish.PublishError as err:
            assert "already published under key" in str(err), err
        # a second key is independent
        run_dir = T / "shared" / "GSE7003-GPL1"
        publish.register(cfg, {"run_dir": str(run_dir), "prefix": "GSE7003-GPL1", "accession": "GSE7003",
                               "title": "other", "samples": 4}, key="other")
        page = net.get_text(cfg.site_url + "/explore.php?key=other", timeout=30, retries=0)
        assert re.findall(r'<option value="([^"]+)"', page) == ["GSE7003-GPL1-2"], page
        assert (site / "index.html").read_text().count("explore.php?key=tregcd") == 1
    finally:
        stop(php)


def t_metadata_plot_check():
    """A dataset whose metadata Hegemon cannot join must be refused, not published."""
    site = php_site()
    from hegemon_nb import native, publish
    from hegemon_nb.config import Config
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    cfg = Config()
    run_dir = T / "nometa" / "GSE8001-GPL1"
    if run_dir.parent.exists():
        shutil.rmtree(str(run_dir.parent))
    _good_trio(run_dir, "GSE8001-GPL1", ("GSM1", "GSM2", "GSM3", "GSM4"))
    # Exactly the professor's failure: survival/ih name samples the expr header does not,
    # so Hegemon shows the dataset but every metadata value comes back blank.
    for kind in ("survival", "ih"):
        f = run_dir / "GSE8001-GPL1-{}.txt".format(kind)
        f.write_text(f.read_text().replace("GSM", "SAMPLE_"))
    native.make_idx(run_dir, "GSE8001-GPL1")
    subprocess.run(["bash", str(HERE / "FAKE_jupyter_gse_processing.sh"), "GSE8001-GPL1"], cwd=str(run_dir),
                   check=True, stdout=subprocess.DEVNULL)
    native.open_permissions(run_dir, T / "nometa")
    conf_before, index_before = (site / "explore.conf").read_bytes(), (site / "index.html").read_bytes()
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    try:
        try:
            publish.register(cfg, {"run_dir": str(run_dir), "prefix": "GSE8001-GPL1", "accession": "GSE8001",
                                   "title": "unplottable", "samples": 4})
            raise AssertionError("expected the metadata check to refuse this dataset")
        except publish.PublishError as err:
            assert "every sample's value came back blank" in str(err), err
            assert "were restored exactly" in str(err), err
    finally:
        stop(php)
    assert (site / "explore.conf").read_bytes() == conf_before
    assert (site / "index.html").read_bytes() == index_before
    # and a dataset with no 'c '/'n ' field at all is refused too
    run2 = T / "nometa" / "GSE8002-GPL1"
    _good_trio(run2, "GSE8002-GPL1", ("GSM1", "GSM2", "GSM3", "GSM4"))
    f = run2 / "GSE8002-GPL1-survival.txt"
    f.write_text(f.read_text().replace("c tissue_ch1", "tissue"))
    native.make_idx(run2, "GSE8002-GPL1")
    subprocess.run(["bash", str(HERE / "FAKE_jupyter_gse_processing.sh"), "GSE8002-GPL1"], cwd=str(run2),
                   check=True, stdout=subprocess.DEVNULL)
    native.open_permissions(run2, T / "nometa")
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    try:
        try:
            publish.register(cfg, {"run_dir": str(run2), "prefix": "GSE8002-GPL1", "accession": "GSE8002",
                                   "title": "no c field", "samples": 4})
            raise AssertionError("expected refusal for a survival file with no 'c ' field")
        except publish.PublishError as err:
            assert "cannot be grouped or plotted" in str(err), err
    finally:
        stop(php)
    # and the case from the lab's review: metadata present, but only titles (all different) and
    # organism (all the same). The page would open with nothing to group by, so it must be refused.
    run3 = T / "nometa" / "GSE8003-GPL1"
    _good_trio(run3, "GSE8003-GPL1", ("GSM1", "GSM2", "GSM3", "GSM4"))
    f = run3 / "GSE8003-GPL1-survival.txt"
    f.write_text("ArrayId\ttime\tstatus\tc title\tc organism_ch1\n" + "".join(
        "GSM{}\t\t\tsample {}\tHomo sapiens\n".format(k, k) for k in range(1, 5)))
    native.make_idx(run3, "GSE8003-GPL1")
    subprocess.run(["bash", str(HERE / "FAKE_jupyter_gse_processing.sh"), "GSE8003-GPL1"], cwd=str(run3),
                   check=True, stdout=subprocess.DEVNULL)
    native.open_permissions(run3, T / "nometa")
    conf_before = (site / "explore.conf").read_bytes()
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    try:
        try:
            publish.register(cfg, {"run_dir": str(run3), "prefix": "GSE8003-GPL1", "accession": "GSE8003",
                                   "title": "titles only", "samples": 4})
            raise AssertionError("expected refusal: title and organism group nothing")
        except publish.PublishError as err:
            assert "no field splits the samples into groups" in str(err), err
            assert "c title (identifier/technical field)" in str(err), err
            assert "were restored exactly" in str(err), err
    finally:
        stop(php)
    assert (site / "explore.conf").read_bytes() == conf_before


def t_batch_shared_key():
    """A batch publishes what works under one key, survives failures, and continues on a second run."""
    site = php_site()
    good = (HERE / "fixture_geo_convert.py").read_text()
    # GSE999002 has no source files -> its build fails; the batch must carry on.
    plan = [{"text": "Rank\n" + json.dumps({"ranked": [{"accession": "GSE999001", "why": "direct match"},
                                                       {"accession": "GSE999002", "why": "second"},
                                                       {"accession": "GSE999003", "why": "third"}]})}]
    reply = "Source: GSM VALUE tables.\n```python\n" + good + "\n```"
    plan += [{"text": reply}] * 8
    ai = _mock(plan)
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    src = T / "fixtures" / "geo"
    fake_search = (
        "import sys; sys.path.insert(0, {agent!r})\n"
        "from hegemon_nb import repos\n"
        "def fake(query, limit=40):\n"
        "    return [repos.Dataset(a, 'GEO', a + ' study', kind='Expression profiling by array')\n"
        "            for a in ('GSE999001', 'GSE999002', 'GSE999003')]\n"
        "repos.search_geo = fake\n"
        "repos.search_arrayexpress = lambda q, limit=40: []\n"
        "_real = repos.fetch_geo\n"
        "def fetch(ds, source, cfg, yes=False, **kw):\n"
        "    import shutil, pathlib\n"
        "    if ds.accession == 'GSE999002':\n"
        "        raise repos.net.NetError('simulated download failure')\n"
        "    source = pathlib.Path(source); source.mkdir(parents=True, exist_ok=True)\n"
        "    shutil.copy({src!r} + '/GSE999001_family.soft.gz', str(source / (ds.accession + '_family.soft.gz')))\n"
        "    gse = repos.read_geo_soft(source / (ds.accession + '_family.soft.gz'), tables=False)\n"
        "    return source / (ds.accession + '_family.soft.gz'), gse, [], ['soft tables']\n"
        "repos.fetch_geo = fetch\n"
        "import hegemon\n"
        "sys.exit(hegemon.main(sys.argv[1:]))\n"
    ).format(agent=str(AGENT), src=str(src))
    runner = T / "batch_runner.py"
    runner.write_text(fake_search)
    try:
        p = subprocess.run([PY, str(runner), "batch", "regulatory T cells", "--key", "tregcd", "--count", "2",
                            "--min-samples", "0"],
                           env=env(), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True, timeout=900)
        out = p.stdout
    finally:
        stop(ai)
        stop(php)
    assert p.returncode == 0, out
    assert "simulated download failure" in out, out
    for word in ("1 of 2 published under 'tregcd'", "2 of 2 published under 'tregcd'",
                 "2 dataset(s) are now on one page under key 'tregcd'"):
        assert word in out, (word, out[-3000:])
    conf = (site / "explore.conf").read_text()
    assert conf.count("key= tregcd") == 2, conf
    assert {v["source"] for v in __import__("hegemon_nb.publish", fromlist=["x"])._sections(conf).values()
            if v.get("key") == "tregcd"} == {"GSE999001", "GSE999003"}
    assert (site / "index.html").read_text().count("explore.php?key=tregcd") == 1
    state = json.loads((T / "runs" / "batches" / "tregcd.json").read_text())
    assert state["query"] == "regulatory T cells" and state["tried"]["GSE999002"]["status"] == "DOWNLOAD_FAILED"
    # second run with the same key: skips what is published, reuses the saved query
    ai = _mock([{"text": "Rank\n" + json.dumps({"ranked": [{"accession": "GSE999002", "why": "retry"}]})}] +
               [{"text": reply}] * 4)
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    try:
        p2 = subprocess.run([PY, str(runner), "batch", "--key", "tregcd", "--count", "1", "--retry-failed",
                             "--min-samples", "0"],
                            env=env(), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True, timeout=900)
    finally:
        stop(ai)
        stop(php)
    assert "2 dataset(s) already published under 'tregcd': GSE999001, GSE999003" in p2.stdout, p2.stdout[-2500:]


def t_arrayexpress_sdrf():
    """The SDRF is mandatory, its per-sample files must all arrive, and its columns are parsed."""
    from hegemon_nb import net, repos
    from hegemon_nb.config import Config
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    cfg = Config()
    # the URL pattern used when BioStudies' info endpoint is unavailable (no network in tests)
    real_json = repos.net.get_json
    def offline(*a, **k):
        raise net.NetError("offline in tests")
    repos.net.get_json = offline
    try:
        assert repos.ae_files_base("E-MTAB-7604") == \
            "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/604/E-MTAB-7604/Files/", repos.ae_files_base("E-MTAB-7604")
        assert repos.ae_files_base("E-MTAB-13084").endswith("/E-MTAB-/084/E-MTAB-13084/Files/")
    finally:
        repos.net.get_json = real_json
    d = T / "ae_serve"
    if d.exists():
        shutil.rmtree(str(d))
    (d / "Files").mkdir(parents=True)
    (d / "Files" / "E-MTAB-9001.sdrf.txt").write_text(
        "Source Name\tCharacteristics[disease]\tComment[ENA_RUN]\tArray Data File\tDerived Array Data File\n"
        "S1\tcontrol\tERR1\tS1.fastq.gz\tS1.count\n"
        "S2\tcrohn\tERR2\tS2.fastq.gz\tS2.count\n")
    for name in ("S1.count", "S2.count"):
        (d / "Files" / name).write_text("GAPDH\t100\nACTB\t200\n")
    (d / "Files" / "S1.fastq.gz").write_bytes(b"raw")
    # counts live in a subfolder, as in E-MTAB-7604: a flat listing cannot see them,
    # so the SDRF's names must fetch them and the folder structure must be kept.
    (d / "Files" / "counts").mkdir()
    (d / "Files" / "counts" / "S3.count").write_text("GAPDH\t300\nACTB\t400\n")
    sdrf = d / "Files" / "E-MTAB-9001.sdrf.txt"
    sdrf.write_text(sdrf.read_text().replace(
        "S2\tcrohn\tERR2\tS2.fastq.gz\tS2.count\n",
        "S2\tcrohn\tERR2\tS2.fastq.gz\tS2.count\nS3\tcrohn\tERR3\tS3.fastq.gz\tcounts/S3.count\n"))
    server = start([PY, "-m", "http.server", "8791", "--bind", "127.0.0.1", "--directory", str(d)], 8791)
    try:
        named = repos.sdrf_data_files(d / "Files" / "E-MTAB-9001.sdrf.txt")
        assert named == ["S1.count", "S2.count", "counts/S3.count"], named  # raw column excluded
        real_base, real_api = repos.ae_files_base, repos._biostudies_files
        real_download = repos.net.download

        def local_only(url, *a, **k):   # the real EBI fallback is unreachable here; fail fast instead of retrying
            if not url.startswith("http://127.0.0.1"):
                raise net.NetError("offline in tests: " + url)
            return real_download(url, *a, **k)
        repos.net.download = local_only
        repos.ae_files_base = lambda acc: "http://127.0.0.1:8791/Files/"
        repos._biostudies_files = lambda acc: []   # no EBI API in tests; FTP listing + SDRF are used
        try:
            out = T / "runs-link" / "ae_out"
            if out.exists():
                shutil.rmtree(str(out))
            notes, got = repos.fetch_arrayexpress(repos.Dataset("E-MTAB-9001", "ArrayExpress", "t"), out, cfg)
            assert (out / "E-MTAB-9001.sdrf.txt").is_file()
            assert {"S1.count", "S2.count"} <= {p.name for p in out.iterdir()}, list(out.iterdir())
            assert (out / "counts" / "S3.count").is_file(), "subfolder structure must be kept"
            assert "S1.fastq.gz" not in {p.name for p in out.rglob("*")}
            assert any("Fetched the 3 processed file(s) the SDRF names; nothing else was downloaded" in n
                       for n in notes), notes
            # a path trying to escape SOURCE_DIR is refused
            assert repos._safe_under(out, "../escape.txt") is None
            assert repos._safe_under(out, "a/../../b") is None
            # an absolute path is treated as relative, so it still lands inside SOURCE_DIR
            assert str(repos._safe_under(out, "/etc/passwd")).startswith(str(out.resolve()) + os.sep)
            assert repos._safe_under(out, "counts/ok.txt") is not None
            # remove one named file, and put a big file nobody named next to it (E-MTAB-184's case):
            # it must say the named file is not served, without downloading the big one first
            (d / "Files" / "S2.count").unlink()
            (d / "Files" / "nonorm_nobkgd_DataFile.txt").write_bytes(b"x" * (2 * 1024 * 1024))
            shutil.rmtree(str(out))
            # the build does not stop there: the problem goes to the AI with the study's other files offered
            ds_ae = repos.Dataset("E-MTAB-9001", "ArrayExpress", "t")
            notes, got = repos.fetch_arrayexpress(ds_ae, out, cfg)
            assert any("S2.count" in f and "does not serve it" in f for f in ds_ae.failures), ds_ae.failures
            assert "nonorm_nobkgd_DataFile.txt" in ds_ae.catalog, list(ds_ae.catalog)
            assert not (out / "nonorm_nobkgd_DataFile.txt").exists(), "downloaded a file the SDRF does not name"
            (d / "Files" / "nonorm_nobkgd_DataFile.txt").unlink()
            # and with no SDRF at all it says so instead of blaming the dataset
            (d / "Files" / "E-MTAB-9001.sdrf.txt").unlink()
            shutil.rmtree(str(out), ignore_errors=True)
            try:
                repos.fetch_arrayexpress(repos.Dataset("E-MTAB-9001", "ArrayExpress", "t"), out, cfg)
                raise AssertionError("expected a refusal with no SDRF")
            except net.NetError as err:
                assert "sdrf.txt" in str(err) and "sample metadata" in str(err), err
        finally:
            repos.ae_files_base, repos._biostudies_files = real_base, real_api
            repos.net.download = real_download
    finally:
        stop(server)


def t_diagnose():
    """diagnose must name the broken dataset on a site that has a good one and a bad one."""
    site = php_site()
    from hegemon_nb import native, publish
    from hegemon_nb.config import Config
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    cfg = Config()
    base = T / "diag"
    if base.exists():
        shutil.rmtree(str(base))
    made = {}
    for name, break_it in (("GSE9001-GPL1", False), ("GSE9002-GPL1", True)):
        run_dir = base / name
        _good_trio(run_dir, name, ("GSM1", "GSM2", "GSM3", "GSM4"))
        if break_it:  # the professor's site: metadata Hegemon cannot join
            for kind in ("survival", "ih"):
                f = run_dir / "{}-{}.txt".format(name, kind)
                f.write_text(f.read_text().replace("GSM", "SAMPLE_"))
        native.make_idx(run_dir, name)
        subprocess.run(["bash", str(HERE / "FAKE_jupyter_gse_processing.sh"), name], cwd=str(run_dir),
                       check=True, stdout=subprocess.DEVNULL)
        native.open_permissions(run_dir, base)
        made[name] = run_dir
        # register by hand: the broken one would never pass the publish checks
        conf = site / "explore.conf"
        conf.write_text(conf.read_text().rstrip("\n") + "\n\n[{0}]\nname= {0} study\nexpr= {1}/{0}-expr.txt\n"
                        "index= {1}/{0}-idx.txt\nsurvival= {1}/{0}-survival.txt\nindexHeader= {1}/{0}-ih.txt\n"
                        "info= {1}/{0}-info.txt\nkey= tregcd\nsource= {2}\n".format(
                            name, run_dir, name.split("-")[0]))
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    try:
        rows = publish.diagnose(cfg, "tregcd")
    finally:
        stop(php)
    by = {r["section"]: r for r in rows}
    assert set(by) == {"GSE9001-GPL1", "GSE9002-GPL1"}, list(by)
    assert not by["GSE9001-GPL1"]["problems"], by["GSE9001-GPL1"]["problems"]
    assert by["GSE9001-GPL1"]["gene"] in ("GAPDH", "ACTB"), by["GSE9001-GPL1"]
    bad = " ".join(by["GSE9002-GPL1"]["problems"])
    assert "survival.txt's ArrayId" in bad, bad
    # and the CLI reports it with a nonzero exit
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    try:
        rc, out = cli("diagnose", "--key", "tregcd", "--min-samples", "0")
    finally:
        stop(php)
    assert rc == 1, out
    assert "1 of 2 dataset(s) are fine" in out and "hegemon unpublish GSE9002-GPL1" in out, out


def _fake_geo_series(root, acc, kind, samples, suppl=None, per_sample=None, title=None):
    """A minimal GEO series layout: soft/<acc>_family.soft.gz plus suppl/ files."""
    import gzip
    series = root / acc
    (series / "soft").mkdir(parents=True, exist_ok=True)
    (series / "suppl").mkdir(parents=True, exist_ok=True)
    lines = ["^SERIES = " + acc, "!Series_title = " + (title or acc + " test"), "!Series_type = " + kind,
             "!Series_platform_id = GPL999"]
    lines += ["!Series_sample_id = " + s for s in samples]
    lines += ["^PLATFORM = GPL999", "!Platform_title = test", "!platform_table_begin", "ID\tSymbol",
              "p1\tGAPDH", "!platform_table_end"]
    for i, s in enumerate(samples):
        lines += ["^SAMPLE = " + s, "!Sample_title = sample " + str(i + 1),
                  "!Sample_characteristics_ch1 = disease: crohn", "!Sample_platform_id = GPL999",
                  "!Sample_data_row_count = 0"]
        for j, name in enumerate(per_sample or [], start=1):
            lines.append("!Sample_supplementary_file_{} = http://127.0.0.1:8792/geo/samples/x/{}_{}".format(
                j, s, name))
    with gzip.open(str(series / "soft" / (acc + "_family.soft.gz")), "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    for name, size in (suppl or {}).items():
        target = series / "suppl" / name
        if name.endswith(".tar"):   # a real tar, so unpacking is exercised for real
            import io
            import tarfile
            payload = b"gene\tvalue\n" + b"GAPDH\t5\n" * 200
            with tarfile.open(str(target), "w") as tf:
                for k in range(max(1, size // len(payload) // 50) or 1):
                    info = tarfile.TarInfo("GSM{}_inside.txt".format(k))
                    info.size = len(payload)
                    tf.addfile(info, io.BytesIO(payload))
        elif name.endswith(".gz"):
            target.write_bytes(gzip.compress(b"x" * size))
        else:
            target.write_bytes(b"x" * size)
    for s in samples:
        for name in (per_sample or []):
            (root / "geo" / "samples" / "x").mkdir(parents=True, exist_ok=True)
            body = b"gene\t1\n" * 20
            (root / "geo" / "samples" / "x" / "{}_{}".format(s, name)).write_bytes(
                gzip.compress(body) if name.endswith(".gz") else body)
    return series


def t_geo_acquisition_plan():
    """It must pick the cheapest source that works, and never pull a bundle it does not need."""
    from hegemon_nb import repos
    from hegemon_nb.config import Config
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    cfg = Config()
    root = T / "geo_serve"
    if root.exists():
        shutil.rmtree(str(root))
    root.mkdir(parents=True)
    seq = "Expression profiling by high throughput sequencing"
    big = 3 * 1024 * 1024  # stands in for the 1.5 GB bundle
    _fake_geo_series(root, "GSE999100", seq, ["GSM1", "GSM2"], suppl={"GSE999100_RAW.tar": big})
    _fake_geo_series(root, "GSE999200", seq, ["GSM3", "GSM4"], suppl={"GSE999200_RAW.tar": big},
                     per_sample=["counts.txt.gz", "aligned.bam"])
    _fake_geo_series(root, "GSE999300", seq, ["GSM5", "GSM6"], suppl={"GSE999300_RAW.tar": big})
    (root / "counts").mkdir()
    import gzip
    (root / "counts" / "GSE999100_raw_counts.tsv.gz").write_bytes(gzip.compress(b"gene\tGSM1\tGSM2\n"))

    server = start([PY, "-m", "http.server", "8792", "--bind", "127.0.0.1", "--directory", str(root)], 8792)
    real_url, real_counts = repos.geo_series_url, repos._ncbi_counts_links
    repos.geo_series_url = lambda acc: "http://127.0.0.1:8792/{}/".format(acc)
    repos._ncbi_counts_links = lambda acc: (
        ["http://127.0.0.1:8792/counts/GSE999100_raw_counts.tsv.gz?acc={}&file=GSE999100_raw_counts.tsv.gz".format(acc)]
        if acc == "GSE999100" else [])
    try:
        # 1. GEO's own counts exist -> only those, never the bundle
        out = T / "runs-link" / "acq1"
        shutil.rmtree(str(out), ignore_errors=True)
        soft, gse, notes, processed = repos.fetch_geo(repos.Dataset("GSE999100", "GEO", "t"), out, cfg)
        names = {p.name for p in out.rglob("*") if p.is_file()}
        assert "GSE999100_raw_counts.tsv.gz" in names, names
        assert "GSE999100_RAW.tar" not in names, "downloaded the bundle when counts were available"
        assert any("were not needed" in n for n in notes), notes
        assert sum(p.stat().st_size for p in out.rglob("*") if p.is_file()) < big

        # 2. no counts, but per-sample files -> fetched individually into GSE###_RAW/, bundle skipped
        out = T / "runs-link" / "acq2"
        shutil.rmtree(str(out), ignore_errors=True)
        soft, gse, notes, processed = repos.fetch_geo(repos.Dataset("GSE999200", "GEO", "t"), out, cfg)
        names = {str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()}
        assert "GSE999200_RAW/GSM3_counts.txt.gz" in names, names
        assert "GSE999200_RAW/GSM4_counts.txt.gz" in names, names
        assert not any("RAW.tar" in n for n in names), "bundle downloaded even though per-sample files exist"
        assert not any(n.endswith(".bam") for n in names), "downloaded a raw alignment file"
        assert any("bundle was not downloaded" in n for n in notes), notes
        assert any("Skipped 2 raw per-sample file" in n for n in notes), notes

        # 3. nothing else available -> the bundle is the last resort, and it is used
        out = T / "runs-link" / "acq3"
        shutil.rmtree(str(out), ignore_errors=True)
        soft, gse, notes, processed = repos.fetch_geo(repos.Dataset("GSE999300", "GEO", "t"), out, cfg)
        names = {p.name for p in out.rglob("*") if p.is_file()}
        assert "GSE999300_RAW.tar" in names, names
        assert any("only processed source" in n for n in notes), notes
        assert [p for p in out.rglob("*") if p.name.endswith("_inside.txt")], "tar was not unpacked"
        # ftp:// URLs in SOFT files are served over https
        assert repos._https("ftp://ftp.ncbi.nlm.nih.gov/geo/x.gz") == "https://ftp.ncbi.nlm.nih.gov/geo/x.gz"

        # 4. a size limit refuses the bundle and says so, instead of looking like an empty dataset
        out = T / "runs-link" / "acq4"
        shutil.rmtree(str(out), ignore_errors=True)
        cfg.max_file_gb = 0.000001
        try:
            soft, gse, notes, processed = repos.fetch_geo(repos.Dataset("GSE999300", "GEO", "t"), out, cfg)
        finally:
            cfg.max_file_gb = 5
        assert not [p for p in out.rglob("*") if p.name.endswith(".tar")], "downloaded past the size limit"
        assert any("NOT downloaded because of size" in n for n in notes), notes
        assert not processed, processed   # so the build reports FILES_TOO_LARGE, not "no data"

        # 5. any dataset at all: given URLs, it downloads and converts them like the rest
        out = T / "runs-link" / "acq5"
        shutil.rmtree(str(out), ignore_errors=True)
        notes, got = repos.fetch_urls(repos.Dataset("MYDATA", "Other", "t"), out, cfg,
                                      ["http://127.0.0.1:8792/counts/GSE999100_raw_counts.tsv.gz"])
        assert got == ["GSE999100_raw_counts.tsv.gz"], got
        assert repos.detect_repository("PXD000001", allow_other=True) == "Other"
        try:
            repos.detect_repository("PXD000001")
            raise AssertionError("expected a clear message for an unknown accession")
        except ValueError as err:
            assert "--url" in str(err) and "--source" in str(err), err
    finally:
        repos.geo_series_url, repos._ncbi_counts_links = real_url, real_counts
        stop(server)


def _variant_soft(acc, drop_characteristics=False, same_source=False):
    import gzip
    src = (T / "fixtures" / "geo" / "GSE999001_family.soft.gz")
    text = gzip.open(str(src), "rt").read().replace("GSE999001", acc)
    lines = text.splitlines()
    if drop_characteristics:
        lines = [l for l in lines if not l.startswith("!Sample_characteristics_ch1")]
    if same_source:
        lines = [("!Sample_source_name_ch1 = human colon" if l.startswith("!Sample_source_name_ch1") else l)
                 for l in lines]
    d = T / "variants" / acc
    d.mkdir(parents=True, exist_ok=True)
    with gzip.open(str(d / (acc + "_family.soft.gz")), "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    return d


def t_dropped_metadata_is_repaired():
    """If the AI's script drops metadata the source has, the agent names the fields and gets them back."""
    good = (HERE / "fixture_geo_convert.py").read_text()
    dropping = good.replace(
        "survival, ih = geo_survival_ih(gse, list(data.columns[2:]))",
        "survival, ih = geo_survival_ih(gse, list(data.columns[2:]))\n"
        "survival = survival[['ArrayId', 'time', 'status', 'c title', 'c source_name_ch1']]")
    assert dropping != good
    ai = _mock([{"text": _reply(dropping)}, {"text": _reply(good)}])
    try:
        rc, out = cli("build", "GSE999001", "--source", str(T / "fixtures" / "geo"), "--no-publish",
                      HEGEMON_BUILTIN="0")
    finally:
        stop(ai)
    assert rc == 0, out
    assert "survival.txt dropped" in out and "c tissue_ch1" in out, out
    calls = [json.loads(l) for l in (T / "ai-log.jsonl").read_text().splitlines()]
    assert len(calls) == 2, len(calls)
    assert "missing metadata fields" in calls[1]["body"]["input"], calls[1]["body"]["input"][-1500:]
    assert "c tissue_ch1" in calls[1]["body"]["input"]
    assert "groups to compare by 'c tissue_ch1': tumor 3 / normal 3" in out, out


def t_no_grouping_metadata():
    """A study whose samples carry nothing that splits them is refused with that reason, after one look
    for a study-supplied sample table - never published as a page with nothing to plot by."""
    src = _variant_soft("GSE999004", drop_characteristics=True, same_source=True)
    good = (HERE / "fixture_geo_convert.py").read_text()
    ai = _mock([{"text": _reply(good)},
                {"text": "UNCONVERTIBLE: NO_GROUPING_METADATA: samples only have titles and one source name"}])
    try:
        rc, out = cli("build", "GSE999004", "--source", str(src), "--no-publish", HEGEMON_BUILTIN="0")
    finally:
        stop(ai)
    assert rc == 2 and "NO_GROUPING_METADATA" in out, out
    assert "asking the AI to look for a sample table" in out, out
    calls = [json.loads(l) for l in (T / "ai-log.jsonl").read_text().splitlines()]
    assert len(calls) == 2 and "study-supplied per-sample table" in calls[1]["body"]["input"]
    assert "c source_name_ch1: the same value for every sample" in calls[1]["body"]["input"], calls[1]["body"]["input"][-2000:]
    receipt = json.loads((last_build("GSE999004-*") / "build.json").read_text())
    assert receipt["status"] == "NO_GROUPING_METADATA"
    # and without an AI (a supplied script) it is refused straight away
    rc, out = cli("build", "GSE999004", "--source", str(src), "--script", str(HERE / "fixture_geo_convert.py"),
                  "--no-publish", OPENAI_API_KEY="")
    assert rc == 2 and "NO_GROUPING_METADATA" in out, out


def t_min_samples():
    """Small studies are left out before anything downloads, and caught after conversion if their size was unknown."""
    rc, out = cli("build", "GSE999001", "--source", str(T / "fixtures" / "geo"),
                  "--script", str(HERE / "fixture_geo_convert.py"), "--no-publish", "--min-samples", "10")
    assert rc == 2 and "TOO_FEW_SAMPLES" in out and "6 samples, fewer than --min-samples 10" in out, out
    # the batch leaves them out from GEO's own sample counts, without downloading anything
    runner = T / "batch_small.py"
    runner.write_text(
        "import sys; sys.path.insert(0, {agent!r})\n"
        "from hegemon_nb import repos\n"
        "repos.search_geo = lambda q, limit=40: [repos.Dataset(a, 'GEO', a, kind='Expression profiling by array', "
        "n_samples=6) for a in ('GSE1', 'GSE2', 'GSE3')]\n"
        "repos.search_arrayexpress = lambda q, limit=40: []\n"
        "def boom(*a, **k):\n"
        "    raise SystemExit('downloaded something it should have skipped')\n"
        "repos.fetch_geo = boom\n"
        "import hegemon\n"
        "sys.exit(hegemon.main(sys.argv[1:]))\n".format(agent=str(AGENT)))
    p = subprocess.run([PY, str(runner), "batch", "tregs", "--key", "small", "--min-samples", "10"], env=env(),
                       stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True, timeout=300)
    assert p.returncode == 1, p.stdout
    assert "left out 3 studies with fewer than 10 samples" in p.stdout, p.stdout
    assert "should have skipped" not in p.stdout
    # --max-samples leaves out the biggest the same way, and a build over it stops before converting
    p = subprocess.run([PY, str(runner), "batch", "tregs", "--key", "big", "--min-samples", "0", "--max-samples", "5"],
                       env=env(), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True, timeout=300)
    assert p.returncode == 1 and "left out 3 studies with more than 5 samples" in p.stdout, p.stdout
    # batches leave out studies over 250 samples unless told otherwise
    big_runner = T / "batch_big.py"
    big_runner.write_text(runner.read_text().replace("n_samples=6", "n_samples=300"))
    p = subprocess.run([PY, str(big_runner), "batch", "tregs", "--key", "big2"], env=env(), stdin=subprocess.DEVNULL,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=300)
    assert p.returncode == 1 and "left out 3 studies with more than 250 samples" in p.stdout, p.stdout
    assert "should have skipped" not in p.stdout
    rc, out = cli("build", "GSE999001", "--source", str(T / "fixtures" / "geo"), "--no-publish", "--max-samples", "5")
    assert rc == 2 and "TOO_MANY_SAMPLES" in out and "6 samples on GPL999, more than --max-samples 5" in out, out
    assert "Built-in conversion" not in out and "Conversion attempt" not in out, "converted before checking the size"
    # and GEO's SOFT file alone decides it: nothing but the SOFT is fetched
    from hegemon_nb import repos
    from hegemon_nb.config import Config
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    root = T / "geo_serve_small"
    shutil.rmtree(str(root), ignore_errors=True)
    root.mkdir(parents=True)
    _fake_geo_series(root, "GSE999500", "Expression profiling by high throughput sequencing", ["GSM1", "GSM2"],
                     suppl={"GSE999500_RAW.tar": 4096})
    server = start([PY, "-m", "http.server", "8793", "--bind", "127.0.0.1", "--directory", str(root)], 8793)
    real_url = repos.geo_series_url
    repos.geo_series_url = lambda acc: "http://127.0.0.1:8793/{}/".format(acc)
    try:
        out_dir = T / "small_out"
        shutil.rmtree(str(out_dir), ignore_errors=True)
        try:
            repos.fetch_geo(repos.Dataset("GSE999500", "GEO", "t"), out_dir, Config(), min_samples=10)
            raise AssertionError("expected TooFewSamples")
        except repos.TooFewSamples as err:
            assert err.count == 2 and err.minimum == 10
        assert [p.name for p in out_dir.iterdir()] == ["GSE999500_family.soft.gz"], list(out_dir.iterdir())
    finally:
        repos.geo_series_url = real_url
        stop(server)


def t_unpublish():
    """Take datasets off the page: by section or by key, with the dead homepage link removed and a backup kept."""
    t_shared_key()  # leaves GSE7001..3 under 'tregcd' and GSE7003 again under 'other'
    site = T / "site"
    from hegemon_nb import net, publish
    from hegemon_nb.config import Config
    cfg = Config()
    before_conf, before_index = (site / "explore.conf").read_text(), (site / "index.html").read_text()
    # asking without --yes and without a terminal changes nothing
    rc, out = cli("unpublish", "--key", "tregcd")
    assert rc == 1 and "Nothing changed" in out, out
    assert (site / "explore.conf").read_text() == before_conf
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    try:
        rc, out = cli("unpublish", "GSE7002-GPL1", "--yes")
        assert rc == 0 and "removed 1 dataset(s)" in out, out
        conf = (site / "explore.conf").read_text()
        assert "[GSE7002-GPL1]" not in conf and "[GSE7001-GPL1]" in conf
        assert (site / "index.html").read_text().count("explore.php?key=tregcd") == 1, "link must stay: others remain"
        page = net.get_text(cfg.site_url + "/explore.php?key=tregcd", timeout=30, retries=0)
        assert set(re.findall(r'<option value="([^"]+)"', page)) == {"GSE7001-GPL1", "GSE7003-GPL1"}
        rc, out = cli("unpublish", "--key", "tregcd", "--yes")
        assert rc == 0 and "removed 2 dataset(s)" in out and "homepage link(s) for: tregcd" in out, out
        conf = (site / "explore.conf").read_text()
        assert "key= tregcd" not in conf and "[GSE7003-GPL1-2]" in conf, conf   # the 'other' key survives
        index = (site / "index.html").read_text()
        assert "key=tregcd" not in index and "key=other" in index and "existing" in index, index
        page = net.get_text(cfg.site_url + "/explore.php?key=other", timeout=30, retries=0)
        assert re.findall(r'<option value="([^"]+)"', page) == ["GSE7003-GPL1-2"], page
        backups = sorted((site / ".hegemon-agent-backups").glob("unpublish-*"))
        assert backups and (backups[-1] / "explore.conf").is_file()
        rc, out = cli("unpublish", "--key", "nosuchkey", "--yes")
        assert rc == 1 and "Nothing in explore.conf matches" in out, out
    finally:
        stop(php)


def _ae_server(root, port):
    from hegemon_nb import net, repos
    server = start([PY, "-m", "http.server", str(port), "--bind", "127.0.0.1", "--directory", str(root)], port)
    saved = (repos.ae_files_base, repos._biostudies_files, repos.net.download)
    real_download = repos.net.download

    def local_only(url, *a, **k):
        if not url.startswith("http://127.0.0.1"):
            raise net.NetError("offline in tests: " + url)
        return real_download(url, *a, **k)
    repos.ae_files_base = lambda acc: "http://127.0.0.1:{}/Files/".format(port)
    repos._biostudies_files = lambda acc: []
    repos.net.download = local_only
    return server, saved


def _ae_restore(server, saved):
    from hegemon_nb import repos
    repos.ae_files_base, repos._biostudies_files, repos.net.download = saved
    stop(server)


def t_arrayexpress_archives_and_limits():
    """The SDRF-named file is found inside the study's processed archive; raw archives and files nobody
    named are never downloaded; a named file over the size limit is FILES_TOO_LARGE, not 'missing'."""
    import io
    import zipfile
    from hegemon_nb import repos
    from hegemon_nb.config import Config
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    cfg = Config()
    root = T / "ae_archive"
    shutil.rmtree(str(root), ignore_errors=True)
    (root / "Files").mkdir(parents=True)
    (root / "Files" / "E-MTAB-9002.sdrf.txt").write_text(
        "Source Name\tCharacteristics[disease]\tArray Data File\tDerived Array Data File\n"
        "A1\tcrohn\tA1.CEL\tfiltData.txt\nA2\tcontrol\tA2.CEL\tfiltData.txt\n")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("filtData.txt", "ID_REF\tA1\tA2\nGAPDH\t9.1\t8.7\n")
    (root / "Files" / "E-MTAB-9002.processed.1.zip").write_bytes(buf.getvalue())
    (root / "Files" / "E-MTAB-9002.raw.1.zip").write_bytes(b"raw" * 1000)
    (root / "Files" / "big_unrelated_DataFile.txt").write_bytes(b"x" * (2 * 1024 * 1024))
    server, saved = _ae_server(root, 8794)
    try:
        out = T / "runs-link" / "ae_archive_out"
        shutil.rmtree(str(out), ignore_errors=True)
        notes, got = repos.fetch_arrayexpress(repos.Dataset("E-MTAB-9002", "ArrayExpress", "t"), out, cfg)
        names = {p.name for p in out.rglob("*") if p.is_file()}
        assert "filtData.txt" in names, names
        assert "E-MTAB-9002.raw.1.zip" not in names, "downloaded the raw archive"
        assert "big_unrelated_DataFile.txt" not in names, "downloaded a file the SDRF does not name"
        assert any("nothing else was downloaded" in n for n in notes), notes

        # the named file over the size limit -> FilesTooLarge, not 'missing'
        (root / "Files" / "E-MTAB-9003.sdrf.txt").write_text(
            "Source Name\tCharacteristics[disease]\tDerived Array Data File\nB1\tcrohn\tbig.txt\n"
            "B2\tcontrol\tbig.txt\n")
        (root / "Files" / "big.txt").write_bytes(b"y" * (3 * 1024 * 1024))
        out = T / "runs-link" / "ae_big_out"
        shutil.rmtree(str(out), ignore_errors=True)
        cfg.max_file_gb = 0.000001
        try:
            repos.fetch_arrayexpress(repos.Dataset("E-MTAB-9003", "ArrayExpress", "t"), out, cfg)
            raise AssertionError("expected FilesTooLarge")
        except repos.FilesTooLarge as err:
            assert "big.txt" in str(err), err
        finally:
            cfg.max_file_gb = 5
    finally:
        _ae_restore(server, saved)


def t_single_cell_skipped():
    """Single-cell studies are skipped before their data is downloaded (no single-cell blueprint yet)."""
    from hegemon_nb import repos
    from hegemon_nb.config import Config
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    cfg = Config()
    assert repos.looks_single_cell("Single-cell atlas of colonic Tregs in Crohn's disease")
    assert repos.looks_single_cell("x", files=["GSM1_barcodes.tsv.gz"])
    assert not repos.looks_single_cell("Bulk RNA-seq of sorted Tregs (Smart-seq2)")
    # GEO: decided from the SOFT file, nothing else fetched
    root = T / "geo_sc"
    shutil.rmtree(str(root), ignore_errors=True)
    root.mkdir(parents=True)
    _fake_geo_series(root, "GSE999600", "Expression profiling by high throughput sequencing", ["GSM1", "GSM2"],
                     suppl={"GSE999600_RAW.tar": 4096}, title="Single-cell RNA-seq of intestinal Tregs")
    server = start([PY, "-m", "http.server", "8795", "--bind", "127.0.0.1", "--directory", str(root)], 8795)
    real_url = repos.geo_series_url
    repos.geo_series_url = lambda acc: "http://127.0.0.1:8795/{}/".format(acc)
    try:
        out = T / "geo_sc_out"
        shutil.rmtree(str(out), ignore_errors=True)
        try:
            repos.fetch_geo(repos.Dataset("GSE999600", "GEO", "t"), out, cfg)
            raise AssertionError("expected SingleCell")
        except repos.SingleCell:
            pass
        assert [p.name for p in out.iterdir()] == ["GSE999600_family.soft.gz"], list(out.iterdir())
    finally:
        repos.geo_series_url = real_url
        stop(server)
    # ArrayExpress: per-cell files in the study folder -> skipped after the SDRF, the h5ad never fetched
    root = T / "ae_sc"
    shutil.rmtree(str(root), ignore_errors=True)
    (root / "Files").mkdir(parents=True)
    (root / "Files" / "E-MTAB-9004.sdrf.txt").write_text(
        "Source Name\tCharacteristics[disease]\nS1\tcrohn\nS2\tcontrol\n")
    (root / "Files" / "atlas.h5ad").write_bytes(b"h5" * 1000)
    server, saved = _ae_server(root, 8796)
    try:
        out = T / "ae_sc_out"
        shutil.rmtree(str(out), ignore_errors=True)
        try:
            repos.fetch_arrayexpress(repos.Dataset("E-MTAB-9004", "ArrayExpress", "IBD atlas"), out, cfg)
            raise AssertionError("expected SingleCell")
        except repos.SingleCell:
            pass
        assert not (out / "atlas.h5ad").exists()
    finally:
        _ae_restore(server, saved)
    # batch: a single-cell title is left out before ranking
    runner = T / "batch_sc.py"
    runner.write_text(
        "import sys; sys.path.insert(0, {agent!r})\n"
        "from hegemon_nb import repos\n"
        "repos.search_geo = lambda q, limit=40: [repos.Dataset('GSE1', 'GEO', 'Single-cell atlas of Tregs', "
        "kind='Expression profiling by high throughput sequencing', n_samples=50)]\n"
        "repos.search_arrayexpress = lambda q, limit=40: []\n"
        "import hegemon\n"
        "sys.exit(hegemon.main(sys.argv[1:]))\n".format(agent=str(AGENT)))
    p = subprocess.run([PY, str(runner), "batch", "tregs", "--key", "sc"], env=env(), stdin=subprocess.DEVNULL,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=300)
    assert "left out 1 single-cell studies" in p.stdout, p.stdout


def t_status_is_metadata():
    from hegemon_nb.checks import comparison_field
    assert comparison_field("c status_ch1", ["CD"] * 4 + ["control"] * 4, 8), "status is disease status, not technical"
    assert comparison_field("c type_ch1", ["tumor"] * 3 + ["normal"] * 3, 6)
    assert not comparison_field("c title", ["a", "b", "c", "d"], 4)


def t_link_opens_on_a_plot():
    """Links name two genes every dataset under the key has, with cmd=explore; requested topic genes win
    when all datasets have them; relink fixes links published before this existed."""
    from hegemon_nb import native, net, publish
    from hegemon_nb.config import Config
    d = T / "genes"
    shutil.rmtree(str(d), ignore_errors=True)
    d.mkdir(parents=True)

    def idx(name, genes):
        f = d / (name + ".txt")
        f.write_text("ProbeID\tPtr\tName\tDescription\n" + "".join(
            "p{}\t{}\t{}\t\n".format(i, i * 10, g) for i, g in enumerate(genes)))
        return f
    a = idx("a", ["FOXP3", "IL2RA", "ACTB", "GAPDH /// GAPDHP1", "TP53"])
    b = idx("b", ["foxp3", "Il2ra", "ACTB", "GAPDH"])
    c = idx("c", ["FOXP3", "ACTB", "GAPDH"])
    assert publish.choose_genes([a, b], ["FOXP3", "IL2RA"])[:2] == ("FOXP3", "IL2RA")
    assert publish.choose_genes([a, b, c], ["FOXP3", "IL2RA"])[:2] == ("FOXP3", "ACTB"), \
        publish.choose_genes([a, b, c], ["FOXP3", "IL2RA"])
    assert publish.choose_genes([a, c])[:2] == ("ACTB", "GAPDH")

    # relink on a site whose key link predates this: it gains the genes and the page draws on load
    site = php_site()
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    cfg = Config()
    run_dir = T / "genes" / "GSE8100-GPL1"
    _good_trio(run_dir, "GSE8100-GPL1", ("GSM1", "GSM2", "GSM3", "GSM4"))
    native.make_idx(run_dir, "GSE8100-GPL1")
    native.open_permissions(run_dir, T / "genes")
    (site / "explore.conf").write_text((site / "explore.conf").read_text() + (
        "\n[OLD1]\nname= old\nexpr= {0}/GSE8100-GPL1-expr.txt\nindex= {0}/GSE8100-GPL1-idx.txt\n"
        "survival= {0}/GSE8100-GPL1-survival.txt\nindexHeader= {0}/GSE8100-GPL1-ih.txt\nkey= oldkey\n"
        "source= GSE8100\n").format(run_dir))
    (site / "index.html").write_text((site / "index.html").read_text().replace(
        "</table>", '<tr>\n<td>\n<a href="explore.php?key=oldkey"> old set </a>\n</td>\n</tr>\n</table>'))
    php = start(["php", "-S", "127.0.0.1:8765", "-t", str(site)], 8765)
    try:
        rc, out = cli("relink", "--key", "oldkey")
        assert rc == 0 and "opens on ACTB vs GAPDH" in out, out
        index = (site / "index.html").read_text()
        assert 'href="explore.php?key=oldkey&amp;A=ACTB&amp;B=GAPDH&amp;cmd=explore"' in index, index
        page = net.get_text(cfg.site_url + "/explore.php?key=oldkey&A=ACTB&B=GAPDH&cmd=explore", timeout=30, retries=0)
        assert 'value="ACTB"' in page and 'value="GAPDH"' in page and "callExplore();" in page, page[-1500:]
        # running it again changes nothing
        before = (site / "index.html").read_text()
        rc, out = cli("relink", "--key", "oldkey")
        assert rc == 0 and (site / "index.html").read_text() == before
    finally:
        stop(php)


def t_builtin_notebook_conversion():
    """GEO microarrays with their tables in the SOFT file are converted by the notebook's fixed cells with
    no AI call, byte-identical to the lab notebook's own output; anything they cannot do goes to the AI."""
    lab = T / "lab" / "GSE999001-GPL999"
    if not lab.is_dir():
        raise Skip("run run_notebook_reference.py geo first")
    rc, out = cli("build", "GSE999001", "--source", str(T / "fixtures" / "geo"), "--no-publish",
                  HEGEMON_OPENAI_URL="http://127.0.0.1:9/v1/responses",   # dead AI: any call would fail the build
                  HEGEMON_PARALLEL_WRITE_CELLS="1")   # and the expr file goes through the parallel writer
    assert rc == 0, out
    assert "converted by the notebook's own cells; no AI call needed" in out, out
    assert "AI use: 0 call(s)" in out and "conversion (built-in, no AI)" in out, out
    assert re.search(r"time: download \S+ \| reading files \S+ \| conversion \(built-in, no AI\) \S+ \| "
                     r"idx \+ lab script \S+ \| total \S+", out), out
    agent = last_build("GSE999001-*") / "GSE999001-GPL999"
    for kind in ("expr", "idx", "ih"):
        a = (agent / "GSE999001-GPL999-{}.txt".format(kind)).read_bytes()
        b = (lab / "GSE999001-GPL999-{}.txt".format(kind)).read_bytes()
        assert a == b, "built-in {} differs from the notebook's".format(kind)
    receipt = json.loads((last_build("GSE999001-*") / "build.json").read_text())
    assert receipt["attempts"] == 0 and receipt["ai_calls"] == 0 and receipt["timing"]["total"] >= 0

    # a platform table without any gene-symbol column: the built-in script cannot name the genes,
    # so the AI takes over (it can pass symbol_col=...), starting fresh
    import gzip
    src = T / "variants" / "GSE999007"
    src.mkdir(parents=True, exist_ok=True)
    text = gzip.open(str(T / "fixtures" / "geo" / "GSE999001_family.soft.gz"), "rt").read().replace("GSE999001", "GSE999007")
    text = text.replace("ID\tGene Title\tGene Symbol", "ID\tProbeTitle\tProbeSym").replace(
        "#Gene Title = Entrez Gene name\n#Gene Symbol = Entrez Gene symbol", "#ProbeTitle = t\n#ProbeSym = s")
    with gzip.open(str(src / "GSE999007_family.soft.gz"), "wt") as fh:
        fh.write(text)
    fixed = (HERE / "fixture_geo_convert.py").read_text().replace(
        "geo_platform_names(gse.gpls[platform])", "geo_platform_names(gse.gpls[platform], 'ProbeSym', 'ProbeTitle')")
    ai = _mock([{"text": _reply(fixed)}])
    try:
        rc, out = cli("build", "GSE999007", "--source", str(src), "--no-publish")
    finally:
        stop(ai)
    assert rc == 0, out
    assert "the built-in conversion did not pass (crashed); the AI takes over" in out, out
    calls = [json.loads(l) for l in (T / "ai-log.jsonl").read_text().splitlines()]
    assert len(calls) == 1 and "WHAT WENT WRONG" not in calls[0]["body"]["input"], "the AI should start fresh"
    assert calls[0]["body"]["reasoning"]["effort"] == "medium"


def t_ranged_download():
    """Big files arrive as 4 byte ranges at once, identical to a single stream, and a range whose
    connection drops mid-way resumes from where it stopped."""
    import hashlib
    from hegemon_nb import net
    d = T / "ranges"
    shutil.rmtree(str(d), ignore_errors=True)
    d.mkdir(parents=True)
    blob = os.urandom(20 * 1024 * 1024)
    (d / "big.bin").write_bytes(blob)
    want = hashlib.sha256(blob).hexdigest()
    server_py = d / "range_server.py"
    server_py.write_text("""
import os, re, sys
from http.server import HTTPServer, SimpleHTTPRequestHandler
DROPPED = {"done": False}
class H(SimpleHTTPRequestHandler):
    def send_head(self):
        path = self.translate_path(self.path)
        if not os.path.isfile(path):
            return SimpleHTTPRequestHandler.send_head(self)
        size = os.path.getsize(path)
        m = re.match(r"bytes=(\\d+)-(\\d*)", self.headers.get("Range", ""))
        f = open(path, "rb")
        if m:
            s = int(m.group(1)); e = int(m.group(2)) if m.group(2) else size - 1
            f.seek(s)
            self.send_response(206)
            self.send_header("Content-Range", "bytes {}-{}/{}".format(s, e, size))
            self.send_header("Content-Length", str(e - s + 1))
            self.left = e - s + 1
            self.cut = (s > 0 and not DROPPED["done"])   # the first non-zero range drops once, half way
            DROPPED["done"] = DROPPED["done"] or self.cut
        else:
            self.send_response(200)
            self.send_header("Content-Length", str(size))
            self.left, self.cut = size, False
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        return f
    def copyfile(self, src, out):
        sent, stop_at = 0, (self.left // 2 if self.cut else None)
        while self.left > 0:
            chunk = src.read(min(1 << 16, self.left))
            if not chunk:
                break
            out.write(chunk)
            self.left -= len(chunk)
            sent += len(chunk)
            if stop_at is not None and sent >= stop_at:
                self.close_connection = True
                return
    def log_message(self, *a):
        pass
os.chdir(sys.argv[2])
HTTPServer(("127.0.0.1", int(sys.argv[1])), H).serve_forever()
""")
    server = start([PY, str(server_py), "8799", str(d)], 8799)
    plain = start([PY, "-m", "http.server", "8800", "--bind", "127.0.0.1", "--directory", str(d)], 8800)
    used = []
    real = net._download_ranged

    def spy(*a, **k):
        used.append(True)
        return real(*a, **k)
    net._download_ranged = spy
    try:
        n = net.download("http://127.0.0.1:8799/big.bin", str(d / "out1.bin"), retries=2, timeout=5)
        assert used and n == len(blob), "expected the ranged path"
        assert hashlib.sha256((d / "out1.bin").read_bytes()).hexdigest() == want, "ranged copy differs"
        used.clear()
        n = net.download("http://127.0.0.1:8800/big.bin", str(d / "out2.bin"))
        assert not used and hashlib.sha256((d / "out2.bin").read_bytes()).hexdigest() == want
        assert not list(d.glob("*.part")), list(d.glob("*.part"))
    finally:
        net._download_ranged = real
        stop(server)
        stop(plain)


_PICKY_SERVER = r"""
import os, sys, threading, time
from http.server import BaseHTTPRequestHandler, HTTPServer
from socketserver import ThreadingMixIn
ROOT, PORT = sys.argv[2], int(sys.argv[1])
LOCK, ACTIVE, SEEN = threading.Lock(), [0], {}
PAGE = b"<!DOCTYPE html>\n<html><body>Too many requests. Try again later.</body></html>\n"
class S(ThreadingMixIn, HTTPServer):
    daemon_threads = True
class H(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass
    def do_HEAD(self):
        self.serve(False)
    def do_GET(self):
        self.serve(True)
    def serve(self, body_wanted):
        name = self.path.split("?")[0].rsplit("/", 1)[-1]
        with LOCK:
            ACTIVE[0] += 1
            SEEN[name] = SEEN.get(name, 0) + (1 if body_wanted else 0)
            n = SEEN[name]
        try:
            time.sleep(0.4)                      # hold the request so parallel ones overlap
            with LOCK:
                busy = ACTIVE[0] > 1
            if name == "_count":
                body, ctype = repr(SEEN).encode(), "text/plain"
            elif name.startswith("always_page") or (name.startswith("parallel_") and busy) \
                    or (name.startswith("once_") and n == 1):
                body, ctype = PAGE, "text/html"    # the way NCBI answers when it refuses: 200 + a page
            elif os.path.isfile(os.path.join(ROOT, name)):
                body, ctype = open(os.path.join(ROOT, name), "rb").read(), "application/octet-stream"
            else:
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            if body_wanted:
                self.wfile.write(body)
        finally:
            with LOCK:
                ACTIVE[0] -= 1
S(("127.0.0.1", PORT), H).serve_forever()
"""


def t_download_integrity():
    """An error page is never saved as data; servers that refuse parallel requests still deliver every
    file; a source with a file that never arrives is not used, and 'no data' is never claimed."""
    import functools
    import gzip
    from hegemon_nb import net, repos
    from hegemon_nb.config import Config
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    d = T / "picky"
    shutil.rmtree(str(d), ignore_errors=True)
    (d / "www").mkdir(parents=True)
    page = b"<!DOCTYPE html>\n<html><body>Too many requests</body></html>\n"
    (d / "ok.tsv.gz").write_bytes(gzip.compress(b"gene\tGSM1\n"))
    (d / "page.gz").write_bytes(page)
    (d / "page.txt").write_bytes(page)
    (d / "table.txt").write_bytes(b"gene\tGSM1\nGAPDH\t5\n")
    (d / "index.html").write_bytes(page)
    assert net.content_problem(d / "ok.tsv.gz") is None
    assert "HTML page" in net.content_problem(d / "page.gz"), net.content_problem(d / "page.gz")
    assert "HTML page" in net.content_problem(d / "page.txt")
    assert net.content_problem(d / "table.txt") is None and net.content_problem(d / "index.html") is None

    www = d / "www"
    for name in ("once_annot.tsv.gz", "always_page_annot.tsv.gz", "counts.tsv.gz") + tuple(
            "parallel_{}.txt.gz".format(i) for i in range(6)):
        (www / name).write_bytes(gzip.compress("{}\tGAPDH\t5\n".format(name).encode() * 50))
    server_py = d / "picky_server.py"
    server_py.write_text(_PICKY_SERVER)
    server = start([PY, str(server_py), "8802", str(www)], 8802)
    real_download = net.download
    try:
        base = "http://127.0.0.1:8802/"
        # 1. an error page first, the file on the retry: the file is what gets saved
        n = net.download(base + "once_annot.tsv.gz", str(d / "got_once.tsv.gz"), retries=2)
        assert gzip.decompress((d / "got_once.tsv.gz").read_bytes()).startswith(b"once_annot"), "saved the page"
        assert n == (www / "once_annot.tsv.gz").stat().st_size
        # 2. never anything but a page: a clear error, and nothing left behind
        try:
            net.download(base + "always_page_annot.tsv.gz", str(d / "never.tsv.gz"), retries=1)
            raise AssertionError("an HTML page was accepted as a .gz file")
        except net.NetError as err:
            assert "HTML page" in str(err), err
        assert not (d / "never.tsv.gz").exists() and not (d / "never.tsv.gz.part").exists()

        # tests do not need to wait through the real backoff between tries
        repos.net.download = functools.partial(real_download, retries=0)
        cfg = Config()
        # 3. a server that answers parallel requests with a page: everything still arrives, one at a time
        out = T / "runs-link" / "picky_parallel"
        shutil.rmtree(str(out), ignore_errors=True)
        out.mkdir(parents=True)
        files = [("parallel_{}.txt.gz".format(i), base + "parallel_{}.txt.gz".format(i), None) for i in range(6)]
        got = repos._download_all(files, out, repos.Budget(cfg, True), require_all=True)
        assert sorted(got) == sorted(f[0] for f in files), got
        for name, _u, _s in files:
            assert gzip.decompress((out / name).read_bytes()).startswith(name.encode()), name

        # 4. GEO's computed counts whose annotation never arrives: that source is dropped, its leftovers
        #    removed, and the per-sample files are used instead
        root = T / "geo_serve_picky"
        shutil.rmtree(str(root), ignore_errors=True)
        root.mkdir(parents=True)
        seq = "Expression profiling by high throughput sequencing"
        _fake_geo_series(root, "GSE999700", seq, ["GSM7", "GSM8"], per_sample=["counts.txt.gz"])
        _fake_geo_series(root, "GSE999701", seq, ["GSM9", "GSM10"])
        geo = start([PY, "-m", "http.server", "8792", "--bind", "127.0.0.1", "--directory", str(root)], 8792)
        real_url, real_counts = repos.geo_series_url, repos._ncbi_counts_links
        repos.geo_series_url = lambda acc: "http://127.0.0.1:8792/{}/".format(acc)
        repos._ncbi_counts_links = lambda acc: [   # shaped like GEO's own links
            base + "counts.tsv.gz?type=rnaseq_counts&acc={}&format=file&file=GSE_raw_counts.tsv.gz".format(acc),
            base + "always_page_annot.tsv.gz?type=rnaseq_counts&format=file&file=Human.annot.tsv.gz"]
        try:
            out = T / "runs-link" / "picky_geo"
            shutil.rmtree(str(out), ignore_errors=True)
            soft, gse, notes, processed = repos.fetch_geo(repos.Dataset("GSE999700", "GEO", "t"), out, cfg)
            names = {str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()}
            assert "GSE999700_RAW/GSM7_counts.txt.gz" in names and "GSE999700_RAW/GSM8_counts.txt.gz" in names, names
            assert not any("counts.tsv.gz" in n or "annot" in n for n in names if "_RAW/" not in n), \
                "kept half of the counts: " + str(names)
            assert any("computed counts could not be downloaded" in n for n in notes), notes
            for p in out.rglob("*.gz"):
                assert net.content_problem(p) is None, p
            # 5. nothing else to fall back on: the failure is handed on (to the AI), never "no processed data"
            out = T / "runs-link" / "picky_geo2"
            shutil.rmtree(str(out), ignore_errors=True)
            ds2 = repos.Dataset("GSE999701", "GEO", "t")
            soft, gse, notes, processed = repos.fetch_geo(ds2, out, cfg)
            assert not processed, processed
            assert any("computed counts" in f and "HTML page" in f for f in ds2.failures), ds2.failures
            assert "Human.annot.tsv.gz" in ds2.catalog and any("Download problems" in n for n in notes), notes
        finally:
            repos.geo_series_url, repos._ncbi_counts_links = real_url, real_counts
            stop(geo)
    finally:
        repos.net.download = real_download
        stop(server)


def t_fast_paths_match():
    """Every speed-up gives exactly what the plain pandas code gives: the parallel expr writer, the
    side-by-side sample matrix, and the quick number check on expr files."""
    import hashlib
    import numpy as np
    import pandas as pd
    from hegemon_nb import checks, helpers
    d = T / "fast_paths"
    shutil.rmtree(str(d), ignore_errors=True)
    d.mkdir(parents=True)
    rng = np.random.default_rng(5)

    # 1. parallel writing: byte-identical to one to_csv, whatever the column types
    n = 3000
    frame = pd.DataFrame(rng.normal(7, 2, (n, 40)), columns=["GSM%d" % i for i in range(40)])
    frame.iloc[3, 2], frame.iloc[4, 2], frame.iloc[5, 2] = np.nan, np.inf, -0.0
    frame.iloc[6, 2], frame.iloc[7, 2], frame.iloc[8, 2] = 1e-7, 1e17, 123456789012345.6
    frame.insert(0, "Name", ['G%d: a "quoted", name' % i for i in range(n)])
    frame.insert(0, "ProbeID", ["%d_at" % i for i in range(n)])
    frame["ints"] = rng.integers(0, 100, n)
    frame["objs"] = ["x"] * (n - 1) + [None]
    frame.to_csv(str(d / "one.txt"), sep="\t", index=False)
    old = os.environ.get("HEGEMON_PARALLEL_WRITE_CELLS")
    os.environ["HEGEMON_PARALLEL_WRITE_CELLS"] = "1"
    try:
        helpers._to_csv(frame, str(d / "par.txt"))
    finally:
        if old is None:
            os.environ.pop("HEGEMON_PARALLEL_WRITE_CELLS")
        else:
            os.environ["HEGEMON_PARALLEL_WRITE_CELLS"] = old
    assert (d / "one.txt").read_bytes() == (d / "par.txt").read_bytes(), "parallel writer changed the bytes"

    # 2. the sample matrix: side by side when probe lists match, pandas' outer merge otherwise
    class _G(object):
        def __init__(self, tables):
            self.gsms = OrderedDict((k, type("M", (), {"table": v, "metadata": {"platform_id": ["P"]}})())
                                    for k, v in tables.items())
    ids = ["p%d" % i for i in range(300)]
    cases = [
        {"G%d" % k: pd.DataFrame({"ID_REF": ids, "VALUE": rng.normal(size=300)}) for k in range(4)},
        {"G%d" % k: pd.DataFrame({"ID_REF": list(range(300)), "VALUE": rng.integers(0, 50, 300)}) for k in range(4)},
        {"G0": pd.DataFrame({"ID_REF": ids, "VALUE": rng.normal(size=300)}),
         "G1": pd.DataFrame({"ID_REF": ids[::-1], "VALUE": rng.normal(size=300)})},
        {"G0": pd.DataFrame({"ID_REF": ids, "VALUE": rng.normal(size=300)}),
         "G1": pd.DataFrame({"ID_REF": ids[:-3], "VALUE": rng.normal(size=297)})},
        {"G%d" % k: pd.DataFrame({"ID_REF": ids, "VALUE": ["1.5"] * 299 + ["null"]}) for k in range(3)},
        {"G%d" % k: pd.DataFrame({"ID_REF": ids[:5] + ids[:5], "VALUE": list(range(10))}) for k in range(3)},
        {"G0": pd.DataFrame({"ID_REF": [1, 2, 3], "VALUE": [1.0, 2.0, 3.0]}),
         "G1": pd.DataFrame({"ID_REF": [1.0, 2.0, 3.0], "VALUE": [4.0, 5.0, 6.0]})},
    ]
    for tables in cases:
        fast = helpers.geo_sample_matrix(_G(tables), "P")
        pieces = [t[["ID_REF", "VALUE"]].drop_duplicates(subset=["ID_REF"], keep="first")
                  .set_index("ID_REF")["VALUE"].rename(k) for k, t in tables.items()]
        plain = pd.concat(pieces, axis=1, join="outer", sort=False).apply(pd.to_numeric, errors="coerce").fillna(0)
        plain.index.name = "ID"
        plain = plain.reset_index()
        assert fast.to_csv(sep="\t", index=False) == plain.to_csv(sep="\t", index=False), list(tables)
        assert list(fast.dtypes) == list(plain.dtypes), list(tables)

    # 3. the number check: the quick parser and the cell-by-cell one agree on every kind of file
    header = "ProbeID\tName\tGSM1\tGSM2\n"
    files = {
        "plain": "p1\tA:a\t1.5\t2.25\np2\tB:b\t-0.5\t7\n",
        "blanks": "p1\tA:a\t\t2\np2\tB:b\t3\t\n",
        "inf": "p1\tA:a\tinf\t2\np2\tB:b\t3\t-Infinity\n",
        "counts": "p1\tA:a\t150\t2\np2\tB:b\t3\t900\n",
        "unlogged": "p1\tA:a\t150.5\t2\np2\tB:b\t3\t900\n",
        "booleans": "p1\tA:a\tTRUE\t2\np2\tB:b\t3\tfalse\n",
        "words": "p1\tA:a\tnan\t2\np2\tB:b\tNA\t1,5\n",
        "spaces": "p1\tA:a\t 5\t2 \np2\tB:b\t  \t4\n",
        "fields": "p1\tA:a\t1\t2\t3\np2\tB:b\t4\n",
        "blank names": "p1\t\t1\t2\np2\t\t3\t4\n",
    }
    real = checks._expr_numbers
    for label, body in files.items():
        folder = d / label.replace(" ", "_")
        folder.mkdir()
        (folder / "X-expr.txt").write_text(header + body)
        (folder / "X-survival.txt").write_text("ArrayId\ttime\tstatus\tc tissue\nGSM1\t\t\ta\nGSM2\t\t\tb\n")
        (folder / "X-ih.txt").write_text("ArrayID\tArrayHeader\tClinicalHeader\nGSM1\tGSM1\ts1\nGSM2\tGSM2\ts2\n")
        quick = checks.check_core(folder, "X", "Other")
        checks._expr_numbers = lambda *a, **k: None
        try:
            slow = checks.check_core(folder, "X", "Other")
        finally:
            checks._expr_numbers = real
        assert (quick.errors, quick.warnings, quick.stats) == (slow.errors, slow.warnings, slow.stats), \
            (label, quick.errors, slow.errors, quick.stats, slow.stats)


_ANY_SITE_SCRIPT = """import os
import pandas as pd
from hegemon_nb.helpers import cpm_log2, write_native
counts = pd.read_csv(os.path.join(SOURCE_DIR, "counts.tsv"), sep="\\t", dtype={"gene": str})
sheet = pd.read_csv(os.path.join(SOURCE_DIR, "samples.tsv"), sep="\\t", dtype=str, keep_default_na=False)
known = [c for c in counts.columns[1:] if c in set(sheet["sample"])]
left_out = [c for c in counts.columns[1:] if c not in set(sheet["sample"])]
values = cpm_log2(counts[known])
expr = pd.concat([counts["gene"].rename("ProbeID"), counts["gene"].rename("Name"), values], axis=1)
sheet = sheet.set_index("sample").loc[known].reset_index()
survival = pd.DataFrame({"ArrayId": known, "time": "", "status": "", "c group": sheet["group"].tolist(),
                         "c patient": sheet["patient"].tolist()})
ih = pd.DataFrame({"ArrayID": known, "ArrayHeader": known, "ClinicalHeader": known})
write_native(PREFIX, expr, survival, ih)
print("matched", len(known), "samples; left out columns without metadata:", left_out)
"""


def t_any_site_fetch():
    """A study on a website the agent has never seen: the page's links become the catalog, the AI asks for
    the files it needs (after being told a missing file is not a reason to stop), extra data columns without
    metadata are left out, and it builds."""
    import gzip
    site = T / "site_any"
    shutil.rmtree(str(site), ignore_errors=True)
    (site / "data").mkdir(parents=True)
    (site / "study.html").write_text(
        "<html><head><title>Psoriasis study 7</title></head><body><p>Bulk RNA-seq of lesional (LS) and "
        "non-lesional (NL) skin.</p><a href='data/counts.tsv'>Gene counts</a> <a href='data/samples.tsv'>Sample "
        "sheet</a> <a href='data/raw_R1.fastq.gz'>raw reads</a></body></html>")
    genes = ["GAPDH", "ACTB", "IL17A", "KRT16", "IL36G", "S100A7", "DEFB4A", "PI3"]
    rows = ["gene\t" + "\t".join("S%d" % i for i in range(1, 8))]
    for k, g in enumerate(genes):
        rows.append(g + "\t" + "\t".join(str(100 + 37 * k + 11 * i) for i in range(1, 8)))
    (site / "data" / "counts.tsv").write_text("\n".join(rows) + "\n")
    (site / "data" / "samples.tsv").write_text(
        "sample\tgroup\tpatient\n" + "".join("S{}\t{}\tP{}\n".format(i, "LS" if i % 2 else "NL", (i + 1) // 2)
                                            for i in range(1, 7)))   # S7 has data but no metadata
    (site / "data" / "raw_R1.fastq.gz").write_bytes(gzip.compress(b"@r1\nACGT\n+\nIIII\n" * 1000))
    web = start([PY, "-m", "http.server", "8803", "--bind", "127.0.0.1", "--directory", str(site)], 8803)
    ai = _mock([{"text": "UNCONVERTIBLE: MISSING_FILES: the sample sheet is not in SOURCE_DIR"},
                {"text": "The counts matrix and the sample sheet are on the study page.\nFETCH: counts.tsv\nFETCH: samples.tsv"},
                {"text": "Source: counts.tsv.\nSamples: matched to samples.tsv by ID; S7 has no metadata.\n"
                         "Normalization: cpm_log2.\n```python\n" + _ANY_SITE_SCRIPT + "```"}])
    try:
        rc, out = cli("build", "MYSTUDY7", "--url", "http://127.0.0.1:8803/study.html", "--no-publish")
    finally:
        stop(ai)
        stop(web)
    assert rc == 0, out
    assert "is a web page: 3 link(s) on it are offered to the AI" in out, out
    assert "the AI said a file is missing; it was told to ask for it instead" in out, out
    assert "The AI asked for 2 more file(s)" in out and "counts.tsv: downloaded" in out, out
    assert "files pass the checks: 6 samples" in out, out
    work = last_build("MYSTUDY7-*")
    assert not (work / "source" / "raw_R1.fastq.gz").exists(), "downloaded raw reads nobody asked for"
    calls = [json.loads(l) for l in (T / "ai-log.jsonl").read_text().splitlines()]
    assert len(calls) == 3, len(calls)
    first, second, third = (c["body"]["input"] for c in calls)
    assert "FILES THIS STUDY OFFERS THAT ARE NOT DOWNLOADED" in first and "Psoriasis study 7" in first
    assert "A missing file is not a reason to stop" in second
    assert "Result of your FETCH request" in third and "counts.tsv: downloaded" in third and "S7" in third
    assert [c["body"]["reasoning"]["effort"] for c in calls] == ["medium"] * 3, "nothing failed, so no high effort"
    receipt = json.loads((work / "build.json").read_text())
    assert receipt["fetched_by_ai"] == ["counts.tsv", "samples.tsv"], receipt.get("fetched_by_ai")
    survival = (work / "MYSTUDY7" / "MYSTUDY7-survival.txt").read_text()
    assert "S7" not in survival and "c group" in survival

    # asking by name pattern, and a URL off the list that points at a private address is refused
    from hegemon_nb import repos
    from hegemon_nb.config import Config
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    (site / "GSE1_RAW").mkdir()
    for g in ("GSM1", "GSM2"):
        (site / "GSE1_RAW" / (g + "_counts.txt.gz")).write_bytes(gzip.compress(b"gene\t1\n"))
    web = start([PY, "-m", "http.server", "8803", "--bind", "127.0.0.1", "--directory", str(site)], 8803)
    try:
        ds = repos.Dataset("GSE1", "GEO", "t")
        for g in ("GSM1", "GSM2"):
            repos.catalog_add(ds, "GSE1_RAW/{}_counts.txt.gz".format(g),
                              "http://127.0.0.1:8803/GSE1_RAW/{}_counts.txt.gz".format(g), None, "per-sample file")
        out_dir = T / "runs-link" / "fetch_pattern"
        shutil.rmtree(str(out_dir), ignore_errors=True)
        out_dir.mkdir(parents=True)
        report, got = repos.fetch_requested(ds, out_dir, Config(), [
            "GSE1_RAW/*_counts.txt.gz", "http://127.0.0.1:8803/data/counts.tsv", "nothere.txt"])
    finally:
        stop(web)
    text = "\n".join(report)
    assert sorted(got) == ["GSE1_RAW/GSM1_counts.txt.gz", "GSE1_RAW/GSM2_counts.txt.gz"], (got, text)
    assert "private or local address" in text and "nothere.txt: not in the list" in text, text


def t_ai_picks_files():
    """When a study offers several big alternative files, their first lines are read and the AI picks what the
    conversion needs; only those are downloaded and the rest stay available. Big studies are capped by default."""
    import gzip
    from hegemon_nb import repos
    from hegemon_nb.config import Config
    sys.path.insert(0, str(AGENT))
    import hegemon
    for k, v in env().items():
        if k.startswith("HEGEMON_"):
            os.environ[k] = v
    cfg = Config()
    seen = []

    def picker(ds, files, what):
        seen.append((what, files))
        return [n for n, _size, _head in files if "raw_counts" in n or n.startswith("GSM#_counts")], "raw counts"

    # 1. ArrayExpress: three matrices named in the SDRF (skin raw, skin normalized, blood raw)
    root = T / "ae_pick"
    shutil.rmtree(str(root), ignore_errors=True)
    (root / "Files").mkdir(parents=True)
    rows = ["Source Name\tCharacteristics[organism part]\tDerived Array Data File"]
    for i in range(1, 9):
        part = "skin" if i <= 4 else "blood"
        rows.append("S{}\t{}\t{}".format(i, part, "Skin_raw_counts.txt" if part == "skin" else "Blood_raw_counts.txt"))
        rows.append("S{}\t{}\t{}".format(i, part, "Skin_norm_counts.txt" if part == "skin" else "Blood_raw_counts.txt"))
    (root / "Files" / "E-MTAB-9010.sdrf.txt").write_text("\n".join(rows) + "\n")
    for name in ("Skin_raw_counts.txt", "Skin_norm_counts.txt", "Blood_raw_counts.txt"):
        (root / "Files" / name).write_text("gene\tS1\tS2\tS3\tS4\n" + "GAPDH\t10\t20\t30\t40\n" * 5000)
    server, saved = _ae_server(root, 8795)
    real_minimum = repos.PICK_ABOVE
    repos.PICK_ABOVE = 1000   # the test files are small
    try:
        out = T / "runs-link" / "ae_pick_out"
        shutil.rmtree(str(out), ignore_errors=True)
        ds = repos.Dataset("E-MTAB-9010", "ArrayExpress", "t")
        notes, got = repos.fetch_arrayexpress(ds, out, cfg, picker=picker)
        names = {p.name for p in out.rglob("*") if p.is_file()}
        assert {"Skin_raw_counts.txt", "Blood_raw_counts.txt"} <= names, names
        assert "Skin_norm_counts.txt" not in names, "downloaded the file the AI did not choose"
        assert "Skin_norm_counts.txt" in ds.catalog and not ds.failures, (list(ds.catalog), ds.failures)
        what, files = seen[-1]
        assert "gene\tS1\tS2" in dict((n, h) for n, _s, h in files)["Skin_norm_counts.txt"], "no first lines shown"
        # without a picker (or under the size threshold) everything named is downloaded, as before
        shutil.rmtree(str(out), ignore_errors=True)
        repos.fetch_arrayexpress(repos.Dataset("E-MTAB-9010", "ArrayExpress", "t"), out, cfg)
        assert (out / "Skin_norm_counts.txt").is_file()
        # the SDRF decides --max-samples before any data file is downloaded
        shutil.rmtree(str(out), ignore_errors=True)
        try:
            repos.fetch_arrayexpress(repos.Dataset("E-MTAB-9010", "ArrayExpress", "t"), out, cfg, max_samples=5)
            raise AssertionError("expected TooManySamples")
        except repos.TooManySamples as err:
            assert err.count == 8 and err.maximum == 5
        assert not (out / "Skin_raw_counts.txt").exists()
    finally:
        repos.PICK_ABOVE = real_minimum
        _ae_restore(server, saved)

    # 2. GEO per-sample files of two kinds: the AI picks the kind
    geo_root = T / "geo_pick"
    shutil.rmtree(str(geo_root), ignore_errors=True)
    geo_root.mkdir(parents=True)
    _fake_geo_series(geo_root, "GSE999800", "Expression profiling by high throughput sequencing", ["GSM81", "GSM82"],
                     per_sample=["counts.txt.gz", "fpkm.txt.gz"])
    geo = start([PY, "-m", "http.server", "8792", "--bind", "127.0.0.1", "--directory", str(geo_root)], 8792)
    real_url, real_counts = repos.geo_series_url, repos._ncbi_counts_links
    repos.geo_series_url = lambda acc: "http://127.0.0.1:8792/{}/".format(acc)
    repos._ncbi_counts_links = lambda acc: []
    try:
        out = T / "runs-link" / "geo_pick_out"
        shutil.rmtree(str(out), ignore_errors=True)
        ds = repos.Dataset("GSE999800", "GEO", "t")
        soft, gse, notes, processed = repos.fetch_geo(ds, out, cfg, picker=picker)
        names = {str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()}
        assert "GSE999800_RAW/GSM81_counts.txt.gz" in names and "GSE999800_RAW/GSM82_counts.txt.gz" in names, names
        assert not any("fpkm" in n for n in names), names
        try:
            repos.fetch_geo(repos.Dataset("GSE999800", "GEO", "t"), T / "runs-link" / "geo_pick_cap", cfg, max_samples=1)
            raise AssertionError("expected TooManySamples")
        except repos.TooManySamples as err:
            assert err.count == 2
    finally:
        repos.geo_series_url, repos._ncbi_counts_links = real_url, real_counts
        stop(geo)

    # 3. the fast model's reply is parsed and only listed names are accepted
    class _FakeLLM(object):
        def ask(self, instructions, text, **kw):
            assert "Skin_norm_counts.txt" in text and "gene\tS1" in text and kw["model"] == cfg.fast_model
            return 'Here: {"download": ["Skin_raw_counts.txt", "not_listed.txt"], "why": "raw counts only"}'
    chosen, why = hegemon.make_picker(cfg, _FakeLLM())(repos.Dataset("E-X", "ArrayExpress", "t"),
                                                       [("Skin_raw_counts.txt", 100, "gene\tS1"),
                                                        ("Skin_norm_counts.txt", 90, "gene\tS1")], "files")
    assert chosen == ["Skin_raw_counts.txt"] and why == "raw counts only", (chosen, why)


def main():
    if not (T / "fixtures" / "geo").is_dir():
        subprocess.run([PY, str(HERE / "make_fixtures.py"), str(T / "fixtures")], check=True)
    for d in ("runs", "home"):
        if (T / d).exists():
            shutil.rmtree(str(T / d))
    tests = [("SOFT reader gives the same objects as GEOparse", t_soft_reader),
             ("helpers match the notebook cells (with the documented fixes)", t_helpers),
             ("AI code screening", t_guard),
             ("checks catch bad files", t_checks),
             ("GEO build -> idx -> lab-script stand-in -> publish -> live PHP checks", t_geo_build_publish),
             ("GEO build vs the lab notebook run on the same files", t_geo_vs_notebook),
             ("publish failure restores explore.conf and index.html", t_publish_rollback),
             ("AI repair loop (mock OpenAI endpoint)", t_ai_repair_loop),
             ("AI provider errors and UNCONVERTIBLE are reported as such", t_ai_outcomes),
             ("ArrayExpress build vs the lab notebook run on the same files", t_ae_build_vs_notebook),
             ("real Hegemon sample data: notebook idx + publish + PHP lookups", t_real_data_through_php),
             ("GEO page parsers and AI reply parsing", t_repo_parsers),
             ("shared key: many datasets on one Hegemon page", t_shared_key),
             ("metadata plot check refuses a dataset Hegemon cannot group", t_metadata_plot_check),
             ("batch: 10-under-one-key, survives failures, continues on re-run", t_batch_shared_key),
             ("ArrayExpress SDRF is mandatory and complete", t_arrayexpress_sdrf),
             ("diagnose names the unplottable dataset already on a site", t_diagnose),
             ("GEO downloads the cheapest source, not the big bundle", t_geo_acquisition_plan),
             ("metadata the AI dropped is named and put back", t_dropped_metadata_is_repaired),
             ("no groupable metadata -> NO_GROUPING_METADATA, never published", t_no_grouping_metadata),
             ("studies under --min-samples are left out before downloading", t_min_samples),
             ("unpublish by section or key, with link cleanup and backup", t_unpublish),
             ("ArrayExpress: named file found in the archive; nothing else downloaded", t_arrayexpress_archives_and_limits),
             ("single-cell studies skipped before download", t_single_cell_skipped),
             ("a field named status counts as metadata", t_status_is_metadata),
             ("links open on a plot with genes every dataset has; relink", t_link_opens_on_a_plot),
             ("built-in notebook conversion: no AI, identical to the notebook; AI takes over", t_builtin_notebook_conversion),
             ("big downloads in 4 parallel ranges, resuming a dropped range", t_ranged_download),
             ("error pages are never saved as data; picky servers still deliver", t_download_integrity),
             ("speed-ups give exactly what plain pandas gives", t_fast_paths_match),
             ("a study on a site never seen before: page links, AI asks for files, extras left out", t_any_site_fetch),
             ("the AI picks between big alternative files before they download", t_ai_picks_files)]
    only = sys.argv[1:]
    for name, fn in tests:
        if not only or any(o in fn.__name__ for o in only):
            run(name, fn)
    print("\n{} passed, {} failed, {} skipped (Python {})".format(
        sum(r[1] == "PASS" for r in RESULTS), sum(r[1] == "FAIL" for r in RESULTS),
        sum(r[1] == "SKIP" for r in RESULTS), sys.version.split()[0]))
    return 1 if any(r[1] == "FAIL" for r in RESULTS) else 0


if __name__ == "__main__":
    sys.exit(main())
