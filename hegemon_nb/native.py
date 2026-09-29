"""The notebooks' fixed final cells, run exactly as the lab runs them.

make_idx         = the idx cell, copied from both notebooks
run_lab_script   = cmd = f'bash {gse_processing} {cwd}', run inside a folder
                   whose name is the file prefix (that is what the notebooks do)
check_support    = reads thr/info/vinfo/bv the way Hegemon's own code does
"""
import hashlib
import os
import subprocess
from collections import deque
from pathlib import Path

from . import ui
from .checks import Result

KINDS = ("expr", "idx", "survival", "ih", "thr", "info", "vinfo", "bv")


def make_idx(run_dir, prefix):
    # ---- copied from the notebooks' "idx file" cell; only the paths are joined to run_dir ----
    ptr = []
    ids = []
    name = []
    desc = []
    pos = 0
    with open(os.path.join(str(run_dir), '{}-expr.txt'.format(prefix)), 'rb') as f:
        for line in f:
            if pos == 0:
                pos += len(line)
            else:
                ptr.append(pos)
                pos += len(line)
                split = line.decode("utf-8").split('\t')
                ids.append(split[0])
                name.append(split[1].split(':')[0])
                desc.append(':'.join(split[1].split(':')[1:]))
    with open(os.path.join(str(run_dir), '{}-idx.txt'.format(prefix)), 'w') as f:
        f.write('ProbeID\tPtr\tName\tDescription\n')
        for i in range(len(ids)):
            f.write('{}\t{}\t{}\t{}\n'.format(ids[i], ptr[i], name[i], desc[i]))
    # ---- end of notebook code ----
    return len(ids)


def check_idx(run_dir, prefix):
    """Every idx row must have 4 fields (Hegemon skips others) and point at its own expr row."""
    run_dir = Path(run_dir)
    with open(str(run_dir / "{}-idx.txt".format(prefix)), "rb") as fi, \
            open(str(run_dir / "{}-expr.txt".format(prefix)), "rb") as fe:
        fi.readline()
        for n, line in enumerate(fi, start=2):
            parts = line.rstrip(b"\r\n").split(b"\t")
            if len(parts) != 4:
                return "idx line {} has {} fields; Hegemon ignores rows without exactly 4".format(n, len(parts))
            fe.seek(int(parts[1]))
            if not fe.readline().startswith(parts[0] + b"\t"):
                return "idx line {}: pointer {} does not land on probe {}".format(n, parts[1].decode(), parts[0].decode())
    return None


def run_lab_script(cfg, run_dir, prefix, log_path):
    """Notebook: gse_processing = '.../jupyter_gse_processing'; cwd = folder name; !bash {gse_processing} {cwd}"""
    run_dir = Path(run_dir)
    if run_dir.name != prefix:
        raise ValueError("The lab script needs the folder to be named {}".format(prefix))
    with open(str(log_path), "w") as log, ui.Heartbeat("jupyter_gse_processing running", every=60):
        try:
            proc = subprocess.run(["bash", str(cfg.lab_script), run_dir.name], cwd=str(run_dir),
                                  stdout=log, stderr=subprocess.STDOUT, timeout=cfg.lab_timeout)
            return proc.returncode
        except subprocess.TimeoutExpired:
            return "timeout after {:.0f}s".format(cfg.lab_timeout)


def tail(path, n=25):
    try:
        with open(str(path), "r", errors="replace") as fh:
            return "".join(deque(fh, maxlen=n))
    except OSError:
        return ""


def _probe_ids(expr_path):
    ids = set()
    with open(str(expr_path), "rb") as fh:
        fh.readline()
        for line in fh:
            ids.add(line.split(b"\t", 1)[0].decode("utf-8", "replace"))
    return ids


def check_support(run_dir, prefix, n_samples):
    """Read thr/info/vinfo/bv: right probes, right shape, made after this run's expr."""
    res = Result()
    run_dir = Path(run_dir)
    expr = run_dir / "{}-expr.txt".format(prefix)
    ids = _probe_ids(expr)
    for kind in ("thr", "info", "vinfo", "bv"):
        path = run_dir / "{}-{}.txt".format(prefix, kind)
        if not path.is_file() or path.stat().st_size == 0:
            res.errors.append("{} is missing or empty".format(path.name))
            continue
        if path.stat().st_mtime + 1 < expr.stat().st_mtime:
            res.errors.append("{} is older than the expr file (left over from an earlier run)".format(path.name))
        rows = matched = foreign = bad = 0
        lengths = set()
        foreign_examples = []
        with open(str(path), "rb") as fh:
            for n, raw in enumerate(fh):
                parts = raw.rstrip(b"\r\n").decode("utf-8", "replace").split("\t")
                if n == 0 and parts[0] not in ids:
                    continue  # header line
                if not parts[0]:
                    continue
                rows += 1
                if parts[0] in ids:
                    matched += 1
                else:
                    foreign += 1
                    if len(foreign_examples) < 3:
                        foreign_examples.append(parts[0])
                if kind == "thr":
                    try:
                        [float(x) for x in parts[1:5]]
                        if len(parts) < 5:
                            raise ValueError
                    except ValueError:
                        bad += 1
                elif kind == "bv":
                    if len(parts) < 3:
                        bad += 1
                    else:
                        lengths.add(len(parts[2]))
                elif kind == "info" and len(parts) < 9:
                    bad += 1
        coverage = matched / float(len(ids) or 1)
        res.stats[kind] = "{} rows, {:.0%} of probes".format(rows, coverage)
        if foreign > max(2, 0.01 * rows):
            res.errors.append("{}: {} rows name probes that are not in this expr file (e.g. {}); "
                              "output from a different dataset?".format(path.name, foreign, ", ".join(foreign_examples)))
        if kind == "thr" and coverage < 0.5:
            res.errors.append("{} covers only {:.0%} of the probes".format(path.name, coverage))
        elif coverage < 0.5:
            res.warnings.append("{} covers {:.0%} of the probes".format(path.name, coverage))
        if kind == "thr" and bad > max(2, 0.01 * rows):
            res.errors.append("{}: {} rows lack 4 numeric thresholds".format(path.name, bad))
        if kind == "info" and bad:
            res.warnings.append("{}: {} rows have fewer than 9 fields (Hegemon skips them)".format(path.name, bad))
        if kind == "bv":
            if bad > max(2, 0.01 * rows) or len(lengths) > 1:
                res.errors.append("{}: bit vectors are malformed or of unequal length".format(path.name))
            elif lengths and lengths != {n_samples}:
                res.warnings.append("{}: bit vectors are {} long but there are {} samples"
                                    .format(path.name, lengths.pop(), n_samples))
    return res


def open_permissions(run_dir, runs_root):
    """The web server must be able to read the files and walk the folders above them."""
    run_dir, runs_root = Path(run_dir), Path(runs_root)
    for p in [run_dir] + list(run_dir.rglob("*")):
        try:
            os.chmod(str(p), 0o755 if p.is_dir() else 0o644)
        except OSError:
            pass
    try:
        run_dir.relative_to(runs_root)
    except ValueError:
        return  # never change permissions outside the agent's runs folder
    parent = run_dir.parent
    while True:
        try:
            os.chmod(str(parent), os.stat(str(parent)).st_mode | 0o055)
        except OSError:
            pass
        if parent == runs_root or parent == parent.parent:
            break
        parent = parent.parent


def file_hashes(run_dir, prefix):
    out = {}
    for kind in KINDS:
        digest = hashlib.sha256()
        with open(str(Path(run_dir) / "{}-{}.txt".format(prefix, kind)), "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                digest.update(chunk)
        out[kind] = digest.hexdigest()
    return out
