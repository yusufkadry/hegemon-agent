"""FAKE thr/info/vinfo/bv for local tests only (formats follow Hegemon's readers, values are simple)."""
import sys
import numpy as np

pre = sys.argv[1]
with open(pre + "-expr.txt") as fh, open(pre + "-thr.txt", "w") as thr, open(pre + "-info.txt", "w") as info, \
        open(pre + "-vinfo.txt", "w") as vinfo, open(pre + "-bv.txt", "w") as bv:
    header = fh.readline().rstrip("\n").split("\t")
    info.write("ProbeID\tName\tthr\tmean\tmean-thr\tperc\tmin\tmax\tsd\n")
    vinfo.write("ProbeID\tName\tvar\n")
    bv.write("ProbeID\tName\tBitVector\n")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        v = np.array([float(x) if x else np.nan for x in p[2:]])
        t = float(np.nanmedian(v)) if np.isfinite(v).any() else 0.0
        thr.write("{}\t{}\t0\t{}\t{}\n".format(p[0], t, t - 0.5, t + 0.5))
        info.write("{}\t{}\t{}\t{}\t0\t0.5\t{}\t{}\t0\n".format(p[0], p[1], t, t, np.nanmin(v), np.nanmax(v)))
        vinfo.write("{}\t{}\t{}\n".format(p[0], p[1], float(np.nanvar(v))))
        bv.write("{}\t{}\t{}\n".format(p[0], p[1], "".join("2" if x > t + 0.5 else ("0" if x < t - 0.5 else "1") for x in v)))
print("FAKE support files written for", pre)
