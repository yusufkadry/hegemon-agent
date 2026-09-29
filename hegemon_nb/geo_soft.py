"""Read a GEO family SOFT file the same way GEOparse does.

The lab notebook calls GEOparse.get_GEO(), which downloads GSExxx_family.soft.gz
and parses it. GEOparse is not installed in the server venv, so this module
reproduces its parsing rules (entry splitting, metadata keys, '#' column
descriptions, tab-separated tables) and exposes the same attributes the notebook
uses: gse.metadata, gse.gsms[name].metadata/.table/.columns, gse.gpls[name].table.
"""
from collections import OrderedDict, defaultdict
from io import BytesIO, StringIO
import gzip
import re
import warnings

import pandas as pd
try:
    from pandas.errors import DtypeWarning
except ImportError:  # very old pandas
    DtypeWarning = Warning


class _Entry(object):
    def __init__(self, name, metadata, columns=None, table=None):
        self.name = name
        self.metadata = metadata
        self.columns = columns if columns is not None else pd.DataFrame(columns=["description"])
        self.table = table if table is not None else pd.DataFrame()
        # Same correction GEOparse applies: rename duplicate column names the
        # way pandas does, then order descriptions like the table header.
        index = self.columns.index.tolist()
        header = [str(c) for c in self.table.columns]
        if index != header and len(self.table.columns):
            if not self.columns.index.is_unique:
                seen, fixed = {}, []
                for name_ in index:
                    if name_ not in seen:
                        seen[name_] = 0
                        fixed.append(name_)
                    else:
                        seen[name_] += 1
                        fixed.append("%s.%i" % (name_, seen[name_]))
                self.columns.index = fixed
            if self.columns.index.tolist() != header and sorted(self.columns.index.tolist()) == sorted(header):
                self.columns = self.columns.loc[self.table.columns]

    def __repr__(self):
        return "<{} {}>".format(type(self).__name__, self.name)


class GSM(_Entry):
    pass


class GPL(_Entry):
    pass


class GSE(object):
    def __init__(self, name, metadata, gsms, gpls):
        self.name = name
        self.metadata = metadata
        self.gsms = gsms
        self.gpls = gpls

    def __repr__(self):
        return "<GSE {}: {} samples, {} platforms>".format(self.name, len(self.gsms), len(self.gpls))


def _parse_entry(line):
    # Identical to GEOparse.__parse_entry
    if line.startswith("!"):
        line = re.sub(r"!\w*?_", "", line)
    else:
        line = line.strip()[1:]
    parts = [part.strip() for part in line.split("=", 1)]
    if len(parts) == 2:
        return parts[0], parts[1]
    return parts[0], ""


def _metadata(lines):
    meta = defaultdict(list)
    for line in lines:
        line = line.rstrip()
        if line.startswith("!"):
            if "_table_begin" in line or "_table_end" in line:
                continue
            key, value = _parse_entry(line)
            meta[key].append(value)
    return dict(meta)


def _columns(lines):
    names, descriptions = [], []
    for line in lines:
        line = line.rstrip()
        if line.startswith("#"):
            key, value = _parse_entry(line)
            names.append(key)
            descriptions.append(value)
    return pd.DataFrame(descriptions, index=names, columns=["description"])


_LINE_END_SPACE = re.compile(r"[^\S\n]+(?=\n)")
_RARE_SPACE = (b"\x0b", b"\x0c", b"\x1c", b"\x1d", b"\x1e", b"\x1f")


def _table(data):
    """GEOparse's parse_table_data, on the table's lines as one block of bytes.

    GEOparse keeps each line that does not start with ^ ! # and is not blank, strips trailing
    whitespace, and reads the result with read_csv. Blank lines are skipped by read_csv itself,
    so the lines only need stripping when some line actually ends in whitespace (rare).
    """
    if not data.endswith(b"\n"):
        data += b"\n"
    plain = data.isascii() and b"\r" not in data and not any(c in data for c in _RARE_SPACE) \
        and b" \n" not in data and b"\t\n" not in data
    with warnings.catch_warnings():
        # GEOparse reads with pandas' defaults; the mixed-type warning it prints is only noise
        warnings.simplefilter("ignore", DtypeWarning)
        if plain:
            if not data.strip():
                return pd.DataFrame()
            return pd.read_csv(BytesIO(data), index_col=None, sep="\t")
        text = data.decode("utf-8", errors="ignore").replace("\r\n", "\n").replace("\r", "\n")
        text = _LINE_END_SPACE.sub("", text)
        if not text.strip():
            return pd.DataFrame()
        return pd.read_csv(StringIO(text), index_col=None, sep="\t")


_SPECIAL_LINE = re.compile(rb"\n[!#^]")
_BLOCK = 1 << 24


class _Reader(object):
    """Lines of a (gzipped) file, plus a fast way to take a run of table lines in one piece.

    Invariant: self.pos is at the start of a line, and the newline before it (if any) is kept
    in self.buf, so a line starting with ^ ! or # is always found as "newline + that character"
    even when it begins exactly at a block boundary.
    """

    def __init__(self, handle):
        self.handle = handle
        self.buf = b""
        self.pos = 0

    def _more(self):
        chunk = self.handle.read(_BLOCK)
        if not chunk:
            return False
        keep_from = self.pos - 1 if self.pos > 0 else 0
        self.buf = self.buf[keep_from:] + chunk
        self.pos -= keep_from
        return True

    def readline(self):
        while True:
            i = self.buf.find(b"\n", self.pos)
            if i >= 0:
                line = self.buf[self.pos:i + 1]
                self.pos = i + 1
                return line
            if not self._more():
                line = self.buf[self.pos:]
                self.pos = len(self.buf)
                return line

    def data_lines(self, keep):
        """The lines from here up to the next line starting with ^ ! or # (or the end of the file),
        as one bytes object when keep is true; otherwise they are skipped without being copied."""
        parts = []
        while True:
            m = _SPECIAL_LINE.search(self.buf, self.pos - 1 if self.pos > 0 else 0)
            if m is not None:
                end = m.start() + 1
                if keep and end > self.pos:
                    parts.append(self.buf[self.pos:end])
                self.pos = max(self.pos, end)
                return b"".join(parts)
            last_nl = self.buf.rfind(b"\n")
            if last_nl + 1 > self.pos:          # every complete line in the buffer is table data
                if keep:
                    parts.append(self.buf[self.pos:last_nl + 1])
                self.pos = last_nl + 1
            if not self._more():                # end of file: the rest is table data too
                if keep:
                    parts.append(self.buf[self.pos:])
                self.pos = len(self.buf)
                return b"".join(parts)


def _open(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def _text_lines(raw):
    """Decode one line the way GEOparse's text-mode file does (utf-8, errors ignored, and a lone
    carriage return also ends a line)."""
    text = raw.decode("utf-8", errors="ignore")
    if "\r" in text:
        return text.replace("\r\n", "\n").replace("\r", "\n").splitlines(True)
    return [text]


def read_geo_soft(path, tables=True, sample_tables=None):
    """Parse GSExxx_family.soft(.gz) into the objects GEOparse gives the notebook.

    tables=False skips every data table (fast). sample_tables=N keeps the platform tables and
    only the first N sample tables of each platform (what a person looks at before converting).
    Tables are found without looking at their lines one by one, so skipping them is quick.
    """
    gsms, gpls = OrderedDict(), OrderedDict()
    series = {"name": None, "meta": {}}
    kept_per_platform = defaultdict(int)
    cur = None   # [kind, name, meta_lines, column_lines, table_parts, has_body, keep_table]

    def finish(entry):
        if entry is None:
            return
        kind, name, meta_lines, col_lines, parts, _body, keep = entry
        if kind == "SERIES":
            series["name"], series["meta"] = name, _metadata(meta_lines)
        elif kind in ("SAMPLE", "PLATFORM"):
            table = (_table(b"".join(parts)) if parts else pd.DataFrame()) if tables and keep else None
            cls = GSM if kind == "SAMPLE" else GPL
            (gsms if kind == "SAMPLE" else gpls)[name] = cls(name, _metadata(meta_lines), _columns(col_lines), table)

    def wants_table(entry):
        kind, meta_lines = entry[0], entry[2]
        if not tables or kind not in ("SAMPLE", "PLATFORM"):
            return False
        if kind == "PLATFORM" or sample_tables is None:
            return True
        platform = "?"
        for line in meta_lines:
            if line.startswith("!Sample_platform_id"):
                platform = _parse_entry(line.rstrip())[1]
                break
        if kept_per_platform[platform] >= sample_tables:
            return False
        kept_per_platform[platform] += 1
        return True

    with _open(path) as handle:
        reader = _Reader(handle)
        while True:
            raw = reader.readline()
            if not raw:
                break
            first = raw[:1]
            if first not in (b"^", b"!", b"#"):
                if cur is None:
                    continue
                cur[5] = True
                if cur[6] is None:
                    cur[6] = wants_table(cur)
                block = raw + reader.data_lines(cur[6])
                if cur[6]:
                    cur[4].append(block)
                continue
            for line in _text_lines(raw):
                c = line[:1]
                if c == "^":
                    if cur is not None and not cur[5]:
                        continue      # GEOparse uses the first of consecutive header lines
                    finish(cur)
                    kind, name = _parse_entry(line)
                    cur = [kind.upper(), name, [], [], [], False, None]
                elif cur is None:
                    continue
                elif c == "!":
                    cur[5] = True
                    cur[2].append(line)
                elif c == "#":
                    cur[5] = True
                    cur[3].append(line)
                elif line.strip():    # a piece of a line split at a lone carriage return
                    cur[5] = True
                    if cur[6] is None:
                        cur[6] = wants_table(cur)
                    if cur[6]:
                        cur[4].append(line.encode("utf-8"))
        finish(cur)
    return GSE(series["name"], series["meta"], gsms, gpls)
