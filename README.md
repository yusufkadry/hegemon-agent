# hegemon-agent

Finds a public dataset (GEO or ArrayExpress) and turns it into a native Hegemon dataset the way the lab's notebooks do. It then publishes it to your Hegemon site and proves it works through Hegemon's own pages.

## Ten datasets on one page (one shared key)

```
hegemon batch "regulatory T cells in Crohn's disease" --key tregcd
```

It searches, ranks, then processes datasets one at a time until 10 are published, and gives every one of them
`key= tregcd`. Hegemon's own `explore.php` lists every dataset carrying a key, so all 10 appear in the dropdown at
a single URL:

```
http://hegemon.ucsd.edu/~yusuf.kadry1/Hegemon/explore.php?key=tregcd
```

The homepage gets **one** link for the whole set, not ten. Nothing in Hegemon changes: the `key=` field already
works this way.

To process the next 10 later, run the same command again:

```
hegemon batch --key tregcd            # reuses the saved query, adds 10 more
hegemon batch --key tregcd --count 5  # or add 5
```

It skips what is already published under that key (it reads `explore.conf`, so that is always accurate) and works
down the ranking for new ones. A different set of datasets goes under a different key:

```
hegemon batch "ulcerative colitis biopsies" --key uc
```

Each dataset is published the moment it passes, so a batch that dies halfway still leaves a usable page. One
dataset failing never stops the batch; the summary at the end lists what did not make it and why. An OpenAI key,
quota or model problem does stop it, because every later dataset would fail the same way.

Options: `--count N` (default 10), `--label "text for the homepage link"`, `--retry-failed` (also retry datasets
that failed in an earlier run for a dataset-specific reason), `--max-samples N` (leave out studies with more than N
samples; default **250** for `batch` and the search menu, `0` for no limit). The lab's `jupyter_gse_processing`
takes several minutes on a study of 1,000 samples, so the default keeps each dataset to a few minutes.

A record of each batch is kept in `<runs>/batches/<key>.json`: every accession tried, its outcome, and its build
folder.

## What gets published: enough samples, and metadata you can group by

Hegemon is for comparing groups of samples. Two kinds of dataset look fine but are useless there, and the agent
refuses both:

- **Too few samples.** Every point on a Hegemon plot is a sample. `--min-samples` (default **20** for `batch` and
  the search menu) leaves smaller studies out, using GEO's own sample count, before anything is downloaded. For
  ArrayExpress the SDRF decides it, before the data files are fetched. `build` of a named accession has no minimum
  unless you give one.
- **Metadata that groups nothing.** Every GEO sample has a `title` (different for each sample) and an `organism`
  (the same for all), so "has metadata" is always true and means nothing. A field only counts if it splits the
  samples into at least two groups with two or more samples each: disease vs control, treatment, tissue,
  genotype, time, response. IDs, titles, organism, platform and file names never count.

That rule runs twice:

1. **When the files are built.** Every metadata field the source holds (each GEO characteristic, each SDRF
   `Characteristics[...]`/`Factor Value[...]` column) must be in `survival.txt`, named the way the notebooks name
   it. If the AI's script dropped any, it is told which ones and puts them back. If nothing in the file groups the
   samples, the AI is asked once to look for a study-supplied sample table (clinical data, sample sheet) and join
   it. If there is none, the build stops with `NO_GROUPING_METADATA`. It is not published.
2. **When it is published.** Through Hegemon's own `getpatientinfojson` / `getpatientdatajson`, every metadata
   field is read back from the live site and the same rule is applied. The groups are printed:

   ```
   OK Hegemon can group samples by 'c disease_ch1': Crohn's 24 / control 18
   ```

   If no field groups the samples, publishing is refused and both site files are restored.

## Taking datasets off the site

```
hegemon diagnose --key tregibd                   # what is wrong with each one
hegemon unpublish GSE327665-GPL24676 GSE123-GPL1  # take off specific ones (section names from diagnose)
hegemon unpublish --key tregibd                   # take off everything under a key
```

It shows what it will remove and asks first (`--yes` to skip). It removes the `explore.conf` sections and, when
no dataset is left under a key, that key's homepage link. The previous `explore.conf` and `index.html` are saved
in `.hegemon-agent-backups/`, and both are restored if `explore.php` stops loading. The data files are not
deleted. After taking a set off, run the batch again and it rebuilds that key from the next good datasets.

## Checking what is already on your site

```
hegemon diagnose             # every dataset in explore.conf
hegemon diagnose --key tregcd
```

For each dataset it checks that the files exist, that Hegemon lists it, that a gene loads, and that the sample
metadata comes back with real values. It names what is wrong with each broken one. Use this on datasets published
by an earlier version of the agent: a dataset that opens but cannot plot shows up here as

```
none of the 18 sample columns in expr.txt appear in survival.txt's ArrayId column,
so every metadata value is blank in Hegemon
```

Rebuild one with `hegemon build <accession>`. Its old entry stays in `explore.conf` until you delete that section.

## What is native and what the AI does

The lab's two notebooks (in `blueprints/`) do the same thing every time:

1. Hand-edited, dataset-specific cells build `PREFIX-expr.txt`, `PREFIX-survival.txt` and `PREFIX-ih.txt`.
2. A fixed idx cell writes `PREFIX-idx.txt`.
3. A fixed cell runs `bash .../jupyter_gse_processing <folder name>`, which makes thr/info/vinfo/bv.

The agent keeps that split:

| Step | Who does it |
|---|---|
| Download files | Fixed code first: it downloads the smallest set that reproduces the notebook's inputs (see below), with no AI. If that is not enough, the AI sees every file the study offers and asks for more (see "Any dataset, any site"). |
| Dataset-specific cells | A GEO microarray whose sample tables are in the SOFT file is the notebook's fixed microarray path, so those cells run as they are, **with no AI call**. If they fail or leave nothing to group by, the AI takes over. Everything else: the AI writes them. The two lab notebooks are in its prompt as the blueprint, and the notebook's fixed cells are available to it as helpers (`hegemon_nb/helpers.py`). |
| Check the three files | Fixed rules (`hegemon_nb/checks.py`): every expr column has exactly one survival row and one ih row, GEO columns are GSMs, nothing is invented, values are numbers on a log scale, ProbeIDs are unique. If a check fails, the exact errors go back to the AI (up to 5 attempts). |
| idx | The notebook cell, copied as-is (`hegemon_nb/native.py`). |
| thr/info/vinfo/bv | The lab's real `jupyter_gse_processing`, called exactly like the notebooks call it. |
| Publish | Fixed code. It adds one `[section]` to `explore.conf` and one link to `index.html`. It then checks through Hegemon's own pages that the dataset is listed with the right n, a gene's row loads, its thresholds load, **and the sample metadata comes back with real per-sample values**. If anything fails, both files are restored byte-for-byte. |

## What it downloads, and why that matters for speed

Doing this by hand is fast because you look at the GEO page, pick the one file with the processed matrix, and
`wget` that. The agent now does the same thing. For a GEO series it takes the first of these that works:

1. The sample tables already inside the family SOFT file (the notebook's microarray path). Nothing else is fetched.
2. GEO's own NCBI-computed counts for a sequencing series — a few MB, and the notebook's "single file" path.
3. Series-level processed files that are not archives (a counts matrix, a series matrix).
4. Each sample's own supplementary file, fetched individually into `GSE###_RAW/`. That is the same content as
   `GSE###_RAW.tar` and the same folder the notebook gets by untarring it, so **the bundle is never downloaded**
   when the per-sample files exist. Raw formats (BAM, FASTQ, CEL, IDAT) are skipped.
5. An archive bundle, only when nothing above yields data.

Files over 16 MB are fetched as 4 byte ranges at once (repositories limit each connection, not the total), and up
to 4 files download at the same time. A range whose connection drops resumes where it stopped; a server that does
not accept ranges gets one normal download.

Some studies offer the same data several times: raw counts, normalized counts, TPM, one file per tissue or
cohort. When several such files add up to more than 100 MB (`HEGEMON_PICK_ABOVE_MB`), the agent reads the first
lines of each one with a single small request and the fast model picks the ones the conversion needs, usually
raw counts for every tissue and cohort, never the same values twice in different units. Only those are
downloaded. The rest stay in the study's file list, so the conversion can still ask for them. Per-sample files
of several kinds (for example `GSM#_counts` and `GSM#_fpkm`) are picked the same way. For E-MTAB-14509 that means
~420 MB of raw counts instead of all ~800 MB.

Every downloaded file is checked to be what its name says (a `.gz` must be gzip data, and no data file may be a web
page). A server that is busy or refusing requests often answers with an error page instead of the file; that page
is never saved as data. The download is retried, and anything that failed while several files were downloading at
once is tried again on its own, because some servers (NCBI's download page among them) refuse parallel requests.
A source is only used when all of its files arrived: GEO's computed counts without their gene annotation file are
dropped and the series' own files are used instead. If nothing arrives, the build ends with `DOWNLOAD_FAILED`,
never with "this study has no data".

Sizes are shown before anything large is fetched. Anything refused by `HEGEMON_MAX_FILE_GB` /
`HEGEMON_MAX_TOTAL_GB` is reported as `FILES_TOO_LARGE`, never as a dataset with no data.

For ArrayExpress it reads the SDRF first, the way you would: the SDRF's derived/processed data-file columns say
which files hold the processed data, and **only those are fetched**. Each is checked to exist before any large
download. A named file that is not served loose is looked for inside the study's processed archives (older
studies package them that way); raw archives and files nobody named are never downloaded. If the named file is
nowhere, the build ends in seconds with `NOT_HOSTED` instead of after downloading everything else first. Only when
the SDRF names no processed file does it fall back to the study's non-raw files, smallest first.

Single-cell studies are skipped (`SINGLE_CELL`) before their data is downloaded: by title in the search, and by
the SOFT file, IDF/SDRF or per-cell files (barcodes, `matrix.mtx`, `.h5ad`, `.loom`) after that. The agent has no
single-cell blueprint yet; adding the lab's `single-cell-process.ipynb` to `blueprints/` is the next step for them.

## How long each dataset took

Every build ends with one line saying where the time went, and the same numbers go into `build.json`:

```
AI use: 0 call(s), 0 input / 0 output tokens (built-in notebook conversion)
time: download … | reading files … | conversion (built-in, no AI) … | idx + lab script … | publish … | total …
```

The batch summary shows each dataset's total. The AI writes its first attempt at medium reasoning, which is several
times faster than high; if that attempt fails the checks, the repairs use high.

Big studies are where the time goes, so the agent's own steps are built for them. Measured on a 434-sample ×
54,675-probe GPL570 series (the size of GSE201827) on a 2-core machine:

| Step | Before | Now |
|---|---|---|
| Reading the SOFT file's metadata (download step) | 25 s | 6 s |
| "Reading the files" step (it read every table a second time) | 55 s | under 0.1 s |
| Reading every sample table (the conversion) | 55 s | 29 s |
| Merging 434 sample tables into one matrix | 7 s | 1 s |
| Writing the expr file | 25 s | 18 s (split over the CPU cores; faster with more cores) |
| Checking the expr file | 20 s | 7 s |

Every one of these gives the same output as before, byte for byte: the SOFT reader is compared with GEOparse
itself (including files with CRLF lines, stray carriage returns, trailing tabs and invalid UTF-8), and the built
GEO files are compared with the lab notebook's own output. `jupyter_gse_processing` is the lab's script and is
called unchanged; on a study this size it takes as long as it takes when the notebook runs it.

## The page opens on a plot

Hegemon's page draws the plot and the **Select Patient Information** dropdown only when both Gene A and Gene B
exist in the selected dataset. Its defaults (TYROBP / FCER1G) are missing from many datasets, and then the page
shows nothing at all, which looks exactly like "no metadata". So every published link names two genes that exist
in **every** dataset under its key, with `cmd=explore` so the page draws the plot on load:

```
http://hegemon.ucsd.edu/~yusuf.kadry1/Hegemon/explore.php?key=tregibd&A=FOXP3&B=IL2RA&cmd=explore
```

The genes are the topic's markers when every dataset has them (the ranking suggests them, e.g. FOXP3 and IL2RA
for regulatory T cells), otherwise housekeeping genes such as ACTB and GAPDH. `--genes A,B` picks them yourself.
Publishing proves the page draws the plot and the dropdown with them before the link is written, and the key's
homepage link is updated whenever a dataset joins the key.

For links published before this existed:

```
hegemon relink                 # every key on the site
hegemon relink --key tregibd   # one key
```

## Any dataset, any site

Nothing after the first download step depends on where a study comes from. Every study gets a list of all the
files it offers:

- GEO: the series folder, each sample's own files, and GEO's computed counts
- ArrayExpress: the study folder (BioStudies and FTP)
- any other site: the links on the study's web page

For GEO and ArrayExpress a fixed plan downloads the cheapest files that usually work, with no AI. When those are
not enough, the AI doesn't stop. It sees the list, with sizes and anything that failed to download, and replies
`FETCH: <name>`, a name pattern (`GSE123_RAW/*_counts.txt.gz`) or a public URL. Those files are downloaded and
checked, archives are unpacked, and the AI gets the evidence again (up to 3 rounds per build). A failed download,
a file the repository lists but doesn't serve, or data in a layout nobody planned for all go through that same
loop. "A file is missing" is never a reason to stop: the AI is told to fetch it.

For a study on any other website, give its page or its files:

```
hegemon build MYSTUDY --url https://some-site.org/study/123          # a web page: its links become the list
hegemon build MYSTUDY --url https://host/path/counts.tsv.gz --url https://host/path/samples.tsv
hegemon build MYSTUDY --source ~/somewhere/with/the/files
```

Samples are matched by what fits: every sample with both expression values and metadata is kept. Data columns
with no metadata, and samples with no data, are left out and counted in the output ("matched 376; left out 26
columns without metadata"). A study fails on matching only when fewer than two samples match. The AI may only
fetch public addresses; links to private or local networks are refused.

The search itself covers GEO and ArrayExpress; a study from anywhere else comes in through `--url`.

## The metadata check

A dataset that opens in Hegemon but cannot plot anything is almost always one thing: Hegemon joins the `ArrayId`
column of `survival.txt` to the sample columns of `expr.txt`, and when they do not match it shows the dataset with
every metadata value blank. Nothing is groupable and no comparison can be made.

The agent refuses to publish that. Before a dataset goes live it asks Hegemon itself, through the same endpoints
the lab's own `HegemonUtil.py` uses:

- `getpatientinfojson` — does Hegemon list the metadata fields?
- `getpatientdatajson` — for each field, do real values come back for the samples?

If every value is blank, publishing is refused and `explore.conf` and `index.html` are restored byte-for-byte. The
build folder is kept so the cause can be found. `hegemon diagnose` runs the same checks on datasets that are
already published.

## Setup (nothing to install)

It uses the Python 3.8 venv already in `~/hegemon-agent-v0/.venv` (pandas and numpy), and the OpenAI key already saved at `~/.config/hegemon-agent/openai.key`.

On your Mac:

```
scp hegemon-agent.zip yusuf.kadry1@hegemon.ucsd.edu:~
ssh yusuf.kadry1@hegemon.ucsd.edu
```

On hegemon:

```
unzip -o ~/hegemon-agent.zip -d ~
~/hegemon-agent/hegemon check
```

`check` makes one real test request to OpenAI per model. If the model name is wrong for your key, set it, for example `export HEGEMON_MODEL=gpt-5.2`.

## First run: prove it against the lab's own files

```
~/hegemon-agent/hegemon build E-MTAB-7604 --no-publish
~/hegemon-agent/hegemon compare <build folder>/E-MTAB-7604  <folder where the lab ran the E-MTAB-7604 notebook>
```

Use the same pair of commands for GSE51984 (its prefix is `GSE51984-<GPL>`). The compare output lists every difference in samples, probes, values, Names and metadata columns.

Expected, deliberate differences:

- **GEO metadata.** The agent never splits a sample title on `:`. It splits `key: value` fields on the first `:` only, so `treatment: drug: 5 mg` keeps its value instead of becoming `c treatmen_ch1 = N/A`.
- **Samples split across several count files.** The notebooks' `df[arr] += df1[arr]` adds files by row position after an outer merge has re-sorted the rows. That can add counts of different genes. The agent adds them gene by gene. The local tests show this happening on a small example (`tests/run_local_tests.py`, ArrayExpress test).

## Daily use

```
~/hegemon-agent/hegemon                                      # search, pick one from the menu, build, publish
~/hegemon-agent/hegemon batch "your topic" --key yourkey     # 10 datasets on one page
~/hegemon-agent/hegemon batch --key yourkey                  # the next 10 under the same key
~/hegemon-agent/hegemon build GSE12345                       # one accession
~/hegemon-agent/hegemon diagnose                             # check what is already on the site
~/hegemon-agent/hegemon unpublish --key yourkey              # take a set off the site (files kept)
~/hegemon-agent/hegemon relink                               # make every key's link open on a plot
~/hegemon-agent/hegemon publish <build folder>               # publish a build made with --no-publish
```

Build options:

| Option | What it does |
|---|---|
| `--platform GPL570` | Pick the platform when a GEO series has several. |
| `--source DIR` | Use files you already downloaded. |
| `--script FILE` | Start from your own conversion script. The AI only steps in if it fails. |
| `--yes` | Skip the size and publish questions. |
| `--min-samples N` / `--max-samples N` | Leave out studies smaller / bigger than N samples (`0` = no limit). |

Everything from a run is kept in `<runs>/builds/<accession>-<time>/`:

- `source/`: the downloaded files
- `evidence.txt`: exactly what the AI was shown
- `ai-reply-N.txt`: each AI reply
- `convert.py`: the script that worked
- `logs/`: output of every attempt and of the lab script
- `build.json`: the receipt, with file hashes and the published URL
- `<PREFIX>/`: the dataset folder

Runs go to `/booleanfs2/sahoo/Data/BooleanLab/<you>/agent_runs` when that folder exists, otherwise to `~/hegemon-agent-runs`.

## When a build fails

The last line says why, and it says whether the cause is the dataset or not.

These are about the dataset:

- `TOO_FEW_SAMPLES`: fewer samples than `--min-samples`.
- `TOO_MANY_SAMPLES`: more samples than `--max-samples` (only when you set it; not remembered as a failure, so a
  later run without the limit tries it again).
- `NOT_HOSTED`: the file that holds the data is listed but can't be downloaded, and no other file holds it
  (decided by the AI after looking at everything the study offers).
- `SINGLE_CELL`: a single-cell study; skipped until the single-cell blueprint is added.
- `NO_GROUPING_METADATA`: the samples carry no annotation that splits them into groups (only IDs, titles, or
  fields that are the same for everyone), and no study-supplied sample table was found.
- `NO_PROCESSED_DATA`: the study posted only raw reads or images, so there is nothing to convert.
- `NO_SAMPLE_MATCH`: fewer than two data columns can be matched to samples (partial matches are converted).
- `NOT_EXPRESSION`: not a gene-by-sample table.

These are not about the dataset, and I want to hear about them:

- `MISSING_FILES`: the AI still says a file is missing after being told to fetch it. Rare now.
- `FILES_TOO_LARGE`: the data exists but the size limits refused it. Raise `HEGEMON_MAX_FILE_GB` and
  `HEGEMON_MAX_TOTAL_GB`, or pass `--yes`, and run it again.
- `AGENT_GAVE_UP`: this agent could not write a working conversion. Its limit, not the dataset's.
- `AI_PROVIDER_*`: an OpenAI key, quota or model problem.
- `LAB_SCRIPT_FAILED` or `SUPPORT_FILES_BAD`: the lab script's output.
- `DOWNLOAD_FAILED`: the metadata file itself could not be downloaded. When a data file fails, the AI is told and
  can use another one.

If a publish fails with "does not list the dataset" or "explore.php crashed", the web server can't read the files. Check the permissions of the runs folder and every folder above it.

## Settings (environment variables)

| Variable | Default |
|---|---|
| `HEGEMON_MODEL` / `HEGEMON_FAST_MODEL` | `gpt-5.6-sol` / `gpt-5.6-luna` (the models V6 used) |
| `HEGEMON_REASONING` | `medium` for the first attempt, `high` for repairs. Set it to force one level for every call. |
| `HEGEMON_BUILTIN=0` | always ask the AI, even for a GEO microarray the notebook's fixed cells can convert |
| `HEGEMON_DOWNLOAD_PARTS` / `HEGEMON_DOWNLOAD_WORKERS` | `4` / `4`: byte ranges per large file, files at once. `1` turns each off. |
| `HEGEMON_PARALLEL_WRITE_CELLS` | `4000000`: expr files with more values than this are written by up to 4 processes |
| `HEGEMON_ATTEMPTS` | `5` |
| `HEGEMON_AI_TIMEOUT` | `900` seconds per AI request |
| `HEGEMON_MAX_FILE_GB` / `HEGEMON_MAX_TOTAL_GB` | `5` / `20`. It asks before downloading more; nothing is skipped silently. |
| `HEGEMON_SITE_ROOT` / `HEGEMON_SITE_URL` | `~/public_html/Hegemon` / `http://hegemon.ucsd.edu/~<you>/Hegemon` |
| `HEGEMON_BOOLEANLAB_SCRIPT` | the path from the notebooks |
| `HEGEMON_GENOME_DIR` | the path from the GSE notebook |
| `HEGEMON_RUNS` | see above |
| `HEGEMON_LIVE_CHECK=list` | only check that the dataset is listed, if Hegemon's JSON endpoints are broken on the server |

## What has and hasn't been tested

`tests/run_local_tests.py` passes 32/32 on Python 3.8.20 with pandas 2.0.3 (the server's versions). What those tests ran for real:

- The SOFT reader gives the same objects as GEOparse.
- A GEO microarray build gives expr, ih and idx files byte-identical to what the lab's GEO notebook cells produce on the same files.
- The ArrayExpress build matches the E-MTAB notebook's output for every sample except the split-file case above.
- The notebook's idx code on Hegemon's real sample dataset (48,702 probes × 256 samples) gives the same ProbeIDs, pointers and Names as the idx Hegemon ships.
- Publishing through Hegemon's real PHP (served locally with a PHP 8 shim) works, and rollback restores both files exactly.
- The AI repair loop and the provider-error reporting work against a local mock of the OpenAI endpoint.
- Three datasets published under one key all appear at `explore.php?key=...` with one homepage link, a duplicate
  is refused, and a second key stays independent.
- A dataset whose `survival.txt` ArrayIds do not match its expr columns is refused by the metadata check, with both
  site files restored; `diagnose` names that same dataset when it is already published.
- A batch publishes what works under one shared key, carries on past a dataset that fails to download, records
  every outcome, and on a second run skips what is already published.
- Against a local stand-in for GEO: when its computed counts exist only those are fetched and the bundle is left
  alone; without them the per-sample files are fetched individually into `GSE###_RAW/` and the bundle is still left
  alone, with raw files skipped; the bundle is used only when nothing else yields data; and a size limit produces
  `FILES_TOO_LARGE` rather than looking like an empty dataset.
- A dataset whose only metadata is `title` and `organism` is refused at publishing, with both site files
  restored; so is one whose survival IDs do not match its expression columns.
- When the AI's script drops metadata the SOFT file has, the next attempt is told exactly which fields
  (`c tissue_ch1`, ...) and the build ends with them present.
- A study whose samples have nothing that groups them ends as `NO_GROUPING_METADATA` after one request for a
  study-supplied sample table, and is never published.
- Studies under `--min-samples` are left out of a batch from GEO's own counts without any download, and a GEO
  download stops after the SOFT file when the study is too small.
- `unpublish` removes a single section (the key's homepage link stays while others remain) or a whole key (link
  removed, other keys untouched, backup kept, `explore.php` still renders).
- ArrayExpress: the SDRF-named file is found inside the study's processed archive while the raw archive and a
  large file nobody named are never downloaded; a named file that is not served anywhere ends as `NOT_HOSTED`
  without the large file being fetched first; a named file over the size limit is `FILES_TOO_LARGE`.
- Single-cell studies are skipped from the GEO SOFT file alone, from per-cell files in an ArrayExpress folder (the
  `.h5ad` is never fetched), and by title before ranking in a batch.
- A characteristic named `status` (disease status) counts as metadata.
- Published links name two genes every dataset under the key has (topic genes when all datasets have them,
  otherwise ACTB/GAPDH), Hegemon draws the plot and the metadata dropdown with them, and `relink` rewrites an
  older link and is safe to run twice.
- The ArrayExpress path requires the SDRF, fetches every per-sample file the SDRF names (including files in
  subfolders such as `counts/`), keeps that folder structure, refuses to build when one is missing, and rejects
  file paths that try to escape the download folder.
- A GEO microarray is converted by the notebook's fixed cells with no AI call (the AI endpoint was dead during the
  test), and its expr, ih and idx files are byte-identical to the lab notebook's. When those cells cannot do a
  dataset (a platform table without a gene-symbol column), the AI takes over with a fresh medium-reasoning
  request, and a repair after a failed check uses high.
- A server that answers with an error page: the page is never saved (a retry gets the file); a server that refuses
  parallel requests still delivers all 6 files, one at a time; GEO counts whose annotation file never arrives are
  dropped (nothing half-downloaded is left) and the per-sample files are used; with nothing else to use, the
  build is a download failure, not "no processed data".
- The parallel expr writer, the side-by-side sample matrix and the quick number check give exactly what plain
  pandas gives, on awkward inputs (NaN, inf, -0.0, quotes, integer columns, reordered or missing probes, duplicate
  IDs, true/false cells, words, wrong field counts).
- `--max-samples` leaves big studies out of a batch before anything downloads, and a build over the limit stops
  before converting.
- A study on a website the agent has never seen, given only its web page: the page's 3 links become the list, the
  AI says a file is missing and is told to fetch it, asks for the counts and the sample sheet, the raw reads are
  never downloaded, a data column with no metadata is left out, and the dataset builds. Name patterns fetch
  several per-sample files at once, and a URL pointing at a private address is refused.
- A GEO download that fails, and an ArrayExpress file the SDRF names but the study doesn't serve, are handed to
  the AI along with everything else the study offers, instead of ending the build.
- When an SDRF names skin raw counts, skin normalized counts and blood raw counts, the first lines of each are
  shown to the picker and only the chosen raw-count files are downloaded; the other stays in the file list. GEO
  per-sample files of two kinds are narrowed to one kind. The SDRF or SOFT file decides `--max-samples` before any
  data file downloads, and a batch leaves out studies over 250 samples by default.
- A 20 MB file downloaded as 4 ranges is byte-identical to the original, including when one range's connection is
  cut half way; a server without range support gets a normal single download.

Not tested locally:

- the lab's `jupyter_gse_processing` (a FAKE stand-in wrote thr/info/vinfo/bv)
- real OpenAI calls
- real downloads from GEO and ArrayExpress (the sandbox can't reach them)

The first real build on the server is the real test of those three.

Not supported yet: single-cell datasets. The lab's `single-cell-process.ipynb` needs to be added to `blueprints/` first.

Also not verified here: the exact column names in a real ArrayExpress SDRF, because this sandbox cannot reach EBI.
The code reads the MAGE-TAB derived/processed data-file columns and falls back to the study's full file list, so a
column name I have not seen means extra files are downloaded, not that files are missed.
