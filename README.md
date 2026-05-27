# Aluminion

**Automated pipeline for bacterial whole-genome sequencing from Oxford Nanopore reads.**

Aluminion takes a `fastq_pass` folder containing demultiplexed Oxford Nanopore
reads from a MinKNOW run and produces polished assemblies, taxonomic classification, 
AMR profiles, mobile-genetic-element detection (plasmids, integrons, prophages, IS elements), 
an interactive HTML report, and two **cumulative lab databases** that grow across sequencing runs.

Target organisms: Enterobacteriaceae. The assembly, AMR, and annotation modules
also work on any bacterial species with a published MLST scheme.

---

## Pipeline overview

```mermaid
flowchart TB
    IN(["MinKNOW run<br/>list_seq.tsv"])
    S1["<b>1 · QC and filtering</b><br/>NanoPlot · Chopper"]
    S2["<b>2 · Assembly</b><br/>Flye · dorado polish<br/>deconcat · dnaapler"]
    S3["<b>3 · Assembly QC</b><br/>Kraken2 - QUAST · Bandage"]
    S4["<b>4 · Annotation, typing and AMR</b><br/>Bakta · Abricate · MOB-suite<br/>GAMBIT · MLST · Kleborate · ECTyper"]
    S5["<b>5 · Mobile genetic elements</b><br/>Copla · Phastest · Integron_Finder · ISfinder BLASTn"]
    S6["<b>6 · Consolidation</b><br/>parser.py → aluminion_reporter.py → lab DB"]
    S7["<b>7 · Cross-run MGE alerts</b><br/>mge_repository → mge_alerts → alerts_reporter"]
    OUT(["Aluminion_Report.html · data_seq.tsv · data_analysis.tsv<br/>alerts.tsv · Alerts_Report.html · repository/"])

    IN --> S1 --> S2 --> S3 --> S4 --> S5 --> S6 --> S7 --> OUT
```

| Stage | Tools                                                            | Conda env(s)                                |
|------:|------------------------------------------------------------------|---------------------------------------------|
| 1     | NanoPlot, Chopper                                                | `aluminion_reads`                           |
| 2     | Kraken2, Flye, Dorado polish, deconcat, dnaapler, QUAST, Bandage | `aluminion_assembly`, `aluminion_circlator` |
| 3     | Bakta, Abricate, MOB-suite, GAMBIT, MLST, Kleborate, ECTyper     | `aluminion_annot`, `aluminion_kleborate`    |
| 4     | Copla, Phastest, Integron_Finder, ISfinder BLASTn, Abricate-VFDB | `aluminion_annot`, `aluminion_integron`     |
| 5     | parser.py, aluminion_reporter.py, lab_db_updater.py              | `aluminion_annot`                           |
| 6     | mge_repository.py, mge_alerts.py, alerts_reporter.py             | `aluminion_annot`                           |

---

## System requirements

| Resource | Minimum               | Recommended                                  |
|----------|-----------------------|----------------------------------------------|
| OS       | Linux (Ubuntu 20.04+) | Ubuntu 22.04 LTS                             |
| CPU      | 16 cores              | 32+ cores                                    |
| RAM      | 64 GB                 | 128 GB (allows the Kraken2 DB in `/dev/shm`) |
| Disk     | 500 GB free           | 1 TB+ free                                   |
| GPU      | —                     | NVIDIA GPU ≥16 GB VRAM (for Dorado polish)   |

macOS is not supported (Docker networking requirements for Phastest).
Windows is community-tested via WSL2 but not officially supported.

> With less than 128 GB of RAM, comment out the `/dev/shm` copy in
> `aluminion.sh` and point `--db` directly to disk.
---

## Installation

### Automated installer (recommended)

```bash
git clone https://github.com/Aluminio-visto/aluminion.git
cd aluminion
chmod +x aluminion.sh install.sh
./install.sh -b /path/to/Databases
```

`install.sh` creates the conda environments, pulls the Docker images, optionally
downloads the databases, and installs the `aluminion` command in `~/.local/bin`.

| Flag           | Effect                          |
|----------------|---------------------------------|
| `--skip-envs`  | Skip conda environment creation |
| `--skip-docker`| Skip Docker image pulls         |
| `--skip-dbs`   | Skip database downloads         |

### Manual installation

```bash
# 1 · Clone
git clone https://github.com/Aluminio-visto/aluminion.git
cd aluminion
chmod +x aluminion.sh

# 2 · Install Mambaforge (skip if already present)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh

# 3 · Create the six conda environments
for env in reads assembly circlator annot integron kleborate; do
    mamba env create -f envs/aluminion_${env}.yml
done

# 4 · Install Docker (required for MOB-suite, Copla, Phastest)
sudo apt-get install -y docker.io docker-compose-plugin
sudo usermod -aG docker $USER && newgrp docker

# 5 · Pull Docker images
docker pull kbessonov/mob_suite:3.0.3
docker pull rpalcab/copla:1.0
```

**Dorado** is a proprietary basecaller / polisher from ONT. Download the latest
Linux binary from <https://github.com/nanoporetech/dorado/releases> and place it
in your `PATH`. GPU polishing requires NVIDIA drivers and CUDA ≥12.

**Phastest** has no public Docker image. Follow <https://phastest.ca> to set up
the local docker-compose, then point Aluminion to it with
`export ALUMINION_PHASTEST_DIR=/path/to/phastest-docker`.

---

## Database setup

All databases live under a single root directory (passed with `-b`).

| Database         | Path                                | Size    | Source                                                            |
|------------------|-------------------------------------|---------|-------------------------------------------------------------------|
| Kraken2 standard | `<db>/Kraken/`                      | ~100 GB | <https://genome-idx.s3.amazonaws.com/kraken/>                     |
| GAMBIT           | `<db>/gambit/`                      | ~1 GB   | <https://github.com/jlumpe/gambit>                                |
| Bakta            | `<db>/bakta/db/`                    | ~30 GB  | `bakta_db download --output <db>/bakta --type full`               |
| ISfinder         | `<db>/ISfinder/ISfinder-nucl.fasta` | ~5 MB   | <https://www.is-finder.org/download.html>                         |
| Abricate         | (auto, in `aluminion_annot`)        | ~150 MB | `abricate-get_db --db {ncbi,resfinder,card,argannot,vfdb}`        |
| MLST (PubMLST)   | (auto, on first run)                | ~10 MB  | downloaded by the `mlst` tool itself                              |



---

## Input files

### `list_seq.tsv` — per-run sample sheet

Tab-separated, six columns. Aluminion searches for it in the following order:

1. The path passed to `-l / --list` (explicit override, wins over anything else).
2. `list_seq.tsv` already present **inside the run folder** — the preferred
   layout when each run is self-contained (see *Batch processing* below).
3. `list_seq.tsv` in the **parent working directory** (`-d`) — the single-run layout.

The first location found is used; the file is copied into the run folder if it
isn't already there. If none is found, Aluminion writes an empty template and
exits with an error.

> **Recurrent / batch analyses:** location 3 is **disabled**. When a single parent
> directory holds many runs, each run **must** carry its own `list_seq.tsv` inside
> its run folder — otherwise the sample IDs and barcodes of different runs would
> collide on a shared parent-level sheet. `aluminion_batch` enforces this by passing
> `--require-run-list` (which drops the parent fallback) and skips any run folder
> that lacks its own `list_seq.tsv`.

| Column        | Description                                                       |
|---------------|-------------------------------------------------------------------|
| `Lab_id`      | Internal lab culture ID                                           |
| `Strain`      | Strain collection code                                            |
| `ID`          | **Unique sample identifier — used as the sample name throughout** |
| `Barcode`     | Barcode number assigned by MinKNOW (`01`, `24`, …)                |
| `DNA_conc`    | DNA concentration (ng/µL), for lab QC purposes                    |
| `is_repeated` | `x` if this is a re-sequencing of a previously failed sample      |

See `examples/list_seq.tsv` for a complete example.

### MinKNOW outputs (copied automatically)

| File                  | Purpose                                              |
|-----------------------|------------------------------------------------------|
| `fastq_pass/`         | Demultiplexed FASTQ reads, one folder per barcode    |
| `final_summary_*.txt` | Run summary (instrument, flow cell, dates, duration) |
| `report_*.json`       | JSON report (pore counts, yield)                     |

Override the MinKNOW data path with `-m` or `$ALUMINION_MINKNOW_DIR`.
If already in your folder, use --

### Cumulative lab databases (auto-created on first run)

| File                | Purpose                                                        |
|---------------------|----------------------------------------------------------------|
| `data_seq.tsv`      | Sequencing QC, flow cell metadata, depth per sample (all runs) |
| `data_analysis.tsv` | Taxonomy, AMR, MGE counts per sample (all runs)                |

These are updated automatically at the end of each run. The first time, they
are created from scratch — no manual setup needed. `--init-db` forces a rebuild.

---

## Usage

### Directory layout

```
/home/user/$project_folder/      ← parent working directory (-d)
├── list_seq.tsv               ← fill this before each run
├── data_seq.tsv               ← cumulative sequencing database
├── data_analysis.tsv          ← cumulative analysis database
├── BAC_2025_NOV_25/           ← run folder, created automatically
│   ├── 01_reads/  02_filter/
│   ├── 03_assemblies/
│   ├── 04_taxonomies/  05_plasmids/
│   ├── 08_Anotacion/  09_phages/  11_integrons/
│   ├── Aluminion_Report.html
│   └── aluminion_YYYYMMDD_HHMMSS.log
└── BAC_2026_FEB_10/
    └── …
```

### Standard run

```bash
cd /home/user/$project_folder
aluminion -r BAC_2025_NOV_25 -b /path/to/Databases -t 30 -l list_seq.tsv
```

A timestamped log is written inside the run folder.

### Flags

| Flag                     | Description                                                                 | Default                       |
|--------------------------|-----------------------------------------------------------------------------|-------------------------------|
| `-r / --run`             | MinKNOW run folder name **(mandatory)**                                     | —                             |
| `-b / --db-dir`          | Path to databases root                                                      | `/home/$user/Databases`     |
| `-t / --threads`         | CPU threads                                                                 | `30`                          |
| `-l / --list`            | Path to `list_seq.tsv`                                                      | —                             |
| `-d / --dir`             | Parent working directory                                                    | `/home/$user/$project_folder` |
| `-p / --phastest-dir`    | Local Phastest docker-compose folder                                        | `~/Programs/phastest-docker`  |
| `-m / --minknow-dir`     | MinKNOW data root                                                           | `/var/lib/minknow/data`       |
| `--init-db`              | Create `data_seq.tsv` / `data_analysis.tsv` from scratch                    | —                             |
| `--resume`               | Resume an interrupted run; skip any step whose output already exists        | —                             |
| `--skip-import-from-minknow` | Skip importing data from the MinKNOW data tree; assume the run folder already contains `fastq_pass/` (and optionally the `final_summary_*.txt` / `report_*.json`). Aborts if `fastq_pass/` is missing. `--no-minknow` is a back-compat alias. | — |
| `--unique-run`           | Self-contained single-run mode. Do not create or touch the shared repository (`<parent>/repository` or `--repo`): no reads are deposited there, `is_repeated` samples are processed as fresh runs, and no cumulative artefacts are shared with other runs. | — |
| `--keep-everything`      | Keep all intermediate files. By default, after a complete run Aluminion prunes Flye staging dirs, the per-assembly BLASTn DB (`assembly.fasta.n*`), Kraken2 `.out` streams, and the per-sample reads in `01_reads/`/`02_filter/` (regenerable from the retained `fastq_pass/` with `--just-preprocessing`). | — |
| `--require-run-list`     | Disable the parent-directory `list_seq.tsv` fallback; require a per-run sheet. Set automatically by `aluminion_batch`. | — |
| `--polish-batchsize <N>` | Override `dorado polish --batchsize` (lower it on `CUDA out of memory`)     | dorado default                |
| `--skip-preprocessing`   | Skip NanoPlot + Chopper (reuse existing `01_reads/`, `02_filter/`). If those reads were pruned by the end-of-run cleanup, regenerate them with `--just-preprocessing` first (this flag aborts with a clear error if they are missing). | — |
| `--skip-kraken`          | Skip Kraken2 read-level classification                                      | —                             |
| `--skip-abr`             | Skip Abricate AMR gene screen                                               | —                             |
| `--skip-typing`          | Skip GAMBIT, MLST, Kleborate, ECTyper                                       | —                             |
| `--skip-integrons`       | Skip Integron_Finder                                                        | —                             |
| `--skip-plasmids`        | Skip Copla plasmid typing (MOB-suite always runs)                           | —                             |
| `--skip-phages`          | Skip Phastest prophage detection                                            | —                             |
| `--repo PATH`            | MGE repository directory for cross-run alerts                               | `$BASE_DIR/repository`        |
| `--init-repo`            | Create the MGE repository structure at `--repo` (idempotent)                | —                             |
| `--alert-new-priority`   | Also alert on first-occurrence MGEs carrying a priority gene                | —                             |
| `--no-alerts`            | Skip the cross-run MGE alerts step                                          | —                             |
| `--just-preprocessing`   | Stop after Stage 1. Output: `02_filter/<sample>.fastq.gz`                   | —                             |
| `--just-assembly`        | Stop after Stage 2. Output: `03_assemblies/<sample>.fasta`                  | —                             |
| `-h / --help`            | Show help and exit                                                          | —                             |

### Parallelism and performance

Several light, independent per-sample steps run as **bounded parallel pools** so they
no longer execute one sample at a time while most cores sit idle. Heavy steps that
already saturate every core (Flye assembly, dorado polish, Bakta, Integron_Finder)
are left sequential on purpose; Phastest is also kept sequential. The fastq import
copy is sequential too, since on a mechanical HDD parallel copies only cause seek
thrashing.

Pool widths are tuned via environment variables (no extra CLI flags). Defaults are
conservative for a single mechanical HDD; raise them on SSD/NVMe or many-core hosts:

| Env var                  | Step                          | Default |
|--------------------------|-------------------------------|---------|
| `ALUMINION_PAR_FILTER`   | Chopper read filtering        | `4`     |
| `ALUMINION_PAR_QC`       | NanoPlot pre/post-filter QC   | `4`     |
| `ALUMINION_PAR_ABRICATE` | Abricate AMR + virulence      | `8`     |
| `ALUMINION_PAR_MOBSUITE` | MOB-suite plasmid recon       | `2`     |
| `ALUMINION_PAR_COPLA`    | Copla plasmid typing          | `3`     |

```bash
# Example: faster filtering and AMR screening on an SSD/64-core host
ALUMINION_PAR_FILTER=8 ALUMINION_PAR_ABRICATE=16 aluminion -r RUN -b /db -l list_seq.tsv
```

`pigz` (parallel gzip) is used automatically for FASTQ compression when it is
installed in the `aluminion_reads` environment, with a transparent fallback to
plain `gzip`.

### Resuming an interrupted run

```bash
aluminion -r BAC_2025_NOV_25 -b /path/to/Databases -l list_seq.tsv --resume
```

`--resume` checks every step's expected output and silently skips completed
work. It can be combined with any `--skip-*` flag.

### Batch processing — many runs in one go

When a project contains many run folders under a single parent, use
`aluminion_batch` to walk a list of runs sequentially. **Each run folder must
contain its own `list_seq.tsv`** (the parent-level fallback is disabled in batch
mode). By default the wrapper imports `fastq_pass/` + reports from the MinKNOW
data tree for each run, exactly like a direct `aluminion` call; pass
`--skip-import-from-minknow` if the run folders already hold their `fastq_pass/`.

```
/seqs/KLEBIRE/
├── runs.tsv                 ← chronological list of run names
├── 2025_NOV_25_BAC/
│   ├── list_seq.tsv         ← mandatory; fastq_pass/ imported from MinKNOW by default
├── 2025_DEC_03_BAC/
│   ├── list_seq.tsv  [fastq_pass/  final_summary_*.txt  report_*.json if self-contained]
└── 2026_FEB_10_BAC/
    └── …
```

`runs.tsv` is a plain text file, one run name per line; `#` introduces a
comment, blank lines are ignored:

```text
# Klebsiella surveillance — 2025-26 season
2025_NOV_25_BAC
2025_DEC_03_BAC
2026_FEB_10_BAC
```

Then:

```bash
aluminion_batch --runlist runs.tsv -d /seqs/KLEBIRE -b /path/to/Databases -t 30 \
    --init-repo -- --resume --skip-kraken
```

`--init-repo` is recognised directly (it is forwarded to the first launched run;
the MGE repository is created automatically on the first run regardless).

> **The standalone `--` is real syntax, not a typo.** The wrapper only understands
> its own options (`--runlist`, `-d`, `-b`, `-t`, `--skip-import-from-minknow`,
> `--init-repo`, `--force`). Every flag meant for `aluminion` itself — `--resume`,
> `--skip-kraken`, `--polish-batchsize`, etc. — must come **after a bare `--`
> separator**, which marks "everything past here is forwarded verbatim to each
> `aluminion` invocation." So in `... -t 30 -- --resume --skip-kraken`, the wrapper
> consumes `-t 30`, then forwards `--resume --skip-kraken` unchanged. Passing
> `--resume` *without* the leading `--` makes the wrapper reject it as an unknown
> argument (`[ERROR] Unknown argument: --resume`).

Per-run logic:

| Condition                                          | Action                                                            |
|----------------------------------------------------|-------------------------------------------------------------------|
| `<run>/` missing or no `list_seq.tsv`              | Skip with a warning (each run needs its own sample sheet).        |
| `<run>/Aluminion_Report.html` exists               | Skip (run already analysed). Override with `--force`.             |
| Default (no `--skip-import-from-minknow`)          | Invoke `aluminion -r <run> -d <parent> --require-run-list <forwarded>`; reads imported from MinKNOW. |
| `--skip-import-from-minknow` set                   | Require a non-empty `<run>/fastq_pass/`; invoke with `--skip-import-from-minknow` added. Skip with a warning if reads are missing. |

A summary table (`processed` / `already done` / `skipped` / `failed`) is
printed at the end. The wrapper returns a non-zero exit code if any run
failed, so it composes cleanly inside CI or a parent script.

### Running only the consolidation stage

If assemblies and annotations are already complete, run the Python parsers
directly:

```bash
conda activate aluminion_annot
python3 scripts/parser.py -i /path/to/run/
python3 scripts/aluminion_reporter.py /path/to/run/
python3 scripts/lab_db_updater.py --input_path /path/to/run/
```

`parser.py` performs a preflight check at startup and lists any missing input
files with the corresponding `--skip-*` flag to bypass them.

### Assembly failure handling

If Flye fails to assemble a sample, Aluminion pauses and offers three choices:

| Choice | Effect                                                               |
|:------:|----------------------------------------------------------------------|
|  `1`   | Skip the sample. The rest of the run continues normally.             |
|  `2`   | Retry with `--meta` (tolerates uneven coverage, high-copy plasmids). |
|  `3`   | Stop the pipeline for manual inspection.                             |

Samples skipped here are removed from the internal `samples` tracking file,
so every downstream module (polishing, Bakta, Kleborate, …) ignores them
automatically.

Optional refinement steps (polishing, deconcat, dnaapler, Bandage, the typing
tools) are **non-fatal**: a failure only emits a warning. The assembly is kept
as-is and the run continues. A consolidated warning summary is printed at the
end listing every sample / tool that failed.

---

## Outputs

| File                     | Description                                                         |
|--------------------------|---------------------------------------------------------------------|
| `Aluminion_Report.html`  | **Interactive HTML report — open in any browser, no server needed** |
| `taxonomy.csv` / `.xlsx` | Kraken2 + GAMBIT + Kleborate + ECTyper + MLST per sample            |
| `AbR_modif.xlsx`         | Abricate AMR genes per sample                                       |
| `VF_modif.xlsx`          | Abricate-VFDB virulence genes per sample (`Virulence_genes` column) |
| `mlst_modif.csv`         | MLST scheme, ST, and allele calls                                   |
| `kraken_mlst.xlsx`       | Merged Kraken2 + MLST quick-reference                               |
| `copla_modif.csv`        | Copla plasmid typing (PTU, MOB, Rep, AMR per plasmid)               |
| `integron_summary.csv`   | Integron_Finder results with cassette gene annotations              |
| `phage_summary.csv`      | Phastest prophage regions with completeness scores                  |
| `kleborate.tsv`          | Full Kleborate output (Enterobacterales loci)                       |
| `data_seq.tsv`           | Updated cumulative sequencing database                              |
| `data_analysis.tsv`      | Updated cumulative analysis database                                |
| `alerts.tsv`             | Cross-run MGE alerts raised this run (recurrences + priority hits)  |
| `Alerts_Report.html`     | Interactive report of the MGE alerts                                |

`data_seq_new.tsv` and `data_analysis_new.tsv` are written instead when the
historical databases already exist, so changes can be reviewed before
overwriting the main files.

---

## Cross-run MGE alerts

After consolidating each run, Aluminion ingests its plasmids and integrons into a
persistent **MGE repository** (default `$BASE_DIR/repository`, override with
`--repo`) and compares them against everything seen in previous runs. The goal is
to surface epidemiologically relevant recurrences — the same resistance plasmid or
integron cassette array reappearing across isolates or over time.

- **Plasmid matching** is ANI-based (via `skani`, ≥99% identity by default) with a
  tuple-level fallback on PTU / Rep / MOB when an assembly is unavailable.
- **Integron matching** uses the Jaccard index over cassette gene sets (≥0.8).
- A curated priority-gene catalog (`scripts/_priority_genes.py`) flags
  carbapenemases, *mcr*, hypervirulence loci, etc. With `--alert-new-priority`,
  even a first-occurrence MGE carrying such a gene raises an alert.

Outputs `alerts.tsv` and an interactive `Alerts_Report.html`. The step is
non-fatal and can be skipped with `--no-alerts`. Initialize the repository
explicitly with `--init-repo` (idempotent; also auto-created on first run).

> This repository-backed system replaces the earlier in-database
> `find_shared_mges` comparison that lived in `lab_db_updater.py`; that exact-tuple
> engine has been retired in favour of the ANI/Jaccard matching here.

---

## Repository structure

```
aluminion/
├── aluminion.sh           # Main pipeline orchestrator
├── aluminion_batch.sh     # Wrapper for sequential multi-run processing
├── install.sh             # Automated installer
├── scripts/               # Python parsers and reporters
├── envs/                  # Six conda environment YAMLs
├── examples/              # Example input / output files
└── tests/                 # Automated tests
```

---

## Tests

```bash
python -m pytest tests/ -v
```

`tests/test_parser.py` uses the example files in `examples/` and covers clean
exit, row counts, duplicate detection, AMR gene presence (OXA-48, VIM-1, KPC-2),
and HTML report content. `tests/test_mge_repository.py` and
`tests/test_mge_alerts.py` cover the repository ingestion and alert-matching logic
(they tolerate `skani` being absent — ANI matching is simply skipped in dev). No
bioinformatics tools required to run the suite.

---

## Troubleshooting

**`dorado not found in PATH`** — Install Dorado manually from
<https://github.com/nanoporetech/dorado/releases>.

**`Docker not found`** — `sudo systemctl start docker` and confirm your user is
in the `docker` group.

**Kraken2 runs out of memory** — Your system has less RAM than the database.
Comment out the `/dev/shm` copy in `aluminion.sh` and point `--db` directly to
disk.

**Phastest produces no output** — Verify `$ALUMINION_PHASTEST_DIR` contains
`docker-compose.yml`, `phastest_inputs/`, and `phastest-app-docker/`.

**`parser.py` reports missing files** — The preflight check lists each missing
file and the matching `--skip-*` flag to bypass it.

**Flye `ERROR: No disjointigs were assembled`** — Sample has very high-copy
elements (large plasmids, expression vectors). Choose `2` in the interactive
menu to retry with `--meta`.

**`Failed to initialize NVML: Driver/library version mismatch`** — The running
NVIDIA kernel module no longer matches the user-space library after an
unattended driver upgrade. Reboot, or reload the modules:

```bash
sudo rmmod nvidia_uvm nvidia_drm nvidia_modeset nvidia
sudo modprobe nvidia
nvidia-smi
```

**`CUDA out of memory` during `dorado polish`** — Lower the inference batch
size with `--polish-batchsize 8` (or `4`). Polishing is non-fatal: the
unpolished assembly is kept and the sample is listed in the final warning
summary.

**`list_seq.tsv` not found** — You probably passed a relative path while
running from inside a subfolder. Either `cd` to the parent directory or pass
an absolute path with `-l`.

---

## Contributing

Issues and pull requests welcome. The pipeline is optimised for
Enterobacteriaceae (*Klebsiella*, *Escherichia*, *Enterobacter*,
*Citrobacter*, …) but the assembly, AMR, and annotation modules work for any
bacterial species with a published MLST scheme.
