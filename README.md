# Aluminion
---
Aluminion is an automated, modular pipeline for bacterial whole-genome sequencing (WGS) analysis using Oxford Nanopore Technology (MinION/Mk1D). It takes raw `fastq_pass` reads from the minknow run folder and produces polished assemblies, taxonomic classification, antimicrobial resistance (AMR) profiling, and detection of mobile genetic elements (plasmids, integrons, prophages), all consolidated in an interactive HTML report and two cumulative lab databases.
The key advance in this pipeline is the possibility to keep an updated database of assemblies across runs in a microbiology laboratory or during a project, allowing the creation of comprehensive reports.  
You may find examples for all input tables in the examples folder in this repo.
---

## Pipeline overview

### `aluminion.sh` — main orchestrator

```
 INPUT
 ─────────────────────────────────────────────────────────────────────────────
 MinKNOW run folder          list_seq.tsv         data_seq.tsv (optional)
 (fastq_pass/ + summaries)   (sample–barcode map)  data_analysis.tsv (optional)




 STAGE 1 · READ QC & FILTERING          [conda: aluminion_reads]
 ─────────────────────────────────────────────────────────────────────────────
  NanoPlot ──────────────── QC stats pre-filter ──────────────► QC_reads.csv (pre)
  Chopper ───────────────── quality & length filtering
  NanoPlot ──────────────── QC stats post-filter ─────────────► QC_reads.csv (post)


 STAGE 2 · ASSEMBLY & POLISHING         [conda: aluminion_assembly + aluminion_circlator]
 ─────────────────────────────────────────────────────────────────────────────
  Kraken2 ────────────────── read-level classification ────────► 04_taxonomies/kraken2/
                                                                  genus.csv · species.csv
  Flye ──────────────────── de novo assembly ─────────────────► 03_assemblies/<sample>/
  Dorado polish ─────────── consensus polishing (ONT, GPU opt.)
  Circlator ─────────────── chromosome recircularization       [aluminion_circlator]
  Quast ──────────────────── assembly metrics ────────────────► QC_assembly.csv
  Bandage ────────────────── assembly graph (visual QC)


 STAGE 3 · ANNOTATION, TAXONOMY & AMR  [conda: aluminion_annot + aluminion_kleborate]
 ─────────────────────────────────────────────────────────────────────────────
  GAMBIT ─────────────────── genome-level species ID ──────────► 04_taxonomies/gambit.csv
  Bakta ──────────────────── genome annotation ────────────────► 08_Anotacion/<sample>/
  Abricate ───────────────── AMR gene detection ───────────────► AbR_report.csv
  MLST ───────────────────── Multi-Locus sequence typing ──────► mlst.csv
  Kleborate ──────────────── Enterobacterales loci ────────────► kleborate.tsv (root copy)
  ECTyper ────────────────── E. coli serotyping ───────────────► 04_taxonomies/ectyper/output.tsv


 STAGE 4 · MOBILE GENETIC ELEMENTS
 ─────────────────────────────────────────────────────────────────────────────
  MOB-suite  (Docker) ─────── plasmid reconstruction ──────────► 08_Anotacion/<sample>/mob_recon/
  Copla      (Docker) ─────── plasmid typing ──────────────────► copla.txt

  Phastest   (Docker Compose) prophage detection ──────────────► 09_phages/<sample>/
    └─ phage_parser.py ──────────────────────────────────────────► phage_summary.csv

  Integron_Finder [conda: aluminion_integron]
    └─ integron_parser.py ──── integron annotation ─────────────► integron_summary.csv

  ISfinder BLASTn ─────────── IS element mapping ─────────────► 08_Anotacion/<sample>/IS_chr.tsv
    └─ IS_parser.py ─────────────────────────────────────────────► IS_chr_out.tsv


 STAGE 5 · CONSOLIDATION & REPORTING   [conda: aluminion_annot]
 ─────────────────────────────────────────────────────────────────────────────
  parser.py ──────────────── merge all tool outputs ───────────► taxonomy.csv / .xlsx
              (preflight check                                    AbR_modif.xlsx
               lists any missing                                  mlst_modif.csv
               input files with                                   kraken_mlst.xlsx
               --skip-* hints)
                └─ copla_parser.py ─── copla.txt → ─────────────► copla_modif.csv

  aluminion_reporter.py ───── interactive HTML report ─────────► Aluminion_Report.html

  Datos_seq_unified2.py ───── update lab databases ────────────► data_seq.tsv
              (--init or auto                                     data_analysis.tsv
               on first run)



 OUTPUT SUMMARY
 ─────────────────────────────────────────────────────────────────────────────
 Aluminion_Report.html   Interactive browser report (no server needed)
 taxonomy.xlsx           Taxonomy · AMR · MLST per sample
 data_seq.tsv            Cumulative sequencing QC database (all runs)
 data_analysis.tsv       Cumulative analysis database (all runs)
```



### Standalone scripts (run independently or after `aluminion.sh`)

```
 scripts/parser.py
   Usage : python3 scripts/parser.py -i <run_dir> [--skip-*]
   Input : AbR_report.csv · mlst.csv · kraken2/{genus,species}.csv
           gambit.csv · ectyper/output.tsv · kleborate output
           copla_modif.csv · phage_summary.csv · integron_summary.csv
   Output: taxonomy.csv/.xlsx · AbR_modif.xlsx · mlst_modif.csv · kraken_mlst.xlsx

 scripts/aluminion_reporter.py
   Usage : python3 scripts/aluminion_reporter.py <run_dir>
   Input : list_seq.tsv · QC_reads.csv · QC_assembly.csv
           taxonomy.xlsx · AbR_modif.xlsx · copla_modif.csv
           integron_summary.csv · phage_summary.csv · kleborate.tsv
   Output: Aluminion_Report.html

 scripts/Datos_seq_unified2.py
   Usage : python3 scripts/Datos_seq_unified2.py --input_path <run_dir> [--init]
   Input : list_seq.tsv · QC_reads.csv · QC_assembly.csv · taxonomy.csv
           final_summary_*.txt · report_*.json  (MinKNOW — skipped gracefully if absent)
           copla_modif.csv · phage_summary.csv · integron_summary.csv
           data_seq.tsv · data_analysis.tsv  (created from scratch with --init or if absent)
   Output: data_seq.tsv  (or data_seq_new.tsv in update mode)
           data_analysis.tsv  (or data_analysis_new.tsv in update mode)
```

---

## System requirements

| Resource | Minimum | Recommended |
|---|---|---|
| OS | Linux (Ubuntu 20.04+) | Ubuntu 22.04 LTS |
| CPU | 16 cores | 32+ cores |
| RAM | 64 GB | 128 GB (see Kraken2 note) |
| Disk | 500 GB free | 1 TB+ free |
| GPU | — | NVIDIA GPU ≥16 GB VRAM (for Dorado polish) |

> **macOS** is not supported due to Docker networking requirements for Phastest.  
> **Windows** requires WSL2 (community-tested, not officially supported).

---

## Installation

### 0. Automated installer (recommended)

```bash
git clone https://github.com/Aluminio-visto/aluminion.git
cd aluminion
chmod +x aluminion.sh install.sh
./install.sh -b /$your_database_folder
```

`install.sh` creates all conda environments, pulls Docker images, optionally downloads databases, and installs the `aluminion` command in `~/.local/bin` so the pipeline can be called from anywhere. Flags:

| Flag | Effect |
|---|---|
| `--skip-envs` | Skip conda environment creation |
| `--skip-docker` | Skip Docker image pulls |
| `--skip-dbs` | Skip database downloads |

### 1. Manual installation

#### Clone the repository

```bash
git clone https://github.com/Aluminio-visto/aluminion.git
cd aluminion
chmod +x aluminion.sh
```

#### Install Conda / Mamba

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```

#### Create the conda environments

Aluminion uses six isolated conda environments. Copla runs entirely via Docker and needs no conda env.

```bash
mamba env create -f envs/aluminion_reads.yml       # NanoPlot, Chopper, pillow, kaleido (pip)
mamba env create -f envs/aluminion_assembly.yml    # Kraken2, Flye, QUAST, Bandage, samtools
mamba env create -f envs/aluminion_circlator.yml   # Circlator (isolated — complex legacy deps)
mamba env create -f envs/aluminion_annot.yml       # Bakta, GAMBIT, Abricate, MLST, BLAST, datamash, Python
mamba env create -f envs/aluminion_integron.yml    # Integron_Finder
mamba env create -f envs/aluminion_kleborate.yml   # Kleborate + ECTyper
```

> **Note on `aluminion_reads`**: this environment includes `kaleido` (installed via pip, not conda-forge — it is not available there) and `pillow`. Both are required for NanoPlot to generate static PNG outputs on headless servers without spawning a browser. See the Troubleshooting section if the environment creation hangs or fails.

All six are required for a full run. Individual environments can be omitted if you skip the corresponding module with a `--skip-*` flag.


#### Install Dorado

Dorado is a proprietary basecaller and polisher from Oxford Nanopore Technologies. It is **not available via conda** and must be installed manually.

```bash
# Download the latest release for Linux x86_64:
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-latest-linux-x64.tar.gz
tar -xzf dorado-latest-linux-x64.tar.gz
sudo mv dorado-*/bin/dorado /usr/local/bin/
dorado --version
```

> See https://github.com/nanoporetech/dorado/releases for the latest URL.
> GPU support requires compatible NVIDIA drivers and CUDA ≥12.0; `dorado polish`
> auto-detects the device (no `--device` flag is passed by Aluminion).
> If the GPU runs out of memory, lower the inference batch size with
> `--polish-batchsize <N>` (e.g. `8` or `4`).


#### Install Docker and Docker Compose

```bash
sudo apt-get update && sudo apt-get install -y docker.io docker-compose-plugin
sudo usermod -aG docker $USER && newgrp docker
docker --version && docker compose version
```


#### Pull Docker images

```bash
docker pull kbessonov/mob_suite:3.0.3   # plasmid reconstruction
docker pull rpalcab/copla:1.0           # plasmid clustering
```


#### Set up Phastest (local Docker Compose — no public image)

Phastest does not have a public Docker image. It requires a local Docker Compose setup. Follow the official instructions at https://phastest.ca. By default, Aluminion expects the Phastest folder at `~/Programs/phastest-docker/`. Override with:

```bash
export ALUMINION_PHASTEST_DIR=/path/to/phastest-docker
# or pass -p /path/to/phastest-docker to aluminion.sh
```

---

## Database setup

All databases should reside under a single root directory (e.g., `/$your_database_folder`). Pass this path with the `-b` flag.


### Kraken2 — standard database (~100 GB)

```bash
mkdir -p /$your_database_folder/Kraken && cd /$your_database_folder/Kraken
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz
tar -xzf k2_standard_20240904.tar.gz && rm k2_standard_20240904.tar.gz
```

> **RAM optimization**: With ≥128 GB RAM, Aluminion copies Kraken2 `.k2d` files to `/dev/shm` (RAM disk) before classification and removes them after, giving ~10× speedup. On systems with less RAM, set `--db` to point directly to the disk path and skip the `/dev/shm` copy (edit lines ~172–174 in `aluminion.sh`).

### GAMBIT — genomic species identification

```bash
mkdir -p /$your_database_folder/gambit && cd /$your_database_folder/gambit
wget https://storage.googleapis.com/gambit-public/gambit-db/gambit-signatures-2024-09-01.h5
wget https://storage.googleapis.com/gambit-public/gambit-db/gambit-metadata-2024-09-01.gdb
```

> Check https://github.com/jlumpe/gambit for the latest database release.

### Bakta — bacterial genome annotation (~30 GB full, ~2 GB light)

```bash
mkdir -p /$your_database_folder/bakta
bakta_db download --output /$your_database_folder/bakta --type full
```

Replace `full` with `light` for a faster, less complete setup.

### ISfinder — insertion sequence database

```bash
mkdir -p /$your_database_folder/ISfinder
wget https://www.is-finder.org/download/IS.fna -O /$your_database_folder/ISfinder/ISfinder-nucl.fasta
```

### Abricate — AMR databases (automatic on first use)

```bash
conda activate aluminion_annot
for db in ncbi resfinder card argannot vfdb; do abricate-get_db --db $db; done
```

### MLST — PubMLST schemas (automatic on first use)

The `mlst` tool downloads schemas automatically from PubMLST the first time it runs. No manual setup required.

---

## Input files

### Required for every run

#### `list_seq.tsv` — per-run sample list

Maps each sample to its barcode for the current sequencing run. Keep this file in your **parent working directory** (the `-d` path) and pass it with `-l list_seq.tsv` when you run from that folder, or with `-l /absolute/path/to/list_seq.tsv` from anywhere. Aluminion copies it into the run subfolder at startup, so editing the parent-folder copy between runs is safe.

| Column | Description |
|---|---|
| `Lab_id` | Internal lab culture ID |
| `Strain` | Strain collection code |
| `ID` | **Unique sample identifier** — used as the sample name throughout the entire pipeline |
| `Barcode` | Barcode number assigned by MinKNOW (e.g., `01`, `24`) |
| `DNA_conc` | DNA concentration (ng/µL) — informational |
| `is_repeated` | Mark `x` if this is a repeat sequencing of a previously failed sample |

See `examples/list_seq.tsv` for a filled example.

> The column `ID` is the key identifier. It must match barcode assignments in MinKNOW and will appear as `Sample` in all output files.

### Automatically provided by MinKNOW (copied at run start)

| File | Description |
|---|---|
| `fastq_pass/` | Demultiplexed FASTQ reads, one subfolder per barcode |
| `final_summary_*.txt` | MinKNOW run summary (instrument, flow cell, dates, duration) |
| `report_*.json` | MinKNOW JSON report (pore counts, yield) |

Aluminion copies these automatically from `$MINKNOW_DIR/$RUN_NAME/`. Override the MinKNOW data path with `-m` or `$ALUMINION_MINKNOW_DIR`.

### Historical lab databases (auto-created on first run)

| File | Description |
|---|---|
| `data_seq.tsv` | Cumulative sequencing run database (QC metrics, flow cell info, depth per sample) |
| `data_analysis.tsv` | Cumulative analysis database (taxonomy, AMR, MGE counts per sample) |

These files are updated automatically at the end of each run by `Datos_seq_unified2.py`. On the **first run ever** (when the files do not yet exist), they are created from scratch — no manual setup required. You can also force creation with `--init-db`.

Empty column-schema templates are available in `examples/` for reference.

---

## Usage

### Directory layout

Aluminion organises output around a **parent working directory** (set with `-d`) and a **per-run subfolder** (named after `-r`). The parent is where you keep `list_seq.tsv` and the cumulative databases; every run creates its own child folder beneath it:

```
/home/user/Seqs/Servicio/          ← parent working directory (-d)
├── list_seq.tsv                   ← fill this before each run
├── data_seq.tsv                   ← cumulative sequencing database (all runs)
├── data_analysis.tsv              ← cumulative analysis database (all runs)
├── BAC_2025_NOV_25/               ← run folder, created automatically by Aluminion
│   ├── fastq_pass/
│   ├── 03_assemblies/
│   ├── 08_Anotacion/
│   ├── Aluminion_Report.html
│   └── aluminion_YYYYMMDD_HHMMSS.log
└── BAC_2026_FEB_10/
    └── …
```

**How to run**: navigate to (or stay in) the parent directory, prepare `list_seq.tsv` there, then call `aluminion` passing the file with `-l`. You can also pass an absolute path from anywhere.

```bash
cd /home/user/Seqs/Servicio
aluminion -r BAC_2025_NOV_25 -b /$your_database_folder -t 30 -l list_seq.tsv
```

A log file (`aluminion_YYYYMMDD_HHMMSS.log`) is written inside each run folder, capturing all pipeline output for debugging.

### Standard run

```bash
aluminion -r BAC_2025_NOV_25 -b /$your_database_folder -t 30 -l /path/to/list_seq.tsv
```

### First run (no historical databases yet)

```bash
aluminion -r BAC_2025_NOV_25 -b /$your_database_folder -t 30 -l /path/to/list_seq.tsv --init-db
```

`--init-db` is optional if `data_seq.tsv` and `data_analysis.tsv` do not exist — Aluminion detects this automatically. Use `--init-db` to make the intent explicit or to force a rebuild.

### All flags

| Flag | Description | Default |
|---|---|---|
| `-r / --run` | MinKNOW run folder name **(mandatory)** | — |
| `-b / --db-dir` | Path to databases root directory | `/home/usuario/Databases` |
| `-t / --threads` | Number of CPU threads | `30` |
| `-l / --list` | Path to `list_seq.tsv` | — |
| `-d / --dir` | Base working directory | `/home/usuario/Seqs/Servicio` |
| `-p / --phastest-dir` | Path to local Phastest docker-compose folder | `~/Programs/phastest-docker` |
| `-m / --minknow-dir` | Path to MinKNOW data root | `/var/lib/minknow/data` |
| `--init-db` | Create `data_seq.tsv` / `data_analysis.tsv` from scratch | — |
| `--resume` | Resume interrupted run: skip any step whose output already exists | — |
| `--skip-preprocessing` | Skip NanoPlot + Chopper (reuse existing 01_reads/, 02_filter/) | — |
| `--skip-kraken` | Skip Kraken2 read-level classification | — |
| `--skip-abr` | Skip Abricate AMR resistance gene screen | — |
| `--skip-typing` | Skip all typing tools: GAMBIT, MLST, Kleborate, ECTyper | — |
| `--skip-integrons` | Skip Integron_Finder and integron parsing | — |
| `--skip-plasmids` | Skip Copla plasmid typing (MOB-suite always runs) | — |
| `--skip-phages` | Skip Phastest prophage detection | — |
| `--polish-batchsize <N>` | Override the `dorado polish --batchsize` value. Lower it (e.g. `8`, `4`) when the GPU runs out of memory. Omit to use dorado's default. | — |
| `--just-preprocessing` | Run only Stage 1 (QC + filtering) and exit. Output: `02_filter/<sample>.fastq.gz` | — |
| `--just-assembly` | Run Stages 1–2 (QC + filtering + assembly + polishing + QUAST) and exit. Output: `03_assemblies/<sample>/assembly.fasta` | — |

### Resuming a partial run

#### `--resume` — automatic resume (recommended)

The simplest way to continue after an interruption. Aluminion checks whether each step's expected output already exists before running it, and silently skips it if so. This works even if the pipeline failed in the middle of a sample loop — samples already processed are skipped, and processing continues from the first incomplete sample.

```bash
aluminion -r BAC_2025_NOV_25 -b /$your_database_folder -t 30 -l /path/to/list_seq.tsv --resume
```

`--resume` can be combined with any `--skip-*` flag. For example, to resume and also skip Kraken2 (already ran on a previous attempt):

```bash
aluminion -r BAC_2025_NOV_25 -b /$your_database_folder -t 30 -l /path/to/list_seq.tsv \
  --resume --skip-kraken
```

**Sentinel files used by `--resume`:**

| Step | Skipped if this exists |
|---|---|
| Read concatenation | `01_reads/<sample>.fastq.gz` |
| Pre-filter NanoPlot | `01_reads/QC/<sample>/NanoStats.txt` |
| Chopper | `02_filter/<sample>.fastq.gz` |
| Post-filter NanoPlot | `02_filter/QC/<sample>/NanoStats.txt` |
| Kraken2 | `04_taxonomies/kraken2/<sample>.report` |
| Flye | `03_assemblies/<sample>/assembly.fasta` |
| Polishing | `03_assemblies/<sample>/.polished` |
| Deconcat | `03_assemblies/<sample>/deconcat/assembly_corr.fasta` |
| Circlator | `03_assemblies/<sample>/.circlator_done` |
| QUAST | `03_assemblies/quast/transposed_report.tsv` |
| Bakta | `08_Anotacion/<sample>/<sample>.gbff` |
| MOB-suite | `08_Anotacion/<sample>/mob_recon/` |
| Abricate | `08_Anotacion/<sample>/abricate/<sample>.tab` |
| Integron_Finder | `11_integrons/<sample>/` |
| Copla | `08_Anotacion/<sample>/copla/` |
| GAMBIT | `04_taxonomies/gambit.csv` |
| MLST | `mlst.csv` |
| Kleborate | `04_taxonomies/kleborate/enterobacterales__species_output.txt` |
| ECTyper | `04_taxonomies/ectyper/output.tsv` |
| Phastest | `09_phages/phastest_deep/<sample>/` |
| IS search (BLASTn) | `08_Anotacion/<sample>/IS_chr_out.tsv` |

> IS.tsv is always rebuilt from existing output files, so IS results for completed samples are always included even when the BLASTn step is skipped.

#### `--skip-*` — manual module skip

For coarser control, use `--skip-*` flags to bypass entire modules regardless of whether their output exists. All flags can be combined freely — execution continues with the next module after any skipped one.

```bash
# Skip only Kraken2 (already done) and re-run everything else
aluminion -r BAC_2025_NOV_25 -b /$your_database_folder -t 30 -l /path/to/list_seq.tsv \
  --skip-preprocessing --skip-kraken

# Only re-run the final parsing and report (all assemblies and annotations done)
aluminion -r BAC_2025_NOV_25 -b /$your_database_folder -t 30 -l /path/to/list_seq.tsv \
  --skip-preprocessing --skip-kraken --skip-abr --skip-typing \
  --skip-integrons --skip-plasmids --skip-phages
```

> **`--skip-preprocessing`** requires that `01_reads/`, `02_filter/`, and the `samples` file already exist in the run folder from a previous run. If the `samples` file is missing, the pipeline will stop with an error.

> All other `--skip-*` flags bypass both the tool execution **and** the corresponding parsing step in `parser.py`. The downstream consolidation scripts (`parser.py`, `aluminion_reporter.py`, `Datos_seq_unified2.py`) always run regardless of which modules were skipped.

### Running only the parsing stage (no reads)

If assemblies and annotations are already complete (e.g., resuming a failed run), you can run the Python parsers directly:

```bash
conda activate aluminion_annot

# Generate all result tables
python3 scripts/parser.py -i /path/to/run/

# Generate interactive HTML report
python3 scripts/aluminion_reporter.py /path/to/run/

# Update historical databases
python3 scripts/Datos_seq_unified2.py --input_path /path/to/run/
# On first run: add --init
python3 scripts/Datos_seq_unified2.py --input_path /path/to/run/ --init
```

`parser.py` performs a **preflight check** at startup and lists any missing input files with the corresponding `--skip-*` flag suggestion.

---

## Assembly failure handling

Flye assembly failures are handled interactively. If Flye cannot assemble a sample (e.g., no disjointigs produced), the pipeline pauses and presents three options:

```
  ┌──────────────────────────────────────────────────────────────────────┐
  │  ASSEMBLY FAILED: <sample>
  │  Flye could not assemble this sample (no disjointigs produced).
  │                                                                      │
  │  Options:                                                            │
  │    1) Skip sample and continue with the rest                        │
  │    2) Retry with --meta (for high-copy or fragmented assemblies)    │
  │    3) Stop pipeline for manual inspection                           │
  └──────────────────────────────────────────────────────────────────────┘
  Choose [1/2/3]:
```

| Choice | Behaviour |
|---|---|
| `1` | Sample is excluded from all downstream steps (removed from the internal `samples` list). The rest of the run continues normally. |
| `2` | Flye is retried with `--meta` mode, which tolerates uneven coverage and high-copy elements (plasmids, expression vectors). If `--meta` also fails, the sample is automatically skipped. |
| `3` | Pipeline stops immediately. Re-run with `--skip-preprocessing` once the issue is resolved to avoid repeating the read QC step. |

> **When to use `--meta`**: standard Flye (`--nano-hq`) requires uniform coverage and a minimum overlap depth. Samples containing high-copy-number plasmids or expression vectors may fail because the extreme coverage imbalance confuses the overlap graph. `--meta` is designed for uneven coverage and often succeeds in these cases, though the assembly may be more fragmented.

> Samples skipped due to assembly failure are removed from the `samples` tracking file. All subsequent steps (polishing, Bakta, Kleborate, etc.) iterate over this file, so skipped samples are automatically excluded from every downstream module without any manual intervention.

---

## Polishing internals (`dorado polish`)

`dorado polish` enforces two requirements on its input BAM:

1. The `@PG` header line must contain `PN:dorado`. This means the alignment **must**
   be produced by `dorado aligner` — `minimap2` fails this check with
   *"Input BAM file was not aligned using Dorado"*.
2. The `@RG` header line must contain `DS:basecall_model=<model>`. `dorado aligner`
   does **not** inject this from a FASTQ file, so Aluminion adds it manually with
   `samtools addreplacerg`.

To stay compatible with both the legacy and current Dorado FASTQ formats, Aluminion
inspects the first read header and routes through one of two streaming pipes:

| FASTQ has `RG:Z:` tag? | Strategy |
|---|---|
| Yes (newer Dorado) | `dorado aligner --add-fastq-rg \| samtools sort \| samtools addreplacerg -w` — `-w` overwrites the existing `@RG` so the `DS:basecall_model=…` field replaces whatever was already there. |
| No (older Dorado / re-basecalled FASTQ) | `dorado aligner \| samtools sort \| samtools addreplacerg` — adds a fresh `@RG` line. |

The basecaller model name is parsed from `basecall_model_version_id=…` in the FASTQ
header, or — if that field is missing — from the legacy `dna_<model>` token in the
read description. If neither can be detected, polishing is skipped for that sample
and the unpolished assembly is kept (the sample is listed under
`Unpolished (dorado polish failed)` at the end of the run).

The streaming pipe avoids writing intermediate sorted BAMs to disk; only the final
read-group-tagged BAM is materialised and indexed for `dorado polish`.

---

## Intermediate files (generated during the run)

These files are created by `aluminion.sh` before `parser.py` runs. They are kept in the working directory and can be used to resume or debug:

| File | Generated by | Used by |
|---|---|---|
| `QC_reads.csv` | NanoPlot + shell | reporter, Datos_seq |
| `QC_assembly.csv` | Quast + `cut` | reporter, Datos_seq |
| `kraken.csv` | Kraken2 + awk | (reference only) |
| `04_taxonomies/kraken2/genus.csv` | awk from Kraken2 reports | parser.py |
| `04_taxonomies/kraken2/species.csv` | awk from Kraken2 reports | parser.py |
| `04_taxonomies/gambit.csv` | GAMBIT | parser.py |
| `04_taxonomies/ectyper/output.tsv` | ECTyper | parser.py |
| `04_taxonomies/kleborate/enterobacterales__species_output.txt` | Kleborate | parser.py |
| `kleborate.tsv` | copied from above | aluminion_reporter.py |
| `AbR_report.csv` | copied from `08_Anotacion/AbR.tab` | parser.py |
| `mlst.csv` | MLST tool | parser.py |
| `copla.txt` | Copla Docker | copla_parser.py |
| `copla_modif.csv` | copla_parser.py | reporter, Datos_seq |
| `phage_summary.csv` | phage_parser.py | reporter, Datos_seq |
| `integron_summary.csv` | integron_parser.py | reporter, Datos_seq |

---

## Output files

| File | Description |
|---|---|
| `taxonomy.csv / .xlsx` | Consolidated taxonomy: Kraken2 + GAMBIT + Kleborate + ECTyper + MLST |
| `AbR_modif.xlsx` | Abricate AMR gene summary, one row per sample |
| `mlst_modif.csv` | MLST scheme, ST, and allele calls |
| `kraken_mlst.xlsx` | Merged Kraken2 + MLST quick-reference table |
| `copla_modif.csv` | Copla plasmid typing (PTU, MOB, Rep, AMR genes per plasmid) |
| `integron_summary.csv` | Integron_Finder results with cassette gene annotations |
| `phage_summary.csv` | Phastest prophage regions with completeness scores |
| `kleborate.tsv` | Full Kleborate output (virulence/resistance loci for Enterobacterales) |
| `Aluminion_Report.html` | **Interactive HTML report** — open in any browser, no server needed |
| `data_seq.tsv` | Updated cumulative sequencing database (all runs, all samples) |
| `data_analysis.tsv` | Updated cumulative analysis database (taxonomy, AMR, MGE counts) |

> `data_seq_new.tsv` and `data_analysis_new.tsv` are written instead when the historical databases already exist (update mode), so you can review the changes before overwriting the main files.

---

## Repository structure

```
aluminion/
├── aluminion.sh                  # Main pipeline orchestrator
├── install.sh                    # Automated installer
├── scripts/
│   ├── parser.py                 # Table consolidation (runs after all tools finish)
│   ├── aluminion_reporter.py     # Interactive HTML report generator
│   ├── Datos_seq_unified2.py     # Historical database updater
│   ├── deconcat.py               # Assembly deconcatenation before Circlator
│   ├── copla_parser.py           # Copla text output → copla_modif.csv (called by parser.py)
│   ├── phage_parser.py           # Phastest output → phage_summary.csv
│   ├── integron_parser.py        # Integron_Finder output → integron_summary.csv
│   └── IS_parser.py              # ISfinder BLAST output → IS_chr_out.tsv
├── envs/
│   ├── aluminion_reads.yml       # NanoPlot, Chopper, pillow, kaleido (via pip)
│   ├── aluminion_assembly.yml    # Kraken2, Flye, QUAST, Bandage, samtools
│   ├── aluminion_circlator.yml   # Circlator (isolated — complex legacy deps)
│   ├── aluminion_annot.yml       # Bakta, GAMBIT, Abricate, MLST, BLAST, datamash, Python
│   ├── aluminion_integron.yml    # Integron_Finder
│   └── aluminion_kleborate.yml   # Kleborate + ECTyper
├── examples/                     # Example input/output files for testing
│   ├── list_seq.tsv
│   ├── data_seq.tsv              # Empty schema template
│   ├── data_analysis.tsv         # Empty schema template
│   └── …
└── tests/
    └── test_parser.py            # Automated tests for parser.py + aluminion_reporter.py
```

---

## Running the tests

No bioinformatics tools required — tests use the example files in `examples/`:

```bash
# With pytest:
python -m pytest tests/test_parser.py -v

# Without pytest:
python tests/test_parser.py
```

The test suite covers: clean exit, row counts, duplicate detection, key column completeness, AMR gene presence (OXA-48, VIM-1, KPC-2), and HTML report content.

---

## Troubleshooting

**`dorado not found in PATH`** — Dorado must be installed manually; see Installation step above.

**`Docker not found`** — Install Docker and ensure the daemon is running (`sudo systemctl start docker`).

**Kraken2 runs out of memory** — Your system has less RAM than the Kraken2 database. Comment out the `/dev/shm` copy lines (~172–174) in `aluminion.sh` and point `--db` directly to the disk database. You can also add `--memory-mapping` to Kraken2 to memory-map from disk instead of loading into RAM.

**Phastest produces no output** — Verify that `$PHASTEST_DIR` contains `docker-compose.yml` and the `phastest_inputs/` and `phastest-app-docker/` subdirectories as described in the Phastest setup guide.

**`parser.py` reports missing files** — The preflight check lists each missing file and the `--skip-*` flag to bypass it. Use the flags during partial runs (e.g., when phage/integron analysis was not performed).

**`[ERROR] Empty 'list_seq.tsv' created`** — Aluminion could not find the file at the path you passed with `-l`. If you used a relative path (e.g., `-l list_seq.tsv`), make sure you are in the parent working directory when you call `aluminion`, not inside a run subfolder. Alternatively, pass an absolute path: `-l /home/user/Seqs/Servicio/list_seq.tsv`.

**`ERROR: No disjointigs were assembled` (Flye)** — Flye failed to assemble the sample. This typically happens with samples containing very high-copy-number elements (plasmids, expression vectors) that create extreme coverage imbalance. When this occurs, Aluminion pauses and offers three options: skip the sample, retry with `--meta` (recommended for high-copy cases), or stop the pipeline. See [Assembly failure handling](#assembly-failure-handling) for details.

**NanoPlot hangs or produces Chrome/choreographer warnings on headless servers** — NanoPlot uses Plotly for interactive plots and tries to render static PNGs via `choreographer`, which spawns a browser. On servers without a display this can stall. The `aluminion_reads` environment ships `kaleido` (pip) and `pillow`, and all NanoPlot calls in `aluminion.sh` set `PLOTLY_RENDERER=kaleido MPLBACKEND=Agg` to route rendering through kaleido instead of the browser. If you see this issue, make sure the `aluminion_reads` environment was created with the current `envs/aluminion_reads.yml` (delete and recreate if it predates this fix: `conda env remove -n aluminion_reads && mamba env create -f envs/aluminion_reads.yml`).

**`data_seq.tsv` not found on first run** — This is expected. Either add `--init-db` to your first run or let Aluminion auto-detect and create the databases from scratch.

**taxonomy.csv has duplicate rows** — This is caused by duplicate entries in `kraken2/genus.csv` or `kraken2/species.csv`, usually from re-running the pipeline without clearing those files. The current version clears these files before re-appending (`> genus.csv`, `> species.csv`) to prevent accumulation.

**`Failed to initialize NVML: Driver/library version mismatch`** — The running NVIDIA kernel module no longer matches the user-space library after an unattended driver upgrade. `dorado polish` will fall back to CPU (very slow) or fail outright. Two ways to fix:

```bash
# Quickest: reboot to load the matching kernel modules
sudo reboot

# Or, unload and reload the NVIDIA modules without rebooting (requires no GPU processes):
sudo rmmod nvidia_uvm nvidia_drm nvidia_modeset nvidia
sudo modprobe nvidia
nvidia-smi   # should now report correctly
```

If the error persists, reinstall the driver matching your current kernel
(`uname -r`) and CUDA version (`nvidia-smi`).

**`CUDA out of memory` / `no kernel image available` during `dorado polish`** — The
inference batch size is too large for the GPU's free VRAM. Lower it with
`--polish-batchsize 8` (or `4`). Polishing is automatically skipped (non-fatal) if
it still fails — the unpolished assembly is kept and the sample is listed in the
final warning summary.

---

## Contributing

Issues and pull requests are welcome. The pipeline is optimized for Enterobacteriaceae (*Klebsiella*, *Escherichia*, *Enterobacter*, *Citrobacter*...) but the assembly, AMR, and annotation modules work for any bacterial species with an appropriate MLST scheme.
