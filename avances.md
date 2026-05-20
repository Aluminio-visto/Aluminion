# Aluminion — Development log and progress summary

---

## 1. Project philosophy and objectives

### Primary objective

Aluminion is a fully automated, modular pipeline for **bacterial whole-genome sequencing (WGS)** analysis from Oxford Nanopore Technology (MinION / Mk1D) reads. The pipeline covers the complete analysis cycle from raw reads to an interactive HTML report, and maintains **cumulative lab databases** (`data_seq.tsv`, `data_analysis.tsv`) that grow across sequencing runs, enabling longitudinal tracking of isolates in a clinical microbiology setting.

**Target organisms:** Enterobacteriaceae — primarily *Klebsiella pneumoniae*, *Escherichia coli*, *Enterobacter* spp., *Citrobacter* spp.
**Primary use case:** Serial surveillance in a hospital microbiology laboratory; detecting and tracking AMR determinants and mobile genetic elements (MGEs) across runs.

### Secondary objectives

- Detect and characterize **mobile genetic elements** (plasmids via MOB-suite + Copla; integrons via IntegronFinder; prophages via Phastest; insertion sequences via BLASTn vs ISfinder) to support outbreak tracking.
- Enable **partial runs** — the user should be able to run only what they need (just filtered reads, just an assembly, just annotation) without running the full pipeline.
- Allow **resuming interrupted runs** at any point without repeating completed work.
- Remain **non-fatal** for optional refinement steps: polishing, deconcatenation, circularization, and MGE detection are valuable but not required to produce a valid assembled genome.
- Be **deployable on a headless Linux server** (SSH, no display, no GPU required).

### Future direction

- Migration to **Nextflow DSL2** for portability and HPC job scheduling.
- Cross-run MGE comparison: detect shared plasmids / integron cassettes across sequential runs (`data_analysis.tsv` as the reference).

---

## 2. Coding norms and conventions

All code and comments in **American English** (the pipeline targets an international scientific audience).

### Bash (`aluminion.sh`)

- `set -euo pipefail` always active — exit on error, unset variable, or pipe failure.
- All paths come from arguments or environment variables; **no hardcoded absolute paths**.
- Timestamped `log()` (green), `error_log()` (red), `warn()` (yellow) functions for consistent output.
- Conda environments are activated inline with `conda activate <env>` — the script sources `conda.sh` at startup.
- Background jobs (`&`) are used only where the next step does not depend on the output (e.g., NanoPlot QC alongside Chopper filtering); `wait` is avoided when it would block indefinitely.
- Sentinel files (`.polished`, `.circlator_done`) mark completion of steps that overwrite their own input, enabling correct `--resume` detection.

### Python scripts (`scripts/`)

- One script per logical task; `parser.py` is the main orchestrator and calls sub-parsers internally.
- All column names and output headers in English.
- No hardcoded paths — all paths received as CLI arguments.
- Graceful handling of missing optional input files (log and skip, do not crash).

### Conda environments

- One environment per tool group to isolate conflicting dependency trees.
- Pin Python to `3.12` in environments that use pip-installed binary packages (kaleido); Python ≥3.13 has free-threading changes that break binary wheels.
- Prefer conda-forge over pip; use pip only when a package is not available in any conda channel.

---

## 3. Confirmed working changes

### 3.1 `aluminion.sh` — orchestrator

#### Relative path resolution (pre-`cd` fix)
`SEQ_LIST_INPUT` and `BASE_DIR` are resolved to absolute paths with `readlink -f` **immediately after argument parsing**, before any `cd "$WORKDIR"` call. This prevents `list_seq.tsv` not found errors when the pipeline is launched from the parent directory with a relative path.

```bash
[ -n "$SEQ_LIST_INPUT" ] && SEQ_LIST_INPUT="$(readlink -f "$SEQ_LIST_INPUT")"
[ -n "$BASE_DIR"       ] && BASE_DIR="$(readlink -f "$BASE_DIR")"
```

#### Pipeline log file
Every run writes a timestamped log to `RUN_NAME/aluminion_YYYYMMDD_HHMMSS.log` by duplicating stdout and stderr via `tee`:

```bash
LOG_FILE="${WORKDIR}/aluminion_$(date +'%Y%m%d_%H%M%S').log"
exec > >(tee -a "$LOG_FILE") 2>&1
```

#### Helper functions
```bash
log()       { echo -e "\n\033[1;32m[$(date +'%Y-%m-%d %H:%M:%S')] $1\033[0m"; }
error_log() { echo -e "\n\033[1;31m[ERROR] $1\033[0m"; }
warn()      { echo -e "\033[1;33m[WARNING] $1\033[0m"; }
resume_done() { [ -n "$RESUME" ] && { [ -f "$1" ] || [ -d "$1" ]; }; }
```

`resume_done <path>` returns true when `--resume` is active **and** the sentinel file or directory exists. Used before every tool call.

#### `--resume` flag — per-sample sentinel checks
Every tool in the pipeline is individually guarded:

| Step | Sentinel |
|------|----------|
| Read concatenation | `01_reads/<sample>.fastq.gz` |
| Pre-filter NanoPlot | `01_reads/QC/<sample>/NanoStats.txt` |
| Chopper | `02_filter/<sample>.fastq.gz` |
| Post-filter NanoPlot | `02_filter/QC/<sample>/NanoStats.txt` |
| Kraken2 | `04_taxonomies/kraken2/<sample>.report` |
| Flye | `03_assemblies/<sample>/assembly.fasta` |
| Dorado polish | `03_assemblies/<sample>/.polished` |
| Deconcat | `03_assemblies/<sample>/deconcat/assembly_corr.fasta` |
| Circlator | `03_assemblies/<sample>/.circlator_done` |
| QUAST | `03_assemblies/quast/transposed_report.tsv` |
| Bakta | `08_Anotacion/<sample>/<sample>.gbff` |
| MOB-suite | `08_Anotacion/<sample>/mob_recon/` |
| Abricate | `08_Anotacion/<sample>/abricate/<sample>.tab` |
| IntegronFinder | `11_integrons/<sample>/` |
| Copla | `08_Anotacion/<sample>/copla/` |
| GAMBIT | `04_taxonomies/gambit.csv` |
| MLST | `mlst.csv` |
| Kleborate | `04_taxonomies/kleborate/enterobacterales__species_output.txt` |
| ECTyper | `04_taxonomies/ectyper/output.tsv` |
| Phastest | `09_phages/phastest_deep/<sample>/` |
| IS BLASTn | `08_Anotacion/<sample>/IS_chr_out.tsv` |

Polishing and circlator use `.polished` / `.circlator_done` **touch-files** (not the assembly itself) because both steps overwrite `assembly.fasta`.

#### Kraken2 resume fix — avoid unnecessary 100 GB `/dev/shm` copy
Before copying the Kraken2 database to RAM disk, the pipeline checks whether any sample still needs classification. If all reports exist, the expensive copy is skipped entirely:

```bash
kraken_needed=false
for i in $(cat samples); do
    resume_done "04_taxonomies/kraken2/${i}.report" || { kraken_needed=true; break; }
done
if [ "$kraken_needed" = true ]; then
    cp ${KRAKEN_DB}/*.k2d /dev/shm/
    # ... run kraken2 per sample ...
    rm -f /dev/shm/*.k2d
else
    log "  [resume] All Kraken2 reports found — skipping database copy."
fi
```

#### `--skip-preprocessing` flag
Skips NanoPlot pre-filter, Chopper, and NanoPlot post-filter. Reads the existing `samples` file; exits with an error if it is missing.

#### Skip flags properly wrap tool execution
`--skip-kraken`, `--skip-abr`, `--skip-typing`, `--skip-integrons`, `--skip-plasmids`, `--skip-phages` each wrap the actual tool execution block (not just the downstream parser call). Tool semantics:

| Flag | Tools skipped |
|------|---------------|
| `--skip-kraken` | Kraken2 classification |
| `--skip-abr` | Abricate AMR screen |
| `--skip-typing` | GAMBIT, MLST, Kleborate, ECTyper |
| `--skip-integrons` | IntegronFinder + integron_parser.py |
| `--skip-plasmids` | Copla plasmid typing (MOB-suite always runs) |
| `--skip-phages` | Phastest + phage_parser.py |

`Bakta` is always executed (core annotation, not skippable).

#### Early-stop flags (`--just-*`)
A single `STOP_AFTER` variable controls clean exit after a named stage:

```bash
--just-preprocessing  →  exit after Chopper (output: 02_filter/<sample>.fastq.gz)
--just-assembly       →  exit after QUAST   (output: 03_assemblies/<sample>.fasta)
```

Implemented as two one-liners at the stage boundaries:
```bash
[ "$STOP_AFTER" = "preprocessing" ] && { log "..."; exit 0; }
[ "$STOP_AFTER" = "assembly" ]      && { log "..."; exit 0; }
```

Compatible with `--resume` and `--skip-kraken`.

#### Flye interactive failure handler
When Flye fails, the pipeline pauses and presents a 3-choice menu instead of crashing:

```
1) Skip sample — continue with remaining samples
2) Retry with --meta (high-copy / fragmented assemblies)
3) Stop pipeline for manual inspection
```

Skipped samples are removed from the `samples` tracking file, so all downstream loops (polishing, Bakta, Kleborate, etc.) ignore them automatically without any manual intervention.

#### Polishing — minimap2 + non-fatal dorado polish
The alignment step was switched from `dorado aligner` to `minimap2` (already in `aluminion_assembly` env):

```bash
minimap2 -ax map-ont -t $THREADS_TOTAL 03_assemblies/${i}/assembly.fasta 02_filter/${i}.fastq.gz \
    | samtools sort -@ $THREADS_TOTAL -o 03_assemblies/${i}/${i}_aligned_reads.bam
samtools index -@ $THREADS_TOTAL 03_assemblies/${i}/${i}_aligned_reads.bam
```

**Reason:** `dorado aligner` requires basecaller model metadata in the BAM header (embedded by `dorado basecall`). Reads processed through Chopper become plain FASTQ and lose that metadata, causing `dorado polish` to fail with "Input BAM file has no basecaller models listed in the header." `minimap2` has no such requirement.

`dorado polish` is wrapped in a non-fatal block with `--device cpu` to avoid GPU NVML errors:

```bash
if dorado polish --threads $THREADS_TOTAL --device cpu ... ; then
    mv polished_assembly.fasta assembly.fasta
    touch .polished
else
    warn "Polishing failed for ${i}. Assembly kept unpolished."
    failed_polish+=("$i")
fi
```

#### Non-fatal deconcat and circlator
Same pattern: failures append the sample to `failed_deconcat[]` or `failed_circlator[]` and the pipeline continues using `assembly.fasta` as-is.

#### Final warning summary
At the end of the pipeline, a yellow warning block lists all samples that completed assembly but failed optional refinement:

```bash
[ ${#failed_polish[@]}    -gt 0 ] && warn "  Unpolished  : ${failed_polish[*]}"
[ ${#failed_deconcat[@]}  -gt 0 ] && warn "  Deconcat    : ${failed_deconcat[*]}"
[ ${#failed_circlator[@]} -gt 0 ] && warn "  Circlator   : ${failed_circlator[*]}"
```

#### NanoPlot — headless server Chrome fix
NanoPlot 1.46.2 uses `choreographer` → Chrome for static PNG rendering. On headless Linux servers, Chrome requires `--no-sandbox`. A wrapper script is created at runtime and exported via `BROWSER_PATH`:

```bash
CHROME_REAL=$(ls "$HOME/mambaforge/envs/aluminion_reads/lib/python"*/site-packages/choreographer/cli/browser_exe/chrome-linux64/chrome 2>/dev/null | head -1)
if [ -n "$CHROME_REAL" ] && [ -x "$CHROME_REAL" ]; then
    CHROME_WRAPPER="${WORKDIR}/.chrome_wrapper"
    printf '#!/bin/bash\nexec "%s" --no-sandbox --disable-gpu --disable-dev-shm-usage "$@"\n' "$CHROME_REAL" > "$CHROME_WRAPPER"
    chmod +x "$CHROME_WRAPPER"
    export BROWSER_PATH="$CHROME_WRAPPER"
fi
```

The glob path (`python*/`) is used instead of calling Python, because `conda activate` does not always update `$PATH` for inline Python calls in non-interactive scripts.

`MPLBACKEND=Agg` is set on all NanoPlot calls to prevent matplotlib from opening a display window.

`wait` was removed from the NanoPlot loops. NanoPlot (reads `01_reads/`) and Chopper (writes to `02_filter/`) operate on different directories — there is no file dependency between them, so NanoPlot can run in the background (`&`) while Chopper processes reads in parallel. The old `set +e`/`set -e` + `wait` pattern was required only when `wait` was present; without it, the `&` exit code never reaches the main shell.

---

### 3.2 `envs/aluminion_reads.yml`

```yaml
name: aluminion_reads
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.12   # Pin: Python ≥3.13 has free-threading changes that break binary wheels
  - nanoplot
  - chopper
  - pillow        # DPI metadata embedding in NanoPlot PNG output
  - pip
  - pip:
    - kaleido     # Plotly static image export; not on conda-forge; NanoPlot 1.46.2 requires ≥1.0.0
```

Key decisions:
- `kaleido` must be installed via **pip** (not available on conda-forge or bioconda).
- `kaleido 1.x` internally uses `choreographer` → Chrome. The Chrome wrapper above is still required.
- `pillow` is needed for DPI metadata embedding in NanoPlot PNG output.
- `python=3.12` pin: Python 3.13/3.14 produced `pysam` GIL warnings and may break kaleido binary wheels.

---

### 3.3 `README.md` — documentation

Updated sections:
- **Directory layout** — ASCII tree showing the parent/child folder structure (`list_seq.tsv` lives in the parent directory; each run creates a subfolder).
- **All flags table** — complete with `--resume`, `--skip-preprocessing`, `--just-preprocessing`, `--just-assembly`, and all `--skip-*` flags with correct descriptions.
- **Resuming a partial run** — full sentinel table + `--skip-*` combination examples.
- **Assembly failure handling** — Flye interactive menu with choice table and guidance on when to use `--meta`.
- **`aluminion_reads` environment note** — explains that kaleido is installed via pip, why, and what to do if the env predates the fix.
- **Troubleshooting** — entries for NanoPlot Chrome hang, Flye failure, `list_seq.tsv` not found, `data_seq.tsv` not found on first run, and Kraken2 RAM usage.

---

### 3.4 `CLAUDE.md` — AI coding context

Added a full `### Aluminion` project-specific section containing:
- Repository file tree with descriptions.
- Run-time folder layout (what `aluminion.sh` creates under `RUN_NAME/`) including all sentinel files.
- Pipeline stages table (stage → folders → tools → conda env → skip flag).
- **Key implementation decisions** (non-obvious choices that must not be reverted):
  - Why `set +e`/`set -e` was abandoned in favor of removing `wait`.
  - Why kaleido is pip-only and what version.
  - Why `.polished` / `.circlator_done` exist instead of using `assembly.fasta` as sentinel.
  - Why `readlink -f` must happen before `cd`.
  - The Chrome wrapper rationale.
  - Why `copla.txt` is only truncated on fresh runs (not on `--resume`).

---

## 4. Pending / next steps

| Priority | Task |
|----------|------|
| High | Verify full pipeline run end-to-end: NanoPlot → Chopper → Assembly → Annotation → Report |
| High | Confirm `dorado polish` behavior: if it consistently fails due to missing model header, evaluate replacing with Medaka |
| Medium | Cross-run MGE comparison script: detect shared plasmids / integron cassettes between runs using `data_analysis.tsv` |
| Low | Consider adding `--just-annotation` early-stop flag |
| Low | Nextflow DSL2 migration (when the bash pipeline is fully stable) |

---

## 5. Architecture reference — conda environments

| Environment | Key tools |
|-------------|-----------|
| `aluminion_reads` | NanoPlot, Chopper, pillow, kaleido (pip) |
| `aluminion_assembly` | Kraken2, Flye, QUAST, Bandage, samtools, minimap2, blast, mafft, emboss, pandas, matplotlib |
| `aluminion_circlator` | circlator |
| `aluminion_annot` | Bakta, ABRicate, BLAST, MOB-suite, GAMBIT, mlst, ECTyper, Python stack |
| `aluminion_integron` | IntegronFinder |
| `aluminion_kleborate` | Kleborate |

Docker images: `kbessonov/mob_suite:3.0.3`, `rpalcab/copla:1.0`, Phastest (local docker-compose).
External binaries: `dorado` (in `$PATH`).

---

## 6. Session log — 2026-05-20 (Opus 4.7 review + bug-fix sweep)

### Context
After switching from Sonnet 4.6 to Opus 4.7, a full pipeline review was requested. A run on the Mutantes Klebsiella dataset surfaced two acute blockers (deconcat failing for all samples, QUAST not found) that pre-empted the planned refactor work. A prioritized task list was agreed (P0 → P4); this section logs what landed and what is still open.

### Bugs fixed in this session (all in repo, ready to test)

| Tag | File(s) | Fix |
|---|---|---|
| **B10** | `aluminion.sh` | `tr -d '\r'` + per-field whitespace trim before parsing `list_seq.tsv`, so Excel-exported TSVs (CRLF) no longer break barcode lookups. |
| **B7** | `aluminion.sh` | IS.tsv rebuild made explicit and resume-safe — the per-sample `IS_chr_out.tsv` files keep their own sentinel; IS.tsv is regenerated from them each run. |
| **B8** | `aluminion.sh` | `[ -s "$file" ]` guards around the IS aggregation (`head -n 1 N_IS_*.tsv`, `tail+wc -l IS_chr_out.tsv`). Previously the `\|\| echo` fallback never fired because `head`/`wc` return 0 on empty files; IS.tsv would contain rows with bare `\t\t`. |
| **B3** | `aluminion.sh` | Chrome (used by NanoPlot's choreographer) is now detected across `~/mambaforge`, `~/miniforge3`, `~/miniconda3` instead of mambaforge only. |
| **B2** | `aluminion.sh` | Auto-detects legacy Spanish column headers in `list_seq.tsv` (`Cultivo`, `Cepa`, …) and rewrites them to the current English schema (`Lab_id`, `Strain`, `ID`, `Barcode`, `DNA_conc`, `is_repeated`). Warns the user when it triggers. |
| **B6** | `scripts/integron_parser.py` | `os.path.basename(os.path.dirname(input_path))` replaces `input_path.split('/')[-2]`. Portable across Windows/POSIX path separators. |
| **B5** | `scripts/integron_parser.py` | `prokka_parse` no longer uses the `df.loc[-1] = …; df.index += 1` antipattern. Rows are accumulated in a plain list and converted with `pd.DataFrame(rows, columns=…)` at the end. Removes the risk of index collisions. |
| **B4** | `scripts/phage_parser.py` | `read_summary` locates the header by scanning for a line starting with `REGION` instead of `skiprows=32`. Tolerant to Phastest changing its preamble length between releases. |
| **B1** | `scripts/parser.py` | Removed the dead `pass` block under `--include-run-info`, the unused `import Datos_seq_unified2 as run_info_parser`, and the unused argparse argument. Per user decision: MinKNOW instrument / flowcell metadata is NOT to be merged into `taxonomy.csv`. |
| **subprocess check** | `scripts/integron_parser.py`, `scripts/phage_parser.py` | `subprocess.run` calls for `prokka`, `abricate`, `makeblastdb` now use `check=True` and capture stderr. `subprocess.run(['cp', …])` replaced with `shutil.copy`. blastn `Popen` now waits and raises on non-zero rc. |

### Critical infrastructure fixes (the acute blockers)

| Fix | What changed | Why |
|---|---|---|
| **deconcat env deps** | `envs/aluminion_assembly.yml` — added `blast`, `mafft`, `emboss`, `pandas`, `matplotlib` | `deconcat.py` was failing for every sample with `ModuleNotFoundError: No module named 'pandas'` (and would have hit `blastn`, `mafft`, `em_cons` next). The deconcat.log confirmed pandas as the first missing import. **Action required:** `mamba env update -f envs/aluminion_assembly.yml --prune` before resuming any run. |
| **deconcat preflight** | `aluminion.sh` — `for dep in blastn makeblastdb mafft em_cons minimap2; do command -v "$dep" \|\| error; done` before the deconcat loop | Surfaces the "missing dep" failure mode loudly instead of letting all samples fall into the `warn` branch with no hint of the root cause. Also adds "See 03_assemblies/deconcat.log" to the warn message. |
| **QUAST env / binary** | `aluminion.sh` — explicit `conda activate aluminion_assembly` before the QC block; `quast.py` → `quast` | After the circlator block leaves us in `aluminion_circlator`, QUAST was never reactivated → `quast.py: orden no encontrada`. QUAST 5.x ships the binary as `quast` (no `.py` suffix) on conda-forge / bioconda. |

### Repo hygiene

- **`.gitignore`** created (was missing). Excludes `__pycache__/`, sequencing data (`.fastq*`, `.fast5`, `.pod5`, `.bam`, `.bai`), editor/IDE files, OS junk, and `.claude/worktrees/`.
- `scripts/__pycache__/` purged from the index via `git rm -rf --cached`.
- **Dead env YMLs removed**: `envs/aluminion_base.yml`, `envs/aluminion_copla.yml`, `envs/environment.yml`. Confirmed unused by current `aluminion.sh`, `install.sh`, and `README.md`. The remaining six envs are the canonical set.

### Memory of clarifications made during the session

- **B5 motivation**: integron_parser handles ≤100 elements; the antipattern was an `IndexError`/collision risk, not a speed problem. Refactor was for correctness, not performance.
- **E5 (Kraken awk → pandas) — declined**: the `.report` indent-based parsing is robust in awk and pandas wouldn't be faster on files this small.
- **E6 (drop failed samples from QUAST) — declined**: the user prefers NA rows in `QC_assembly.csv` so they're visible in the HTML report and `data_seq.tsv` — visual signal for "re-sequence this".
- **E7 (Kraken /dev/shm OOM risk) — declined**: 128 GB RAM machine has never OOM'd here; no change.
- **B9 (Copla contig name extracted by character position `${i:(-11):5}`) — deferred**, working in practice.
- **C5 (logging)**: confirmed `logging` is Python stdlib, no env changes needed.
- **C6 (cumulative repositorio + MGE alerts) — deferred to its own session**: needs design discussion (match criteria for alerts, symlink vs copy policy for reads, alert format, lookup strategy for prior `data_seq.tsv` / `data_analysis.tsv`).

### Outstanding from the agreed task list

Priority order to resume in next session, after the user verifies the deconcat + QUAST fixes work end-to-end.

**Tier A — quick wins (small, safe, high leverage):**
- **E2** — Add configurable `--batchsize` to `dorado polish`. If the 4090 chokes at the chosen default, fall back; user explicitly asked to be reminded if this happens.
- **E4** — Stream `dorado aligner | samtools sort | samtools addreplacerg` in a single pipe; eliminates one BAM round-trip per sample.
- **README updates:**
  - Remove the now-stale `--device cpu` note (dorado auto-detects GPU).
  - Document the two-path `@RG` injection logic for polishing (RG:Z: present vs absent).
  - Troubleshooting entry for `NVML driver/library version mismatch` → `sudo reboot` (or `rmmod nvidia_* && modprobe`).
  - Update CUDA version note if outdated.
- **C9** — Document `--init-db` in `--help` and README consistently.

**Tier B — refactor (more code, higher payoff):**
- **C2** — Rename Spanish variables in `parser.py` (`cabecera`, `muestra`, `posibles`, `genes_seguros`, `provis1/2`, `intermedio`, `resultado`, `ultimo_df`) to English. **Critical**: must verify every caller/script that references these — `aluminion_reporter.py`, `Datos_seq_unified2.py` — stays coherent. No half-renaming.
- **C8** — Refactor the opaque names that survive C2 (`intermedio`, `ultimo_df`, `provis1`, `provis2`) into descriptive ones.
- **C5** — Migrate ANSI-coded `print()` calls in all parsers to `logging` (stdlib, no env change).
- **C4** — Extract `safe_read_csv` to `scripts/_utils.py` and import from there.
- **E3** — Refactor MLST processing in `parser.py` from line-by-line file I/O to a single DataFrame build + `.to_csv()` (≤200 rows so speed isn't the win — clarity is).
- **C3** — Move magic constants (`135M` filter, Chopper `-q 12 -l 300 --headcrop 20`, Abricate `--minid 75 --mincov 75`) to top-of-file in `aluminion.sh` and expose as CLI flags (`--min-read-mb`, `--chopper-q`, …). Add to `--help` and README.

**Tier C — rename + docs (low risk, high readability):**
- **C1** — Rename `Datos_seq_unified2.py` → `lab_db_updater.py`. Update import in `parser.py`, call in `aluminion.sh`, references in README and `--help`. Also rename the function `parse_minion_sum` if it survives.
- **README — `repositorio/` folder**: explain its purpose (cross-run cumulative reads/assemblies/MGEs storage). Add to the directory layout diagram and to the abstract at the top of the README. Prepares ground for C6.
- **README — Datos_seq_unified2 / lab_db_updater config**: surface the hard-coded `Depth > 30.0` cutoff and `DNeasy Blood & Tissue` extraction kit. Either make them CLI args or document them explicitly.
- **README — MGE comparison**: cross-reference C6 design (when implemented).

**Tier D — deferred to a dedicated session (architectural):**
- **C6 — Cumulative lab repository + MGE alert system**. Four sub-items, all need design discussion before code:
  - C6.1: explain the cumulative repository idea in the README abstract.
  - C6.2: locate previous `data_seq.tsv` / `data_analysis.tsv` from the last run and seed the current one. Decide: lookup strategy (most recent by mtime in parent? env var? CLI flag?).
  - C6.3: copy assembled FASTAs / plasmids / MGEs into `repositorio/`. Reads as symlinks (size). Decide: naming convention to avoid collisions across runs.
  - C6.4: alert system. Triggers on PTU / MOB / MPF / Rep / AMR-gene / virulence-gene matches against `data_analysis.tsv`. Decide: match criteria (exact tuple? fuzzy?), output channel (`alerts.tsv`? HTML section? terminal?), what to show (which prior runs/samples matched, dates).

### State of the working tree at end of session

- **Branch:** main
- **Staged:** none (all changes applied via `Edit` to existing tracked files; new `.gitignore` is untracked).
- **Untracked new files:** `.gitignore`, `examples/list_seq_template.tsv` (from prior session).
- **Deleted** (via `git rm -f`): `envs/aluminion_base.yml`, `envs/aluminion_copla.yml`, `envs/environment.yml`, `scripts/__pycache__/*`.
- **No commit was made** in this session — user has not yet asked for it. Run a verification round first, then commit per the user's "vamos por partes" convention.

### Immediate next step on resume

1. `mamba env update -f envs/aluminion_assembly.yml --prune` to pick up blast/mafft/emboss/pandas/matplotlib.
2. Verify polished assemblies aren't empty: `ls -la 03_assemblies/<sample>/assembly.fasta`.
3. `aluminion -r <run> -b … -l … --resume` and confirm deconcat → circlator → QUAST → annotation finishes cleanly.
4. Once green, decide whether to commit the bug-fix tranche as a single commit or split.
5. Resume the task list with Tier A items.
