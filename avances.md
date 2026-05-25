# Aluminion — Resume-from-scratch handoff (Claude-only working doc)

> This file is **for Claude Code only** — a self-contained briefing to resume work
> in a fresh session after the machine is powered off. It is **not** pushed to GitHub.
> The canonical repo/runtime layout and "do-not-revert" decisions live in `CLAUDE.md`;
> this file adds the data-flow detail and the list of **open** work. Resolved bugs
> from earlier sessions have been pruned — only what still matters is kept.

Last updated: 2026-05-25 (after merging the MGE alert system into `main`).

---

## 1. What Aluminion is (scope)

A fully automated, modular pipeline for **bacterial whole-genome sequencing (WGS)**
from Oxford Nanopore (MinION / Mk1D) reads. It goes from raw `fastq_pass/` to an
interactive HTML report, and maintains **cumulative lab databases**
(`data_seq.tsv`, `data_analysis.tsv`) plus a persistent **MGE repository** that grow
across runs for longitudinal surveillance.

- **Target organisms:** Enterobacteriaceae (*Klebsiella pneumoniae*, *E. coli*,
  *Enterobacter*, *Citrobacter*, …). Assembly/AMR/annotation work for any bacterium
  with a published MLST scheme.
- **Use case:** serial surveillance in a hospital micro lab; tracking AMR + mobile
  genetic elements (MGEs) across isolates and over time.
- **Design constraints:** partial runs (`--just-*`, `--skip-*`), resumable
  (`--resume`), non-fatal optional steps (a failed polish/typing/MGE step only warns),
  headless-Linux-server friendly (no display, GPU optional).
- **Orchestration:** Bash (`aluminion.sh`) drives CLI tools across 6 conda envs +
  3 Docker images; Python scripts in `scripts/` parse tool outputs into tables.
- **Future direction:** Nextflow DSL2 migration (not started).

---

## 2. Workflow — the 7 stages

`aluminion.sh` is called from a **parent dir** holding `list_seq.tsv`; it creates a
`RUN_NAME/` subfolder and works inside it. (Full runtime tree + sentinel files are in
`CLAUDE.md`.)

| Stage | What | Tools | Conda env(s) | Skip flag |
|------:|------|-------|--------------|-----------|
| 1 | Read QC + filtering | NanoPlot (pre), Chopper, NanoPlot (post) | `aluminion_reads` | `--skip-preprocessing` |
| 2 | Assembly | Kraken2, Flye, dorado polish, deconcat, dnaapler | `aluminion_assembly`, `aluminion_circlator` | — |
| 3 | Assembly QC | QUAST, Bandage | `aluminion_assembly` | — |
| 4 | Annotation, typing, AMR, MGEs | Bakta, MOB-suite, Abricate (AMR + VFDB), GAMBIT, MLST, Kleborate, ECTyper, Copla, IntegronFinder, Phastest, ISfinder BLASTn | `aluminion_annot`, `aluminion_integron`, `aluminion_kleborate` | `--skip-abr`, `--skip-typing`, `--skip-plasmids`, `--skip-integrons`, `--skip-phages` |
| 5 | Consolidation | parser.py → aluminion_reporter.py → lab_db_updater.py | `aluminion_annot` | — |
| 6 | Cross-run MGE alerts | mge_repository.py + mge_alerts.py → alerts_reporter.py | `aluminion_annot` | `--no-alerts` |

Early-stop: `--just-preprocessing` (after Chopper), `--just-assembly` (after QUAST).
`--init-db` builds the cumulative DBs from scratch; `--init-repo` builds the MGE repo.
`--no-minknow` / `--unique-run` control MinKNOW copy and `../repositorio/` interaction.
Bakta and Flye are the only non-skippable / fatal steps; everything else is non-fatal.

---

## 3. Scripts and what they do

All in `scripts/`. `parser.py` is the consolidation orchestrator and **imports**
`phage_parser`, `integron_parser`, `copla_parser` (calls their `run_parsing()`),
so those three must NOT be invoked separately from bash (doing so double-parses).

| Script | Role |
|--------|------|
| `parser.py` | Master consolidator → `taxonomy.csv/.xlsx`, `AbR_modif.xlsx`, `VF_modif.xlsx`, `mlst_modif.csv`, `kraken_mlst.xlsx`. Imports the 3 sub-parsers. |
| `copla_parser.py` | `copla.txt` → `copla_modif.csv` (one row per plasmid). |
| `integron_parser.py` | IntegronFinder output → `integron_summary.csv`. |
| `phage_parser.py` | Phastest `summary.txt` → `phage_summary.csv`. |
| `IS_parser.py` | BLASTn vs ISfinder per sample → `IS_chr_out.tsv` (aggregated into `IS.tsv` by bash). |
| `deconcat.py` | Deconcatenate a circular assembly before dnaapler reorientation. |
| `aluminion_reporter.py` | Builds `Aluminion_Report.html` from the consolidated tables. |
| `lab_db_updater.py` | Updates cumulative `data_seq.tsv` / `data_analysis.tsv`. (Was `Datos_seq_unified2.py`.) |
| `mge_repository.py` | Persistent MGE repository: indices + plasmid FASTAs + integron gene-sets. `Repository.init(root)`, ingestion, ANI/Jaccard matching. |
| `mge_alerts.py` | Ingests this run's plasmids/integrons into the repo, matches vs history, emits `alerts.tsv`. CLI: `--run-dir --repo --run-name [--alert-new-priority] [--no-ingest] [--ani-threshold --jaccard-threshold --size-tolerance --min-plasmid-size]`. |
| `alerts_reporter.py` | `Alerts_Report.html` from alert records. |
| `_priority_genes.py` | Curated priority resistance/virulence gene catalog (carbapenemases, mcr, hypervirulence). `classify_priority()`. |
| `_log.py` | Shared stderr logger (colour on TTY, respects `NO_COLOR`). Used by parser/integron/copla parsers. |

---

## 4. Data-table transformations (how reports get built)

### Bash pre-processing (in `aluminion.sh`, stage 5 head)
- `04_taxonomies/kraken2/<s>.report` → per-sample awk filters genus (`$4=="G"`,
  ≥20%) and species (`$4=="S"`, ≥4%) rows, emits `<sample>\t<percent>\t<taxon>` with
  **`-v OFS='\t'` + manual name concatenation** (do NOT revert to comma-print — that
  reintroduces the whitespace bug that emptied taxonomy columns).
- NanoStats → `QC_reads.csv` (pre+post pasted); QUAST `transposed_report.tsv` →
  `QC_assembly.csv` with header column **`Samples`** (plural — both reporter and
  lab_db_updater key on that exact name).
- Per-sample IS counts → `IS.tsv`.

### `parser.py` (the core merge)
Reads each tool output via `safe_read_csv` (returns empty DF with required cols on
missing/empty file), reshapes, then left-joins everything on **`Sample`**:
- **Abricate** `AbR_report.csv`: per-gene presence matrix → single `Resistance_genes`
  comma-string per sample → `AbR_modif.xlsx`.
- **VFDB** `VF_report.csv`: same shape → `Virulence_genes` → `VF_modif.xlsx`
  (consumed by `mge_alerts.py`, NOT merged into data_analysis.tsv).
- **Kraken** genus/species csv → `Majority_genus`, `Majority_species`,
  `Contaminants` (first row = majority, rest packed comma-joined).
- **GAMBIT** `gambit.csv`: `query`→`Sample`, `closest.description`→`Subspecies`
  (bracket/paren qualifiers stripped).
- **MLST** `mlst.csv` → `parse_mlst()` builds a fixed-13-column DataFrame
  (`MLST_HEADER`); when ST is unresolved it looks up candidate STs/alleles from the
  PubMLST schema (`--pubMLST_database`). Sample = basename of the FASTA path.
- **Kleborate** → `KO_locus`, `Carbapenemase`, `ESBL`, `Other_resistance`,
  `N_AMR_genes`, `AMRscore`, `VIRscore` (version suffixes `.v1^`/`^` stripped, deduped).
- **ECTyper** → `Serotype`.
- Final: merge kraken→gambit→mlst→kleborate→ectyper on `Sample` → `final_columns` →
  `taxonomy.csv` + `taxonomy.xlsx`; a slimmer view → `kraken_mlst.xlsx`.

### `aluminion_reporter.py`
Loads `QC_reads.csv`, `QC_assembly.csv` (key `Samples`), `taxonomy.xlsx`,
`AbR_modif.xlsx`, `copla_modif.csv`, `integron_summary.csv`, `phage_summary.csv`,
`kleborate.tsv` → multi-card Bootstrap/DataTables HTML (`Aluminion_Report.html`):
QC card, Taxonomy & AMR card, Plasmids card, Integrons card, Phages card.

### `lab_db_updater.py`
Reads `list_seq.tsv`, `QC_reads.csv`, `QC_assembly.csv`, `taxonomy.csv`, MinKNOW
`final_summary_*.txt` / `report_*.json` (optional), `copla_modif.csv`,
`phage_summary.csv`, `integron_summary.csv`. Computes Depth, assembly score, MGE
counts; merges on **`ID`** into cumulative `data_seq.tsv` (sequencing/QC) and
`data_analysis.tsv` (taxonomy/AMR/MGE). On first run (or `--init`) creates them;
otherwise writes `*_new.tsv` for review. CLI cutoffs: `--extraction-kit`,
`--depth-threshold`.

### `mge_alerts.py` + `mge_repository.py` (stage 6)
Ingests this run's plasmids (`08_Anotacion/<s>/mob_recon/plasmid_<cluster>.fasta` +
`copla_modif.csv`) and integrons (`integron_summary.csv`) into the repo at `--repo`
(default `$BASE_DIR/repository`). Matches against history: **plasmids by skani ANI
≥99%** (tuple fallback on PTU/Rep/MOB when no assembly), **integrons by Jaccard ≥0.8**
over cassette gene-sets. Host metadata joined from `data_analysis.tsv` (by `ID`) +
`taxonomy.csv` (by `Sample`) + per-sample virulence from `VF_modif.xlsx`. Emits
`alerts.tsv` + `Alerts_Report.html`. `--alert-new-priority` also flags
first-occurrence MGEs carrying a priority gene.

---

## 5. Environments & external dependencies

| Conda env | Key tools |
|-----------|-----------|
| `aluminion_reads` | NanoPlot, Chopper, pillow, kaleido (pip) |
| `aluminion_assembly` | Kraken2, Flye, QUAST, Bandage, samtools, minimap2, blast, mafft, emboss, pandas, matplotlib |
| `aluminion_circlator` | **dnaapler** (env name kept for back-compat; legacy circlator was EOL) + python=3.12 |
| `aluminion_annot` | Bakta, ABRicate, BLAST, MOB-suite, GAMBIT, mlst, ECTyper, datamash, pandas, openpyxl, biopython, bcbio-gff, **skani**, Python stack |
| `aluminion_integron` | IntegronFinder |
| `aluminion_kleborate` | Kleborate + ECTyper |

- Docker: `kbessonov/mob_suite:3.0.3`, `rpalcab/copla:1.0`, Phastest (local compose).
- External binary: `dorado` in `$PATH`. GPU optional.
- **Phastest runs as root inside the container** (its perl scripts are in
  `/root/...` with 700 perms); the call is wrapped in `bash -c "phastest …; chown -R
  <host_uid>:<host_gid> JOBS/<s>; exit $rc"` so outputs are deletable by the host user.
  Do NOT add `--user` or `--phage-only` — both silently break prophage detection.

---

## 6. Durable, non-obvious decisions (do not revert)

(Most are also in `CLAUDE.md`; the merge-specific ones are here.)

- `parser.py` whitespace handling: explicit per-column strip loops, **not**
  `DataFrame.apply(lambda x: x.str.strip())` — the apply idiom is silently a no-op in
  the installed pandas. Plus the OFS=tab awk above. Both defend the `Sample` merge key.
- `mlst_modif.csv` is now a padded 13-column DataFrame (`parse_mlst()`), not the old
  line-by-line writer; reporter coerces MLST with `pd.isna()` before `.endswith('.0')`.
- The cross-run MGE comparison is **only** `mge_alerts.py` + `mge_repository.py`. The
  old in-DB engine in `lab_db_updater.py` (`build_mge_table`/`find_shared_mges`,
  `data_mge.tsv`/`mge_shared.tsv`, `--alert-all-mge`) was **retired** in the merge.
- `scripts/__pycache__/*.pyc` are untracked (gitignored). Don't re-add.

---

## 7. Repo / git / machine state

- **Branch `main`** has the alert system merged in (`8fc0db4` integrate + `921f0dc`
  merge `--no-ff`). **Pushed** — as of 2026-05-25, `origin/main == main == 921f0dc`
  (verified via fetch). `avances.md`/`CLAUDE.md` rewrites are intentionally left
  uncommitted / not pushed.
- Rollback anchors (delete only after production validation): tag
  `premerge-main-20260525`, branch `backup/alert-system` (tip `71ff52a`, identical to
  the now-deleted `feature/alert-system` — preserves all 7 alert-system commits).
- Post-push cleanup **DONE** (2026-05-25): worktree `../aluminion-alerts` pruned,
  branches `integration/alert-merge` (was merged) and `feature/alert-system` (force,
  preserved by `backup`) deleted.
- **Still dangling:** locked worktree `worktree-agent-a8d9f6244aec47eb2` at
  `Proyectos/Aluminion/.claude/worktrees/` (a previous Claude agent's worktree, in a
  different path). Removal needs `git worktree remove --force`; left untouched pending
  user confirmation.
- **Two machines:**
  - *Production server* (`usuario@usuario-System-Product-Name`): full conda envs,
    Docker, dorado, databases. Real runs + end-to-end validation happen here.
  - *Desktop* (`DESKTOP-QSN0516`): two shells seen. Via Git Bash (MSYS2/MINGW64) the
    Windows miniconda `myenv` env runs code-level checks (Python 3.10.9, pandas 1.5.3,
    openpyxl, pytest 8.4.1, biopython) — this is what the Claude Bash tool reaches. A
    WSL side with the `aluminion_test` env (noted in earlier sessions) may also exist
    separately. No bioinformatics tools, Docker, or DBs here.
- Git identity on the WSL desktop set **locally** to `Aluminio-visto
  <jorgergrande@gmail.com>` (matches existing commits).

---

## 8. OPEN bugs and pending tasks

### 🔴 BUG — typing columns empty in taxonomy/data_analysis (highest priority)
`Subspecies` (GAMBIT) and `MLST` / `MLST_scheme` (mlst) come out **empty** in
`taxonomy.csv`, `data_analysis.tsv`, and the HTML report, even though `gambit.csv`
and `mlst.csv` contain the data — while `Majority_genus`/`Majority_species` (kraken,
same `Sample` join) ARE populated. Reproduced on **clean `main`** (pre-dates the alert
merge; not a regression) via the two failing tests:
`tests/test_parser.py::TestParserPy::test_taxonomy_key_columns_populated` (19/24
Subspecies empty) and `…TestReporterPy::test_html_contains_key_data` (`ecloacae`
MLST scheme absent). **Hypothesis:** `Sample`-key mismatch between the kraken-derived
base and the gambit/mlst frames (gambit `query` and mlst path-basename may differ in
form from the kraken sample names), and/or the `examples/gambit.csv` covers fewer
samples than `examples/species.csv`. **Not yet diagnosed/fixed.** Affects real output,
worth fixing. Start by printing the distinct `Sample` values of `kraken_df`,
`gambit_df`, `mlst_df` right before the merges in `parser.py`.

### 🟡 Pending validation
- **End-to-end run of the alert system on the production server** (skani installed via
  `conda env update -f envs/aluminion_annot.yml`; Docker/dorado/DBs present). Do a run
  with `--init-repo`, then a second run reusing the repo, confirming `VF_report.csv`,
  `VF_modif.xlsx`, `alerts.tsv`, `Alerts_Report.html`, `repository/index_plasmids.tsv`
  are produced; and that `--skip-abr` skips VFDB cleanly (mge_alerts tolerates missing
  `VF_modif.xlsx`). Only code-level tests have run so far (36/38; the 2 fails = the bug above).

### 🟢 Lower-priority / deferred
- Expose `mge_alerts.py` matching thresholds (`--ani-threshold`, `--jaccard-threshold`,
  `--size-tolerance`, `--min-plasmid-size`) as `aluminion.sh` flags if the lab wants to
  tune stringency.
- README pipeline mermaid (7 nodes) vs stage table (6 rows) group stages slightly
  differently — cosmetic mismatch, harmless.
- Migrate the bare `print()` calls in `aluminion_reporter.py` and `lab_db_updater.py`
  to `_log` (no ANSI there, just lots of prints). Low value.
- C4: extract `safe_read_csv` to `scripts/_utils.py` — declined (only parser.py uses it;
  YAGNI). Reconsider if another script needs defensive CSV loading.
- Nextflow DSL2 migration (long-term).

---

## 9. How to resume quickly

1. On the **production server**, `git pull origin main` once the user has pushed.
2. To re-run **only the consolidation** on an already-computed run (the most common
   iteration — see the saved Claude memory `reconsolidate-aluminion-run`): regenerate
   the kraken csvs, then `parser.py` → `aluminion_reporter.py` → `lab_db_updater.py`
   (+ `mge_alerts.py`) without relaunching the whole pipeline.
3. Tackle the typing-columns bug (§8) — it's the clearest open defect and has a
   reproducible test to fix against.
