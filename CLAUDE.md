# Bioinformatics Analysis Projects — Lab Context

## Antes de comenzar
Antes de rellenar huecos de información avisa de la información faltante. 
Lee los ficheros del entorno de trabajo. Si encuentras problemas a la hora de leer o escribir en un fichero (por ejemplo por problemas de permisos) avísame antes de continuar el razonamiento para que yo lo solucione antes de continuar. 
Comprueba que las funciones de abrir y escribir ficheros están refiriendo a archivos que están presentes en el repo, si faltan, haz una lista de los ficheros que te hacen falta para trabajar. 
Si hay que realizar imports o se llama a herramientas que no sean propias del sistema, genera un archivo YML que permita ejecutar el código creando un entorno conda.
Genera un README.md que explique cuál es el objetivo general del proyecto, qué archivos abre, qué archivos genera, cómo instalarlo, qué otras herramientas o bases de datos hace falta instalar, dónde deberían guardarse (si hay direcciones hardcodeadas, que no debería ser así), qué entornos o containers son necesarios y cómo instalarlos.
Genera un plan de trabajo antes de comenzar a escribir código en el que queden claras cuáles son las prioridades: 
Clasificación de scripts por nivel de funcionamiento: qué scripts ya funcionan a la perfección, cuáles de ellos tienen pequeños errores de estilo que se podrían mejorar para utilizar menos librerías o hacer el código más simple y legible, qué scripts es posible que no funcionen o directamente aquellos que no podrían funcionar de ningún modo por tener errores evidentes de código. 
Valoración de tareas por prioridad: Organiza el trabajo de manera que nos centremos en que la maquinaria principal funcione para que podamos añadir más adelante nuevos scripts que aporten funciones extra (como plots con R, tablas extra, o mejorar el README) cuando ya tengamos un suelo seguro sobre el que trabajar.
Genera código comentado en inglés americano técnico porque el público destino generalmente serán otros científicos de todo el mundo y queremos que pueda ejecutarlo cualquiera. Esto aplica también para las funciones de Python que vayamos definiendo, para el README, para las cabeceras de todas las tablas y títulos/ejes de los plots de todo lo que vayamos generando, etc.
Cuando termines de realizar cambios en el código o los READMEs genera un texto para el commit a GitHub como por ejemplo: git commit -m "Added new --resume flag to continue stopped executions" o git commit -m "Added logs and improved preprocessing steps". Puedes extenderte más si quieres. 


## Overview

These repositories contain pipelines and scripts for two main types of bioinformatic analyses:

1. **Bacterial genomics** — primarily Oxford Nanopore long-read sequencing data
2. **16S rRNA metataxonomy** — microbiome analysis from clinical human samples

All code must be written and commented in **English**. Variable names, function names,
and documentation must also be in English. All plots must be publication-ready (high
resolution, clean aesthetics, suitable for international peer-reviewed journals).

---

## Project type 1: Bacterial genomics

### Overview

Whole-genome sequencing analysis of bacterial isolates, primarily using **Oxford Nanopore
Technology (ONT) long reads**, with occasional Illumina short reads for hybrid assembly.

### Typical pipeline steps

| Step | Tools (preferred → alternative) |
|------|----------------------------------|
| Preprocessing / QC | NanoPlot, Chopper |
| Assembly | Flye  |
| Assembly polishing | dorado / deconcat |
| Assembly QC | QUAST, CheckM / Busco / Bandage |
| Annotation | Prokka / Bakta |
| MLST typing | mlst (Torsten Seemann) / PubMLST API |
| Species/subspecies detection | Kleborate (Klebsiella), ECTyper (E. coli) |
| AMR detection | ABRicate, RGI (CARD) / ResFinder |
| Virulence genes | ABRicate (VFDB), BLAST |
| Plasmid detection | Platon, MOB-suite / copla |
| Integron detection | IntegronFinder |
| Prophage detection | PHASTER/PHASTEST, PhiSpy / Vibrant |
| Pan-genome | Roary / Panaroo / chewbacca |
| Phylogenetics | IQ-TREE2, FastTree / Snippy (SNP-based) |

### Technology stack

- **Primary language**: Bash (shell scripts) with output parsing via Python
- **Future direction**: Nextflow (actively migrating — add Nextflow versions when possible)
- **Containerization**: Docker preferred over Conda for new tools; Conda still in use for legacy tools
- **Avoid**: hardcoded absolute paths — use variables or arguments for all paths

### Bash scripting conventions

```bash
#!/usr/bin/env bash
set -euo pipefail   # Always: exit on error, undefined vars, pipe failures

# --- Constants and paths (from arguments, never hardcoded) ---
readonly THREADS="${THREADS:-8}"
readonly INPUT_DIR="${1:?Usage: script.sh <input_dir> <output_dir>}"
readonly OUTPUT_DIR="${2:?Usage: script.sh <input_dir> <output_dir>}"

# --- Logging ---
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
log "Starting assembly pipeline for sample ${SAMPLE}"
```

- Always use `set -euo pipefail` at the top of every script
- All scripts must accept input/output as positional arguments — no hardcoded paths
- Use a `log()` function for timestamped logging to stderr
- Validate that input files exist before processing (`[[ -f "$FILE" ]] || { log "ERROR: ..."; exit 1; }`)
- Clean up temporary files using `trap cleanup EXIT`
- Comment each major step explaining *why*, not just *what*

### Nextflow conventions (for new pipelines)

- Use DSL2 syntax exclusively
- Each process in its own `.nf` file under `modules/local/` or `modules/nf-core/`
- Use Docker containers per process via the `container` directive
- All parameters defined in `nextflow.config` and overridable via `params.yaml`
- Include `publishDir` with `mode: 'copy'` for all final outputs
- Follow nf-core module structure when possible for reusability

---

## Project type 2: 16S rRNA metataxonomy

### Overview

Microbiome analysis from clinical human samples (nasopharyngeal, intestinal, stool,
vaginal, etc.) using amplicon sequencing of the 16S rRNA gene.

### Sequencing technologies and target regions

| Technology | 16S region | Typical use case |
|------------|------------|------------------|
| Illumina | V3-V4 | Most projects, high sample throughput |
| PacBio | Full-length 16S | Higher taxonomic resolution projects |
| Oxford Nanopore | Variable | Emerging, case-by-case basis |

### Bioinformatics stack

- **Denoising**: QIIME2 with DADA2 plugin
- **Taxonomy database**: SILVA (weighted classifier, latest available version)
- **Downstream analysis**: R (primary language for all statistics and visualization)

#### Key R packages

```r
# Core microbiome data handling
library(qiime2R)       # Import QIIME2 artifacts into R
library(phyloseq)      # Core microbiome data structure and methods
library(microbiome)    # Utilities and transformations for phyloseq objects

# Diversity and community ecology
library(vegan)         # Diversity indices, ordination, permanova

# Statistics
library(lme4)          # Linear mixed models
library(nlme)          # Nonlinear mixed effects models
library(emmeans)       # Marginal means and pairwise contrasts
library(rstatix)       # Tidy-friendly statistical tests

# Visualization
library(ggplot2)       # All figures
library(ggpubr)        # Publication-ready ggplot2 extensions
library(patchwork)     # Combining multiple ggplot panels
library(RColorBrewer)  # Color palettes
library(viridis)       # Colorblind-friendly continuous palettes
library(ggsignif)      # Significance brackets on ggplot figures

# Data manipulation
library(tidyverse)     # dplyr, tidyr, readr, stringr, purrr, forcats

# Path management
library(here)          # Reproducible relative paths (always use here::here())
```

### Interactive pipeline design — mandatory

**These pipelines must interrogate the user** before running key steps. Each pipeline
script must present interactive menus for the decisions below. Never assume defaults
silently — always show what the default is and ask for confirmation.

Decisions that require user input:

- **Taxonomic level**: Genus / Family / Order / Class / Phylum (or combinations)
- **Clinical variable selection**: show available metadata columns, let user pick which
  to use for grouping, coloring, and statistical comparisons
- **Filtering thresholds**: minimum read count per sample, minimum ASV prevalence
- **Rarefaction**: apply or not, at what sequencing depth
- **Alpha diversity metrics**: Shannon, Simpson, Chao1, Observed, Faith's PD, etc.
- **Beta diversity metrics**: Bray-Curtis, weighted/unweighted UniFrac, Jaccard
- **Ordination method**: PCoA, NMDS, t-SNE, UMAP
- **Statistical framework**: parametric vs non-parametric (check normality first and
  suggest the appropriate test)

#### Standard interactive menu pattern

```r
# Present a numbered menu with clear descriptions
prompt_menu <- function(title, options, default = 1) {
  cat(sprintf("\n=== %s ===\n", title))
  for (i in seq_along(options)) {
    marker <- if (i == default) " [default]" else ""
    cat(sprintf("  %d. %s%s\n", i, options[i], marker))
  }
  choice <- readline(prompt = sprintf("Enter choice [1-%d] (Enter for default): ",
                                      length(options)))
  if (nchar(trimws(choice)) == 0) return(default)
  choice <- suppressWarnings(as.integer(choice))
  if (is.na(choice) || choice < 1 || choice > length(options)) {
    message(sprintf("Invalid input. Using default: %s", options[default]))
    return(default)
  }
  return(choice)
}

# Usage example
tax_choice <- prompt_menu(
  title   = "Select taxonomic level for analysis",
  options = c("Genus", "Family", "Order", "Class", "Phylum"),
  default = 1
)
tax_level <- c("Genus", "Family", "Order", "Class", "Phylum")[tax_choice]
message(sprintf("Proceeding with taxonomic level: %s", tax_level))
```

### R coding conventions

- All code and all comments in **English**
- Use `tidyverse` style: native pipe `|>`, tidy data principles, no `$` chaining
- Never use `attach()`, `setwd()`, or `source()` with absolute paths
- Use `here::here()` for every file path without exception
- Each script must start with a structured header:

```r
# =============================================================================
# Script:  alpha_diversity_analysis.R
# Purpose: Compute and visualize alpha diversity metrics; test associations
#          with clinical variables via linear mixed models
# Input:   data/phyloseq_filtered.rds
#          data/metadata/clinical_metadata.csv
# Output:  figures/alpha_diversity_boxplots.pdf
#          results/diversity/alpha_stats_summary.csv
# Author:  [Lab name / PI]
# Date:    YYYY-MM-DD
# =============================================================================
```

### Publication-quality plots — mandatory standards

All figures must be exportable at ≥300 dpi as PDF (vector) and PNG (raster preview).

```r
# Reusable theme — apply to every ggplot in the project
theme_publication <- function(base_size = 12, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(colour = "grey92", linewidth = 0.3),
    panel.border       = element_rect(colour = "black", linewidth = 0.8),
    axis.ticks         = element_line(colour = "black", linewidth = 0.5),
    axis.text          = element_text(colour = "black", size = base_size - 2),
    axis.title         = element_text(colour = "black", size = base_size,
                                      face = "bold"),
    legend.background  = element_blank(),
    legend.key         = element_blank(),
    legend.title       = element_text(face = "bold"),
    strip.background   = element_rect(fill = "grey92", colour = "black",
                                      linewidth = 0.8),
    plot.title         = element_text(size = base_size + 1, face = "bold"),
    plot.subtitle      = element_text(colour = "grey40", size = base_size - 1)
  )
}

# Standard save function — always export both formats
save_publication_plot <- function(plot, filename, width = 7, height = 5) {
  ggsave(
    filename = here("figures", paste0(filename, ".pdf")),
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = 300,
    device   = cairo_pdf   # preserves vector fonts
  )
  ggsave(
    filename = here("figures", paste0(filename, ".png")),
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = 300
  )
  message(sprintf("Saved: figures/%s (.pdf + .png)", filename))
}
```

**Color guidelines:**
- Use colorblind-safe palettes by default: `viridis`, `RColorBrewer "Set2"` or `"Paired"`
- Never use red/green combinations
- Use consistent color assignments for the same clinical groups across all figures
  in the same project (define a named color vector at the top of each script)

---



**Critical rules:**
- `data/raw/` is strictly read-only — never write there under any circumstance
- All outputs go to `results/` or `figures/`
- Never commit raw sequencing data: `.fastq.gz`, `.fast5`, `.pod5`, `.bam` are in `.gitignore`

---

## Literature and repository search — apply proactively

When implementing a new analysis step, choosing parameters, or selecting a tool,
**always search for prior work** before writing code. This is not optional.

### What to search

1. **GitHub** — existing pipelines with similar tool combinations
   - Search for: workflows combining the same tools, Nextflow nf-core modules,
     Snakemake workflows (e.g., snakemake-workflows organization)
   - Note: parameter choices, filtering thresholds, output formats used by others

2. **PubMed / Google Scholar** — methods validation and benchmark papers
   - Check if the chosen tool/parameter combination is validated in the literature
   - Search for benchmarking studies (e.g., assembler comparisons for ONT data,
     DADA2 vs Deblur for 16S, classifier benchmarks for SILVA vs GTDB)
   - Look for recent publications (last 2 years) suggesting improved approaches

### What to report after searching

After any literature/repository search, always summarize:

- **What exists**: similar pipelines or published methods and their key design decisions
- **Alignment check**: whether our current approach aligns with or deviates from
  established best practices, and why any deviation is justified
- **Gaps**: parameters, QC steps, or statistical tests that others consistently include
  and that we might be missing
- **Recent developments**: any tools or methods published in the last 1-2 years that
  could improve our results (even if we decide not to adopt them, flag them)

This is especially important for: filtering thresholds for ONT reads, assembly
parameters (Flye genome size and mode), DADA2 truncation lengths, SILVA classifier
training parameters, rarefaction decisions, and choice of diversity metrics.

---

## Code quality — universal rules (both project types)

- **Comments**: every non-obvious block must have a comment explaining *why*, not *what*
- **No magic numbers**: all thresholds and parameters go in named constants or config files
- **Error handling**: validate inputs, check expected file formats, fail with clear messages
- **Reproducibility**: log all software versions at the start of every pipeline run
- **Modularity**: small, single-purpose functions and scripts over monolithic ones
- **No silent failures**: every tool call must check exit status (Bash: `set -e`; R: check return values)

---

## Project-specific notes

### Aluminion — ONT bacterial WGS pipeline

Target organisms: Enterobacteriaceae (K. pneumoniae, E. coli, etc.).
Sequencing: Oxford Nanopore MinION / Mk1D, R10.4.1 flowcells, basecalled with Dorado (sup model).
Pipeline language: Bash orchestrator (`aluminion.sh`) + Python parsers (`scripts/`).

---

### Repository layout

```
Aluminion/                        ← git root
├── CLAUDE.md                     ← AI coding context and lab conventions (this file)
├── README.md                     ← Full user-facing documentation
├── aluminion.sh                  ← Main pipeline orchestrator
├── install.sh                    ← Installation (conda envs, Docker images, databases)
│
├── envs/                         ← One conda environment per tool group
│   ├── aluminion_reads.yml       ← NanoPlot, Chopper, pillow, kaleido (>=1.0.0)
│   ├── aluminion_assembly.yml    ← Flye, QUAST, dorado (binary in PATH)
│   ├── aluminion_circlator.yml   ← dnaapler (env name retained; legacy circlator was EOL)
│   ├── aluminion_annot.yml       ← Bakta, ABRicate, BLAST, MOB-suite, GAMBIT, mlst, ECTyper
│   ├── aluminion_integron.yml    ← IntegronFinder
│   └── aluminion_kleborate.yml   ← Kleborate
│
├── scripts/                      ← Python output parsers
│   ├── parser.py                 ← Main aggregator — called at the end of aluminion.sh
│   │                               Accepts --skip-kraken / --skip-abr / --skip-typing /
│   │                               --skip-phages flags; internally calls copla_parser.py
│   ├── copla_parser.py           ← Copla plasmid classification parser (called by parser.py)
│   ├── IS_parser.py              ← BLAST vs ISfinder output → IS_chr_out.tsv per sample
│   ├── integron_parser.py        ← IntegronFinder output → integron_summary.csv
│   ├── phage_parser.py           ← Phastest output → phage_summary.csv
│   ├── lab_db_updater.py         ← Builds data_seq.tsv / data_analysis.tsv (cumulative lab DB)
│   ├── aluminion_reporter.py     ← Generates Aluminion_Report.html from final tables
│   └── deconcat.py               ← Deconcatenates reads for dorado polish step
│
├── docs/
│   └── Dorado Polish Documentation.html   ← Reference docs for dorado polish
│
├── examples/                     ← Real output files from a completed run (reference)
│   ├── list_seq.tsv              ← Sample metadata input format
│   ├── data_seq.tsv / data_analysis.tsv
│   ├── QC_reads.csv / QC_assembly.csv
│   ├── taxonomy.csv / taxonomy.xlsx / kraken.csv / kraken_report.csv
│   ├── mlst.csv / mlst_modif.csv / kleborate.tsv
│   ├── AbR_report.csv / AbR_modif.xlsx
│   ├── IS.tsv / integron_summary.csv / phage_summary.csv
│   ├── copla.txt / copla_modif.csv
│   ├── gambit.csv / genus.csv / species.csv / enterobacterales__species_output.txt
│   └── output.tsv / kraken_mlst.xlsx
│
├── test/                         ← Partial test run (mirrors run-time folder layout)
│   ├── list_seq.tsv
│   ├── 01_reads/QC/              ← NanoPlot pre-filter outputs
│   ├── 02_filter/QC/             ← NanoPlot post-filter outputs
│   ├── 03_assemblies/quast/      ← QUAST assembly QC
│   ├── 04_taxonomies/
│   │   ├── kraken2/              ← Kraken2 reports
│   │   ├── gambit.csv
│   │   ├── kleborate/
│   │   └── ectyper/
│   ├── 08_Anotacion/AbR.tab      ← ABRicate AMR output
│   └── [final output CSVs/XLSXs]
│
└── tests/
    └── test_parser.py            ← Unit tests for parser.py
```

---

### Run-time folder layout (created by aluminion.sh)

The pipeline is called from a **parent directory** containing `list_seq.tsv`.
Each run creates a subdirectory `RUN_NAME/`:

```
parent_dir/
├── list_seq.tsv                         ← Sample sheet (tab-separated; must be here)
└── RUN_NAME/                            ← WORKDIR; created by aluminion.sh
    ├── samples                          ← Plain text list of sample IDs (generated)
    ├── aluminion_YYYYMMDD_HHMMSS.log    ← Full run log (stdout + stderr via tee)
    │
    ├── 01_reads/
    │   ├── {sample}.fastq.gz            ← Concatenated raw reads per sample
    │   └── QC/{sample}/
    │       ├── NanoStats.txt            ← Resume sentinel for pre-filter NanoPlot
    │       ├── LengthvsQualityScatterPlot_loglength_kde.png
    │       └── *.html *.png             ← Other NanoPlot outputs
    │
    ├── 02_filter/
    │   ├── {sample}.fastq.gz            ← Chopper-filtered reads (resume sentinel)
    │   └── QC/{sample}/
    │       ├── NanoStats.txt            ← Resume sentinel for post-filter NanoPlot
    │       └── *.html *.png
    │
    ├── 03_assemblies/{sample}/
    │   ├── assembly.fasta               ← Final polished + circularized assembly
    │   ├── assembly_graph.gfa           ← Flye assembly graph
    │   ├── .polished                    ← Sentinel: dorado polish complete
    │   └── .circlator_done             ← Sentinel: dnaapler reorientation complete (name retained)
    │
    ├── 04_taxonomies/
    │   ├── kraken2/{sample}.report      ← Resume sentinel for Kraken2
    │   ├── gambit.csv                   ← GAMBIT species typing (all samples)
    │   ├── mlst/{sample}.tsv
    │   ├── kleborate/{sample}.tsv
    │   └── ectyper/{sample}/
    │
    ├── 05_IS/{sample}/
    │   └── IS_chr_out.tsv               ← BLAST vs ISfinder (resume sentinel)
    │
    ├── 06_integrons/{sample}/           ← IntegronFinder output directory
    │
    ├── 07_phages/{sample}/              ← Phastest output directory
    │
    ├── 08_Anotacion/{sample}/
    │   ├── {sample}.gbff                ← Bakta annotation (resume sentinel)
    │   ├── {sample}.faa / .gff3 / .tsv ← Bakta accessory outputs
    │   ├── mob_recon/                   ← MOB-suite plasmid reconstruction
    │   ├── AbR.tab                      ← ABRicate AMR per-sample results
    │   └── copla/                       ← Copla plasmid classification
    │
    └── [final output tables — generated by parser.py + aluminion_reporter.py]
        ├── data_seq.tsv                 ← Sequencing run metadata
        ├── QC_reads.csv                 ← Per-sample read QC summary
        ├── QC_assembly.csv             ← Per-sample assembly QC (QUAST)
        ├── taxonomy.csv / .xlsx         ← Species typing (GAMBIT + Kraken2)
        ├── kraken.csv / kraken_report.csv
        ├── mlst.csv / mlst_modif.csv
        ├── kleborate.tsv
        ├── AbR_report.csv / AbR_modif.xlsx
        ├── IS.tsv
        ├── integron_summary.csv
        ├── phage_summary.csv
        ├── copla.txt / copla_modif.csv
        ├── data_analysis.tsv            ← Master table (all results merged)
        └── Aluminion_Report.html        ← Interactive HTML summary report
```

---

### Pipeline stages

| Stage | Folder | Tools | Conda env | Skip flag |
|-------|--------|-------|-----------|-----------|
| 1 — Read QC & filtering | `01_reads/`, `02_filter/` | NanoPlot, Chopper | `aluminion_reads` | `--skip-preprocessing` |
| 2 — Assembly & polishing | `03_assemblies/` | Flye, dorado polish, dnaapler | `aluminion_assembly`, `aluminion_circlator` | — |
| 3 — Annotation & AMR | `08_Anotacion/` | Bakta, ABRicate, MOB-suite, Copla | `aluminion_annot` | `--skip-abr` (ABRicate only) |
| 4 — Taxonomy & typing | `04_taxonomies/` | Kraken2, GAMBIT, mlst, Kleborate, ECTyper | various | `--skip-kraken`, `--skip-typing` |
| 4 — MGEs | `05_IS/`, `06_integrons/`, `07_phages/` | BLAST/ISfinder, IntegronFinder, Phastest | `aluminion_annot`, `aluminion_integron` | `--skip-phages` |
| 5 — Consolidation | RUN_NAME root | parser.py, aluminion_reporter.py | `aluminion_annot` | — |

Docker images required: `kbessonov/mob_suite:3.0.3`, `rpalcab/copla:1.0`, phastest (local compose).

---

### Key implementation decisions (non-obvious — read before editing)

- `set +e` / `set -e` brackets around NanoPlot loops: choreographer (NanoPlot's Plotly
  dependency) spawns child browser processes that become direct bash children; `wait`
  catches their non-zero exit under `set -e`, hanging the pipeline. The brackets + `& done; wait`
  pattern parallelises NanoPlot across samples while ignoring choreographer exit codes.
- `MPLBACKEND=Agg PLOTLY_RENDERER=kaleido` on all NanoPlot calls: forces non-interactive
  matplotlib backend and routes Plotly static exports through kaleido, bypassing
  choreographer/Chrome entirely on headless servers.
- `kaleido` is NOT available on conda-forge or bioconda — must be installed via pip with no
  version pin. NanoPlot 1.46.2 requires kaleido>=1.0.0; pip installs the 1.x Go-based binary
  (~small, fast). The old 0.2.1 bundled a full Chromium (~200 MB) and hung during install.
- `readlink -f` is applied to `SEQ_LIST_INPUT` and `BASE_DIR` immediately after argument
  parsing, before any `cd "$WORKDIR"`, to prevent relative path breakage.
- Resume sentinels: `.polished` and `.circlator_done` are `touch`-created files because
  both steps overwrite `assembly.fasta` (can't use the assembly as its own sentinel).
  Note: `.circlator_done` is kept as the sentinel name even though the tool is now
  dnaapler (legacy circlator was EOL with broken libcrypto.so.1.0.0 dep).
- Flye failure: interactive 3-choice menu (skip sample / retry with `--meta` / stop pipeline).
  Skipped samples are removed from the `samples` file so downstream loops ignore them.
- `copla.txt` is only truncated (`> copla.txt`) on a fresh run, not on `--resume`, to
  allow appending results for samples that weren't done yet.
- **Polishing — two-step @RG injection**: `dorado polish` enforces two checks on the
  input BAM: (1) `@PG` must contain `PN:dorado` — meaning `dorado aligner` must be used,
  **not minimap2** (minimap2 fails this check with "Input BAM file was not aligned using
  Dorado"); (2) `@RG` must contain `DS:basecall_model=<model>`. `dorado aligner` from
  FASTQ input writes the `PN:dorado` `@PG` but does NOT inject `@RG DS:basecall_model`.
  Fix: extract `basecall_model_version_id` from the FASTQ read header, align with
  `dorado aligner`, then inject `@RG\tID:1\tDS:basecall_model=<model>` via
  `samtools addreplacerg`. Two FASTQ formats handled: with `RG:Z:` tag → use
  `--add-fastq-rg` + `addreplacerg -w` (overwrite); without `RG:Z:` (current standard
  Dorado output) → `dorado aligner` alone + `addreplacerg`. Neither `--ignore-read-groups`
  nor `--device cpu` are used: the former was a workaround for missing @RG (now fixed),
  the latter prevented GPU use (dorado now auto-detects).
