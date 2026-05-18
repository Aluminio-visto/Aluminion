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

- **Primary language**: Bash (shell scripts)
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

## Directory structure (both project types)

```
project/
├── CLAUDE.md               ← this file
├── README.md               ← project overview, sample info, how to run
├── config/
│   ├── config.yaml         ← pipeline parameters and tool versions
│   └── samples.tsv         ← sample sheet: sample_id, condition, file paths
├── data/
│   ├── raw/                ← READ-ONLY — never modify original sequencing data
│   ├── processed/          ← intermediate outputs (trimmed, aligned, denoised)
│   └── metadata/           ← clinical variables and sample metadata tables
├── results/
│   ├── qc/                 ← quality control reports
│   ├── assembly/           ← (genomics) assembled genomes
│   ├── annotation/         ← (genomics) Prokka/Bakta output
│   ├── typing/             ← (genomics) MLST, Kleborate, ECTyper
│   ├── amr/                ← (genomics) AMR and virulence genes
│   ├── mge/                ← (genomics) plasmids, integrons, prophages
│   ├── phylogeny/          ← (genomics) trees and pan-genome
│   ├── diversity/          ← (metataxonomy) alpha and beta diversity
│   ├── taxonomy/           ← (metataxonomy) ASV tables, taxonomy assignments
│   └── differential/       ← (metataxonomy) differential abundance results
├── figures/                ← publication-ready plots (PDF + PNG)
├── scripts/
│   ├── bash/               ← shell scripts for genomics pipelines
│   ├── nextflow/           ← .nf modules and workflows
│   │   ├── modules/
│   │   └── workflows/
│   ├── r/                  ← R analysis scripts (metataxonomy and genomics stats)
│   └── python/             ← utility and helper Python scripts
├── notebooks/              ← exploratory R Markdown or Jupyter notebooks
├── docker/                 ← Dockerfiles for custom container images
├── envs/                   ← conda environment .yaml files (legacy tools)
├── tests/                  ← validation scripts for pipeline outputs
└── logs/                   ← timestamped log files (gitignored for large runs)
```

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

<!-- Add project-specific context here as the analysis develops. -->
<!-- Examples:

### Genomics project — K. pneumoniae outbreak 2024
- Organism: Klebsiella pneumoniae, clinical isolates, ICU outbreak
- Sequencing: MinION R10.4.1 flowcells, basecalled with Dorado v0.7.2
- Special: three samples are suspected co-infections — flag but do not exclude
- Reference: K. pneumoniae NTUH-K2044 (GenBank AP006725.1)
- Clinical metadata columns: sample_id, ward, date_collection, sequence_type, outcome

### Metataxonomy project — neonatal gut microbiome
- Sample type: stool, neonatal ICU patients
- Technology: Illumina MiSeq, V3-V4 region, paired-end 2×300 bp
- Clinical variables: gestational_age, delivery_mode, antibiotic_exposure, feeding_type
- Exclusion criteria: samples with < 5000 reads after DADA2 denoising
- Rarefaction depth: TBD after reviewing read count distribution
-->
