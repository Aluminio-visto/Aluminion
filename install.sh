#!/bin/bash
# ==============================================================================
# Aluminion — Installation script
# Creates conda environments, pulls Docker images, and downloads databases.
# Usage: ./install.sh -b /path/to/Databases [options]
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Defaults
# ------------------------------------------------------------------------------
DB_DIR=""
SKIP_ENVS=false
SKIP_DOCKER=false
SKIP_DBS=false
KRAKEN_DATE="20240904"   # Update to the latest date from genome-idx.s3.amazonaws.com

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
log()   { echo -e "\n\033[1;32m[$(date +'%H:%M:%S')] $1\033[0m"; }
warn()  { echo -e "\033[1;33m[WARNING] $1\033[0m"; }
error() { echo -e "\033[1;31m[ERROR] $1\033[0m"; exit 1; }

show_help() {
    cat << EOF

Aluminion install script

Usage: ./install.sh -b <db_dir> [options]

Required:
  -b, --db-dir    Path where databases will be downloaded (e.g. /mnt/data/Databases)

Optional:
  --skip-envs     Skip conda environment creation
  --skip-docker   Skip Docker image pulls
  --skip-dbs      Skip database downloads
  -h, --help      Show this message

EOF
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--db-dir)    DB_DIR="$2"; shift ;;
        --skip-envs)    SKIP_ENVS=true ;;
        --skip-docker)  SKIP_DOCKER=true ;;
        --skip-dbs)     SKIP_DBS=true ;;
        -h|--help)      show_help; exit 0 ;;
        *) error "Unknown parameter: $1" ;;
    esac
    shift
done

if [ -z "$DB_DIR" ]; then
    error "The -b / --db-dir parameter is required."
fi

# Resolve the script directory so relative paths to envs/ always work
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# ==============================================================================
# 1. CONDA ENVIRONMENTS
# ==============================================================================
if [ "$SKIP_ENVS" = false ]; then
    log "Creating conda environments..."

    # Detect mamba or fall back to conda
    if command -v mamba >/dev/null 2>&1; then
        PKG_MANAGER="mamba"
    elif command -v conda >/dev/null 2>&1; then
        PKG_MANAGER="conda"
        warn "mamba not found — using conda (slower). Consider installing Mambaforge."
    else
        error "Neither mamba nor conda found. Install Mambaforge first: https://github.com/conda-forge/miniforge/releases"
    fi

    for env_file in aluminion_reads aluminion_assembly aluminion_circlator aluminion_annot aluminion_integron aluminion_copla aluminion_kleborate; do
        yml="${SCRIPT_DIR}/envs/${env_file}.yml"
        if [ ! -f "$yml" ]; then
            warn "Environment file not found, skipping: $yml"
            continue
        fi
        env_name=$(grep '^name:' "$yml" | awk '{print $2}')
        if conda env list | grep -q "^${env_name} "; then
            warn "Conda environment '${env_name}' already exists. Skipping (remove it manually to reinstall)."
        else
            log "Creating environment: ${env_name}..."
            $PKG_MANAGER env create -f "$yml"
        fi
    done

    log "All conda environments created."
fi

# ==============================================================================
# 2. DOCKER IMAGES
# ==============================================================================
if [ "$SKIP_DOCKER" = false ]; then
    log "Pulling Docker images..."

    if ! command -v docker >/dev/null 2>&1; then
        warn "Docker not found. Skipping Docker image pulls. Install Docker to use MOB-suite, Copla, and Phastest."
    else
        docker pull kbessonov/mob_suite:3.0.3
        docker pull rpalcab/copla:1.0

        log "Docker images pulled: mob_suite, copla."
        echo ""
        warn "Phastest has no public Docker image. Set it up manually following the instructions"
        warn "at https://phastest.ca and place it at \$ALUMINION_PHASTEST_DIR (default: ~/Programs/phastest-docker)."
        warn "Then set: export ALUMINION_PHASTEST_DIR=/path/to/phastest-docker"
    fi
fi

# ==============================================================================
# 3. DATABASES
# ==============================================================================
if [ "$SKIP_DBS" = false ]; then
    mkdir -p "$DB_DIR"

    # --------------------------------------------------------------------------
    # 3a. Kraken2 standard database (~100 GB)
    # --------------------------------------------------------------------------
    KRAKEN_DIR="${DB_DIR}/Kraken"
    if [ -f "${KRAKEN_DIR}/hash.k2d" ]; then
        warn "Kraken2 database already found at ${KRAKEN_DIR}. Skipping download."
    else
        log "Downloading Kraken2 standard database (~100 GB) — this will take a while..."
        mkdir -p "$KRAKEN_DIR"
        KRAKEN_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_${KRAKEN_DATE}.tar.gz"
        wget --progress=dot:giga -O "${KRAKEN_DIR}/k2_standard.tar.gz" "$KRAKEN_URL" \
            || error "Failed to download Kraken2 database. Check the URL or try a different date: https://genome-idx.s3.amazonaws.com/kraken/"
        tar -xzf "${KRAKEN_DIR}/k2_standard.tar.gz" -C "$KRAKEN_DIR"
        rm "${KRAKEN_DIR}/k2_standard.tar.gz"
        log "Kraken2 database ready at ${KRAKEN_DIR}."

        # RAM usage note
        echo ""
        echo "  ┌──────────────────────────────────────────────────────────────────────┐"
        echo "  │  KRAKEN2 RAM NOTE                                                    │"
        echo "  │  The standard database (~100 GB) is copied to /dev/shm (RAM disk)   │"
        echo "  │  before each run. This requires ≥128 GB of available RAM and         │"
        echo "  │  reduces classification time ~10×.                                   │"
        echo "  │                                                                      │"
        echo "  │  If your system has <128 GB RAM, edit aluminion.sh:                 │"
        echo "  │    - Comment out line: cp \${KRAKEN_DB}/*.k2d /dev/shm/             │"
        echo "  │    - Comment out line: rm -f /dev/shm/*.k2d                         │"
        echo "  │    - Change --db /dev/shm to --db \${KRAKEN_DB}                     │"
        echo "  └──────────────────────────────────────────────────────────────────────┘"
        echo ""
    fi

    # --------------------------------------------------------------------------
    # 3b. GAMBIT database
    # --------------------------------------------------------------------------
    GAMBIT_DIR="${DB_DIR}/gambit"
    if ls "${GAMBIT_DIR}"/*.h5 2>/dev/null | head -1 | grep -q '.'; then
        warn "GAMBIT database already found at ${GAMBIT_DIR}. Skipping."
    else
        log "Downloading GAMBIT database..."
        mkdir -p "$GAMBIT_DIR"
        # GAMBIT provides signatures (.h5) and metadata (.gdb) files
        # Check https://github.com/jlumpe/gambit/releases for current URLs
        GAMBIT_BASE="https://storage.googleapis.com/gambit-public/gambit-db"
        wget --progress=dot:mega -P "$GAMBIT_DIR" "${GAMBIT_BASE}/gambit-signatures-2024-09-01.h5" || \
            warn "Could not download GAMBIT signatures. Check https://github.com/jlumpe/gambit for the current URL."
        wget --progress=dot:mega -P "$GAMBIT_DIR" "${GAMBIT_BASE}/gambit-metadata-2024-09-01.gdb" || \
            warn "Could not download GAMBIT metadata. Check https://github.com/jlumpe/gambit for the current URL."
        log "GAMBIT database ready at ${GAMBIT_DIR}."
    fi

    # --------------------------------------------------------------------------
    # 3c. Bakta database (~30 GB full, ~2 GB light)
    # --------------------------------------------------------------------------
    BAKTA_DIR="${DB_DIR}/bakta/db"
    if [ -d "${BAKTA_DIR}" ] && [ "$(ls -A ${BAKTA_DIR} 2>/dev/null)" ]; then
        warn "Bakta database already found at ${BAKTA_DIR}. Skipping."
    else
        log "Downloading Bakta full database (~30 GB) — this will take a while..."
        mkdir -p "${DB_DIR}/bakta"
        if command -v bakta_db >/dev/null 2>&1; then
            bakta_db download --output "${DB_DIR}/bakta" --type full
        elif conda run -n aluminion_annot bakta_db --version >/dev/null 2>&1; then
            conda run -n aluminion_annot bakta_db download --output "${DB_DIR}/bakta" --type full
        else
            warn "bakta_db not found. Activate aluminion_annot and run:"
            warn "  bakta_db download --output ${DB_DIR}/bakta --type full"
        fi
        log "Bakta database ready at ${BAKTA_DIR}."
    fi

    # --------------------------------------------------------------------------
    # 3d. ISfinder — insertion sequence nucleotide FASTA
    # --------------------------------------------------------------------------
    ISFINDER_DIR="${DB_DIR}/ISfinder"
    if [ -f "${ISFINDER_DIR}/ISfinder-nucl.fasta" ]; then
        warn "ISfinder database already found. Skipping."
    else
        log "Downloading ISfinder nucleotide database..."
        mkdir -p "$ISFINDER_DIR"
        wget --progress=dot:mega \
            "https://www.is-finder.org/download/IS.fna" \
            -O "${ISFINDER_DIR}/ISfinder-nucl.fasta" \
            || warn "Could not download ISfinder FASTA. Download manually from https://www.is-finder.org/download.html and save to ${ISFINDER_DIR}/ISfinder-nucl.fasta"
        log "ISfinder database ready at ${ISFINDER_DIR}."
    fi

    # --------------------------------------------------------------------------
    # 3e. Abricate databases (automatic, requires aluminion_annot)
    # --------------------------------------------------------------------------
    log "Updating Abricate databases..."
    if conda run -n aluminion_annot abricate --version >/dev/null 2>&1; then
        for db in ncbi resfinder card argannot vfdb; do
            conda run -n aluminion_annot abricate-get_db --db "$db" 2>/dev/null || \
                warn "Could not download Abricate database: $db"
        done
        log "Abricate databases updated."
    else
        warn "aluminion_annot environment not found. Activate it and run: abricate-get_db --db ncbi (and resfinder, card, argannot, vfdb)"
    fi

fi

# ==============================================================================
# SUMMARY
# ==============================================================================
log "Installation complete."
echo ""
echo "  Next steps:"
echo "  1. Install Dorado manually: https://github.com/nanoporetech/dorado/releases"
echo "     Place the binary in your PATH (e.g. /usr/local/bin/dorado)"
echo ""
echo "  2. Set up Phastest docker-compose: https://phastest.ca"
echo "     Then set: export ALUMINION_PHASTEST_DIR=/path/to/phastest-docker"
echo ""
echo "  3. Set the MinKNOW data path if different from the default:"
echo "     export ALUMINION_MINKNOW_DIR=/var/lib/minknow/data"
echo ""
echo "  4. Run the pipeline:"
echo "     ./aluminion.sh -r RUN_NAME -b ${DB_DIR} -t 30 -l /path/to/list_seq.tsv"
echo ""
