#!/bin/bash

# ==============================================================================
# Aluminion — Batch wrapper for sequential run processing
# ==============================================================================
#
# Loops over a list of run names and invokes `aluminion` on each. Each run lives
# in its own subfolder <parent>/<run>/ which MUST contain its own list_seq.tsv
# (recurrent multi-folder analyses cannot share a single parent-level sheet, or
# sample IDs / barcodes would collide). The wrapper enforces this by passing
# --require-run-list to every invocation.
#
# By default the wrapper behaves exactly like a direct `aluminion` call: it lets
# aluminion import fastq_pass/ + reports from the MinKNOW data tree
# (/var/lib/minknow/data by default). Pass --skip-import-from-minknow if the run
# folders are already self-contained (fastq_pass/ copied in beforehand); only then
# is a populated fastq_pass/ required up front.
#
# Per-run logic:
#   - If <parent>/<run>/ is missing or has no list_seq.tsv, the run is skipped.
#   - If <parent>/<run>/Aluminion_Report.html exists, the run is skipped
#     (unless --force is passed).
#   - With --skip-import-from-minknow, a populated fastq_pass/ is also required.
#
# The MGE repository (for the cross-run alert system) is created automatically on
# the first run; pass --init-repo to force/declare it explicitly (forwarded to the
# first launched run only). Extra flags after `--` are forwarded verbatim, e.g.
#     aluminion_batch.sh --runlist runs.tsv -b /db -d /seqs/project -- --resume
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "$(readlink -f "${BASH_SOURCE[0]}")" )" &> /dev/null && pwd )"
ALUMINION="${ALUMINION_BIN:-${SCRIPT_DIR}/aluminion.sh}"

RUNLIST=""
PARENT_DIR=""
FORCE=""
SKIP_IMPORT=""        # when set, require a pre-populated fastq_pass/ and forward the skip flag
INIT_REPO_FLAG=""     # when set, "--init-repo" is forwarded to the first launched run
PASSTHROUGH=()

show_help() {
    # Quoted 'EOF' so the backtick-quoted examples below are printed literally;
    # an unquoted heredoc would command-substitute them (e.g. run `--` as a command).
    cat << 'EOF'

Aluminion batch wrapper: process multiple MinKNOW runs sequentially.

Usage: aluminion_batch.sh --runlist <file> -d <parent_dir> [options] [-- <aluminion options>]

Mandatory options:
  --runlist <file>     Plain-text file listing run names, one per line.
                       Lines starting with '#' and blank lines are ignored.
                       Runs are processed in the order they appear.
  -d, --dir <path>     Parent directory containing the run subfolders.

Additional options:
  -b, --db-dir <path>  Path to the databases root (forwarded to aluminion).
  -t, --threads <N>    Thread count (forwarded to aluminion).
  --skip-import-from-minknow
                       Do NOT import reads from the MinKNOW data tree; assume each
                       run folder already contains a populated fastq_pass/. Without
                       this flag the wrapper imports from MinKNOW just like a direct
                       aluminion call (the default).
  --init-repo          Force creation of the MGE repository (forwarded to the first
                       launched run). The repository is created automatically anyway
                       on the first run; this just declares it explicitly.
  --force              Re-run even if <run>/Aluminion_Report.html exists.
  -h, --help           Show this message.

Each run folder MUST contain its own list_seq.tsv (enforced via --require-run-list).

Forwarding aluminion options (the standalone "--" is REAL syntax, not a typo):
  This wrapper only understands the options listed above. Every flag meant for
  aluminion itself (e.g. --resume, --skip-kraken, --polish-batchsize) must be placed
  AFTER a bare "--" separator. The "--" marks "forward everything past here verbatim
  to each aluminion invocation". Passing such a flag without the leading "--" is
  rejected as an unknown argument.
    Correct:   aluminion_batch.sh --runlist runs.tsv -d /seqs -t 30 -- --resume
    Wrong:     aluminion_batch.sh --runlist runs.tsv -d /seqs -t 30 --resume

Example:
  aluminion_batch.sh --runlist runs.tsv -d /seqs/KLEBIRE -b /db -t 30 --init-repo \
      -- --resume --skip-kraken

EOF
}

log()       { echo -e "\n\033[1;32m[$(date +'%Y-%m-%d %H:%M:%S')] $1\033[0m"; }
error_log() { echo -e "\n\033[1;31m[ERROR] $1\033[0m"; }
warn()      { echo -e "\033[1;33m[WARNING] $1\033[0m"; }

# Manual arg parsing: stop at `--` and forward the rest to aluminion.
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --runlist) RUNLIST="$2"; shift 2 ;;
        -d|--dir) PARENT_DIR="$2"; shift 2 ;;
        -b|--db-dir) PASSTHROUGH+=("-b" "$2"); shift 2 ;;
        -t|--threads) PASSTHROUGH+=("-t" "$2"); shift 2 ;;
        --skip-import-from-minknow|--no-minknow) SKIP_IMPORT=true; shift ;;
        --init-repo) INIT_REPO_FLAG="--init-repo"; shift ;;
        --force) FORCE=true; shift ;;
        -h|--help) show_help; exit 0 ;;
        --) shift; PASSTHROUGH+=("$@"); break ;;
        *) error_log "Unknown argument: $1"; show_help; exit 1 ;;
    esac
done

# --- Argument validation -----------------------------------------------------
[ -n "$RUNLIST" ]   || { error_log "--runlist is mandatory."; show_help; exit 1; }
[ -n "$PARENT_DIR" ] || { error_log "-d / --dir is mandatory."; show_help; exit 1; }
[ -f "$RUNLIST" ]   || { error_log "Run list not found: $RUNLIST"; exit 1; }
[ -d "$PARENT_DIR" ] || { error_log "Parent directory not found: $PARENT_DIR"; exit 1; }
[ -x "$ALUMINION" ] || { error_log "Aluminion script not executable: $ALUMINION"; exit 1; }

RUNLIST="$(readlink -f "$RUNLIST")"
PARENT_DIR="$(readlink -f "$PARENT_DIR")"

log "Batch wrapper starting"
log "  Run list   : $RUNLIST"
log "  Parent dir : $PARENT_DIR"
log "  Forwarded  : ${PASSTHROUGH[*]:-(none)}"

# --- Counters for the end-of-batch summary -----------------------------------
processed=()
skipped_done=()
skipped_no_inputs=()
failed=()

# Read the run list, stripping comments / blank lines.
while IFS= read -r run_name || [ -n "$run_name" ]; do
    run_name="${run_name%%#*}"            # drop trailing comments
    run_name="$(echo "$run_name" | xargs)" # trim whitespace
    [ -z "$run_name" ] && continue

    run_dir="${PARENT_DIR}/${run_name}"

    log "─── ${run_name} ────────────────────────────────────────────────"

    if [ ! -d "$run_dir" ]; then
        warn "Run folder not found, skipping: $run_dir"
        failed+=("$run_name (folder missing)")
        continue
    fi

    # Each run must carry its own list_seq.tsv (no shared parent-level sheet, or
    # sample IDs / barcodes would collide across runs).
    if [ ! -f "${run_dir}/list_seq.tsv" ]; then
        warn "  No list_seq.tsv inside ${run_dir} — skipping."
        warn "  Each run folder must contain its own list_seq.tsv for recurrent analyses."
        skipped_no_inputs+=("$run_name")
        continue
    fi

    if [ -z "$FORCE" ] && [ -f "${run_dir}/Aluminion_Report.html" ]; then
        log "  Already analysed (Aluminion_Report.html present). Skipping. Use --force to re-run."
        skipped_done+=("$run_name")
        continue
    fi

    # Per-run aluminion flags. --require-run-list forbids the parent list_seq.tsv fallback.
    run_flags=(--require-run-list)

    # By default aluminion imports fastq_pass/ from MinKNOW. Only when the user asked to
    # skip that import do we require a pre-populated fastq_pass/ and forward the flag.
    if [ -n "$SKIP_IMPORT" ]; then
        if [ ! -d "${run_dir}/fastq_pass" ] || [ -z "$(ls -A "${run_dir}/fastq_pass" 2>/dev/null)" ]; then
            warn "  --skip-import-from-minknow set but no populated fastq_pass/ in ${run_dir} — skipping."
            skipped_no_inputs+=("$run_name")
            continue
        fi
        run_flags+=(--skip-import-from-minknow)
    fi

    # Forward --init-repo to the first launched run only (Repository.init is idempotent,
    # and aluminion auto-creates the repo on the first run regardless).
    if [ -n "$INIT_REPO_FLAG" ]; then
        run_flags+=("$INIT_REPO_FLAG")
        INIT_REPO_FLAG=""
    fi

    log "  Launching aluminion for ${run_name}..."
    # Continue the batch even if a single run fails — record it and move on.
    if "$ALUMINION" -r "$run_name" -d "$PARENT_DIR" "${run_flags[@]}" "${PASSTHROUGH[@]}"; then
        processed+=("$run_name")
    else
        rc=$?
        error_log "  Aluminion exited with code ${rc} for ${run_name}. Continuing with next run."
        failed+=("$run_name (exit ${rc})")
    fi
done < "$RUNLIST"

# --- Final summary -----------------------------------------------------------
echo ""
log "════════ Batch summary ════════"
log "  Processed   : ${#processed[@]}"
[ ${#processed[@]}        -gt 0 ] && log "     ${processed[*]}"
log "  Already done: ${#skipped_done[@]}"
[ ${#skipped_done[@]}     -gt 0 ] && log "     ${skipped_done[*]}"
log "  Skipped (no inputs): ${#skipped_no_inputs[@]}"
[ ${#skipped_no_inputs[@]} -gt 0 ] && warn "     ${skipped_no_inputs[*]}"
log "  Failed      : ${#failed[@]}"
[ ${#failed[@]}           -gt 0 ] && error_log "     ${failed[*]}"

# Non-zero exit if any run failed, so CI / wrapper scripts can detect it.
[ ${#failed[@]} -eq 0 ]
