#!/bin/bash

# ==============================================================================
# Aluminion — Batch wrapper for sequential run processing
# ==============================================================================
#
# Loops over a list of run names and invokes `aluminion` on each. Designed for
# the workflow where many MinKNOW runs have already been organised under a
# single parent directory, each one self-contained (its own fastq_pass/,
# list_seq.tsv, and MinKNOW reports already in place).
#
# Per-run logic:
#   - If <parent>/<run>/Aluminion_Report.html exists, the run is skipped
#     (unless --force is passed).
#   - If <parent>/<run>/fastq_pass/ exists, aluminion is invoked with
#     --no-minknow (no copy from /var/lib/minknow/data).
#   - Otherwise the run is skipped with a warning.
#
# Extra flags after `--` are forwarded verbatim to aluminion. Example:
#     aluminion_batch.sh --runlist runs.tsv -b /db -d /seqs/project -- --resume
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "$(readlink -f "${BASH_SOURCE[0]}")" )" &> /dev/null && pwd )"
ALUMINION="${ALUMINION_BIN:-${SCRIPT_DIR}/aluminion.sh}"

RUNLIST=""
PARENT_DIR=""
FORCE=""
PASSTHROUGH=()

show_help() {
    cat << EOF

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
  --force              Re-run even if <run>/Aluminion_Report.html exists.
  -h, --help           Show this message.

Any flag after `--` is forwarded verbatim to each aluminion invocation, e.g.
`-- --resume --skip-kraken --polish-batchsize 8`.

Example:
  aluminion_batch.sh --runlist runs.tsv -d /seqs/KLEBIRE -b /db -t 30 \\
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

    if [ -z "$FORCE" ] && [ -f "${run_dir}/Aluminion_Report.html" ]; then
        log "  Already analysed (Aluminion_Report.html present). Skipping. Use --force to re-run."
        skipped_done+=("$run_name")
        continue
    fi

    if [ ! -d "${run_dir}/fastq_pass" ] || [ -z "$(ls -A "${run_dir}/fastq_pass" 2>/dev/null)" ]; then
        warn "  No populated fastq_pass/ in ${run_dir} — skipping."
        warn "  Either copy MinKNOW outputs into the run folder or run aluminion directly without --no-minknow."
        skipped_no_inputs+=("$run_name")
        continue
    fi

    log "  Launching aluminion for ${run_name}..."
    # Continue the batch even if a single run fails — record it and move on.
    if "$ALUMINION" -r "$run_name" -d "$PARENT_DIR" --no-minknow "${PASSTHROUGH[@]}"; then
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
