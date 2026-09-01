#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aluminion_aggregate.py — build a batch-level ("all runs") view of an Aluminion
project.

`aluminion.sh` produces one `Aluminion_Report.html` and one set of analysis
tables per *run*. When several runs belong to the same project (a batch driven
by `aluminion_batch.sh --runlist runs.tsv`), there is no single place to look at
every strain at once: the only cross-run artefacts are `data_seq.tsv` /
`data_analysis.tsv` (cumulative lab DB, one row per sample, no per-element
detail) and `repository/` (plasmid/integron FASTAs only).

This script fills that gap. It concatenates every per-run analysis table into an
aggregate directory and then runs `aluminion_reporter.py` on it, yielding a
single `Aluminion_Report.html` covering all strains in the batch.

WHY THE SAMPLE KEY IS SUFFIXED WITH THE RUN NAME
------------------------------------------------
Sample IDs are *not* unique across runs — a resequenced isolate deliberately
keeps its ID (that is what `Barcode_rep1/2` in the cumulative DB is for). In the
Pantoea project, ID 143 appears in Ia, II and III. A naive concat would collide
those into indistinguishable rows and every downstream merge in
`aluminion_reporter.py` (which merges the nine input tables on the sample key)
would fan out into a cartesian product.

So the key is rewritten as ``<ID>_<run>``, the same convention
`mge_repository.py` already uses for its host/plasmid UIDs
(``143_III__pl__AC137.fasta``). This keeps `aluminion_reporter.py` completely
unmodified — it still sees one row per key — and makes the run visible in the
report's first column. A `Run` column is added to every table as well, for
whoever loads the aggregated TSV/CSV directly.

Tables whose key is already run-qualified (`alerts.tsv`) are concatenated
without rewriting. The per-run `data_seq.tsv` / `data_analysis.tsv` snapshots
are deliberately NOT aggregated here: the cumulative versions maintained by
`lab_db_updater.py` at the project root are already the batch-level view of
those two, and rewriting their keys would desynchronise them from the
repository's host UIDs.

Usage
-----
    conda activate aluminion_annot
    python3 scripts/aluminion_aggregate.py -d /path/to/project [-l runs.tsv]

    # tables only, skip the HTML report
    python3 scripts/aluminion_aggregate.py -d /path/to/project --no-report
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _log import get_logger  # noqa: E402

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Table registry
#
# mode='suffix' -> rewrite the sample key to '<value>_<run>' (unique per batch)
# mode='plain'  -> concatenate as-is, only add the Run column
#
# `keys` lists candidate column names in priority order: the same table can
# carry a different header depending on which parser wrote it (and legacy
# Spanish headers still turn up in hand-maintained list_seq.tsv sheets).
# ---------------------------------------------------------------------------
TABLES = [
    # --- the nine inputs aluminion_reporter.py reads (must keep name+format) --
    ("list_seq.tsv",         "\t",  ["ID", "ID único", "ID unico"], "suffix"),
    ("QC_reads.csv",         "\t",  ["Sample"],                     "suffix"),
    ("QC_assembly.csv",      "\t",  ["Samples"],                    "suffix"),
    ("taxonomy.xlsx",        "xlsx", ["Sample"],                    "suffix"),
    ("AbR_modif.xlsx",       "xlsx", ["#FILE"],                     "suffix"),
    ("copla_modif.csv",      ",",   ["Sample"],                     "suffix"),
    ("integron_summary.csv", ",",   ["Sample"],                     "suffix"),
    ("phage_summary.csv",    ",",   ["Sample"],                     "suffix"),
    ("kleborate.tsv",        "\t",  ["strain"],                     "suffix"),
    # --- other per-run analysis tables, aggregated for direct inspection -----
    ("AbR_report.csv",       "\t",  ["#FILE"],                      "suffix"),
    ("VF_report.csv",        "\t",  ["#FILE"],                      "suffix"),
    ("VF_modif.xlsx",        "xlsx", ["#FILE"],                     "suffix"),
    ("taxonomy.csv",         ",",   ["Sample"],                     "suffix"),
    ("kraken.csv",           "\t",  ["Sample"],                     "suffix"),
    ("kraken_mlst.xlsx",     "xlsx", ["Sample"],                    "suffix"),
    ("mlst_modif.csv",       "\t",  ["Sample"],                     "suffix"),
    ("IS.tsv",               "\t",  ["sample"],                     "suffix"),
    ("samplesheet.tsv",      "\t",  ["ID"],                         "suffix"),
    ("alerts.tsv",           "\t",  [],                             "plain"),
]

# Deliberately excluded — see module docstring.
#   data_seq.tsv, data_analysis.tsv : cumulative versions already exist at root
#   mlst.csv                        : raw headerless `mlst` output, superseded
#                                     by mlst_modif.csv


def _read(path: Path, sep: str) -> pd.DataFrame:
    """Read a table, forcing every column to str.

    dtype=str is not cosmetic. pd.read_csv infers the key column's dtype per
    file, so a run whose IDs are all numeric reads back as int64 while a run
    containing one alphanumeric ID reads back as object — concatenating those
    two produces a mixed-dtype column in which '143' and 143 are different
    values. This is the same class of bug that silently duplicated rows and
    zeroed the MGE counts in lab_db_updater.py (commit 6168f57); pinning the
    dtype at every load point is the fix that generalises.
    """
    if sep == "xlsx":
        return pd.read_excel(path, dtype=str)
    return pd.read_csv(path, sep=sep, dtype=str)


def _write(df: pd.DataFrame, path: Path, sep: str) -> None:
    if sep == "xlsx":
        df.to_excel(path, index=False)
    else:
        df.to_csv(path, sep=sep, index=False)


def aggregate_table(runs_dir: Path, runs, fname: str, sep: str, keys, mode: str):
    """Concatenate `fname` across `runs`. Returns (df, n_runs_found, n_missing)."""
    frames, missing = [], []
    for run in runs:
        path = runs_dir / run / fname
        if not path.is_file():
            missing.append(run)
            continue
        try:
            df = _read(path, sep)
        except Exception as exc:                                  # noqa: BLE001
            log.warning("%s/%s unreadable (%s: %s) — skipped.",
                        run, fname, type(exc).__name__, exc)
            missing.append(run)
            continue

        if df.empty:
            # An empty table is legitimate (no integrons found, no alerts
            # raised). Keep its header so the aggregate still has the columns
            # aluminion_reporter.py expects, but contribute no rows.
            frames.append(df.assign(Run=pd.Series(dtype=str)))
            continue

        df = df.copy()
        df["Run"] = run

        if mode == "suffix":
            key = next((k for k in keys if k in df.columns), None)
            if key is None:
                log.warning("%s/%s has none of the expected key columns %s "
                            "(found: %s) — concatenated without run-qualifying "
                            "the key, rows may collide.",
                            run, fname, keys, ", ".join(map(str, df.columns[:6])))
            else:
                # Blank/NaN keys must not become the string 'nan_Ia'.
                mask = df[key].notna() & (df[key].astype(str).str.strip() != "")
                df.loc[mask, key] = df.loc[mask, key].astype(str).str.strip() + f"_{run}"

        frames.append(df)

    if not frames:
        return None, 0, missing

    out = pd.concat(frames, ignore_index=True, sort=False)
    # Run first: it is the one column a reader of the aggregate always wants.
    cols = ["Run"] + [c for c in out.columns if c != "Run"]
    return out[cols], len(frames) - sum(1 for f in frames if f.empty), missing


def link_bandage_plots(runs_dir: Path, runs, out_dir: Path) -> int:
    """Mirror 03_assemblies/<sample>.png as <sample>_<run>.png.

    aluminion_reporter.py embeds these Bandage graphs as base64 hover images
    keyed on the sample name (`inject_hover`), resolved relative to its working
    directory. Without this mirror the aggregate report loses every hover plot.
    """
    dest = out_dir / "03_assemblies"
    dest.mkdir(parents=True, exist_ok=True)
    n = 0
    for run in runs:
        src_dir = runs_dir / run / "03_assemblies"
        if not src_dir.is_dir():
            continue
        for png in src_dir.glob("*.png"):
            target = dest / f"{png.stem}_{run}.png"
            if target.exists() or target.is_symlink():
                target.unlink()
            # Copy rather than symlink: the aggregate directory should stay
            # readable after the run folders are pruned by end-of-run cleanup.
            shutil.copy2(png, target)
            n += 1
    return n


def resolve_runs(runs_dir: Path, runlist: Path | None):
    """Authoritative run list. Never glob the project directory.

    Globbing would pick up sibling clutter that is not a run: `repository/`,
    `backup_pre_recovery_*/`, and stray directories left by mis-invoked runs
    (the Pantoea project has an empty `IV/Ia/` from an aborted 2026-08-27 run).
    """
    if runlist is not None:
        if not runlist.is_file():
            log.error("Run list not found: %s", runlist)
            sys.exit(1)
        runs = [ln.split("\t")[0].strip()
                for ln in runlist.read_text(encoding="utf-8").splitlines()
                if ln.strip() and not ln.startswith("#")]
    else:
        default = runs_dir / "runs.tsv"
        if default.is_file():
            log.info("Using run list: %s", default)
            return resolve_runs(runs_dir, default)
        log.error("No run list given and no runs.tsv in %s. Pass -l/--runlist.",
                  runs_dir)
        sys.exit(1)

    valid, skipped = [], []
    for r in runs:
        (valid if (runs_dir / r).is_dir() else skipped).append(r)
    if skipped:
        log.warning("Listed but not present as directories, skipped: %s",
                    ", ".join(skipped))
    if not valid:
        log.error("None of the listed runs exist under %s.", runs_dir)
        sys.exit(1)
    return valid


def main():
    ap = argparse.ArgumentParser(
        description="Aggregate every per-run Aluminion analysis table into a "
                    "single batch-level directory and build one "
                    "Aluminion_Report.html covering all strains.")
    ap.add_argument("-d", "--runs-dir", required=True, type=Path,
                    help="Project directory containing the run folders.")
    ap.add_argument("-l", "--runlist", type=Path, default=None,
                    help="TSV whose first column lists the run folder names "
                         "(default: <runs-dir>/runs.tsv).")
    ap.add_argument("-o", "--out", type=Path, default=None,
                    help="Output directory (default: <runs-dir>/ALL_RUNS).")
    ap.add_argument("--no-report", action="store_true",
                    help="Write the aggregated tables but skip "
                         "aluminion_reporter.py.")
    args = ap.parse_args()

    runs_dir = args.runs_dir.resolve()
    if not runs_dir.is_dir():
        log.error("Not a directory: %s", runs_dir)
        sys.exit(1)

    runs = resolve_runs(runs_dir, args.runlist)

    out_dir = (args.out or runs_dir / "ALL_RUNS").resolve()
    if out_dir in [(runs_dir / r).resolve() for r in runs]:
        log.error("Output directory %s is one of the run folders. Refusing to "
                  "overwrite a run's own tables.", out_dir)
        sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("Aggregating %d run(s) into %s: %s",
             len(runs), out_dir, ", ".join(runs))

    written, samples = [], set()
    for fname, sep, keys, mode in TABLES:
        df, n_found, missing = aggregate_table(runs_dir, runs, fname, sep,
                                               keys, mode)
        if df is None:
            log.warning("%-22s absent from every run — not written.", fname)
            continue
        _write(df, out_dir / fname, sep)
        written.append(fname)
        note = f" (absent in: {', '.join(missing)})" if missing else ""
        log.info("%-22s %4d rows from %d run(s)%s", fname, len(df), n_found, note)

        if fname == "list_seq.tsv":
            key = next((k for k in keys if k in df.columns), None)
            if key:
                samples = set(df[key].dropna())

    n_png = link_bandage_plots(runs_dir, runs, out_dir)
    log.info("Mirrored %d Bandage plot(s) into %s/03_assemblies/", n_png, out_dir)

    if samples:
        log.info("Aggregate covers %d unique sample(s) across %d run(s).",
                 len(samples), len(runs))

    if args.no_report:
        log.info("--no-report given; skipping HTML generation.")
        return

    reporter = Path(__file__).resolve().parent / "aluminion_reporter.py"
    log.info("Running aluminion_reporter.py on the aggregate...")
    proc = subprocess.run([sys.executable, str(reporter), str(out_dir)],
                          capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        log.error("aluminion_reporter.py failed (exit %d). The aggregated "
                  "tables in %s are still valid — inspect them directly.",
                  proc.returncode, out_dir)
        sys.exit(proc.returncode)

    report = out_dir / "Aluminion_Report.html"
    if report.is_file():
        log.info("Batch-level report: %s", report)
    else:
        log.warning("Reporter exited 0 but %s was not created.", report)


if __name__ == "__main__":
    main()
