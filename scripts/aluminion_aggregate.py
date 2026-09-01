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

DEFAULT: ONE ROW PER STRAIN, THE DEFINITIVE RESULT
--------------------------------------------------
Sample IDs are *not* unique across runs — a resequenced isolate deliberately
keeps its ID (that is what `Barcode_rep1/2` in the cumulative DB is for). In the
Pantoea project, ID 143 appears in Ia, II and III.

Those three appearances are not three independent results. `aluminion.sh`
maintains a per-sample read accumulator in the shared repository
(`$REPO/01_reads/<id>.fastq.gz`): for a sample flagged `is_repeated` in
`list_seq.tsv`, the prior reads are *prepended* to the new run's barcode FASTQs
before assembly, and the accumulator is refreshed
(`aluminion.sh`, read-preparation loop). So each successive run assembles a
strictly larger read set, and the appearance from the LAST run is the definitive
one — the assembly of all reads collected for that strain. The earlier
appearances are superseded intermediates.

Hence the default: one row per strain, taken from the last run in which it
appears, with the bare ID as the key and no run suffix. `Run` is kept as a
column so it is visible which run produced the definitive result, and
`list_seq.tsv` additionally carries `N_runs` / `Runs_all`.

"Last" is the last occurrence in the run list (`runs.tsv`), not the latest
calendar date. That is deliberate: the accumulator is filled in batch
*processing* order, so the run that assembled the full read set is the one
processed last, whatever its sequencing date.

The pipeline only merges reads when the sample is flagged `is_repeated` in that
run's `list_seq.tsv`. If the definitive appearance is *not* flagged while
earlier appearances exist, its assembly is NOT cumulative and this collapse
silently discards real data — so that case is reported as a warning instead of
passing quietly.

OPT-OUT: --keep-run-suffix
--------------------------
`--keep-run-suffix` restores the per-appearance view: every (sample, run) pair
is kept and the key is rewritten to ``<ID>_<run>``, the same convention
`mge_repository.py` uses for its host/plasmid UIDs
(``143_III__pl__AC137.fasta``). The rewrite is what makes that mode safe — a
naive concat would collide the repeated IDs and every merge inside
`aluminion_reporter.py` (nine tables joined on the sample key) would fan out
into a cartesian product. Use it to compare runs, or to inspect how a strain's
assembly improved as reads accumulated.

Either way `aluminion_reporter.py` stays completely unmodified: it still sees
exactly one row per key in each of its nine input tables.

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
    ("AbR_report.csv",       "\t",  ["Sample"],                     "suffix"),
    ("VF_report.csv",        "\t",  ["Sample"],                     "suffix"),
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


# Raw ABRicate summaries key their rows by the input PATH
# ('08_Anotacion/143/abricate/143.tab'), not by sample id — unlike their parsed
# counterparts AbR_modif.xlsx / VF_modif.xlsx, whose '#FILE' is the bare id.
# Without deriving a real key here these two tables silently escape both the
# collapse and the suffixing, so a resequenced strain keeps one row per run.
PATH_KEYED = {"AbR_report.csv": "#FILE", "VF_report.csv": "#FILE"}


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


def build_definitive_map(runs_dir: Path, runs):
    """Decide, for each strain ID, which run holds its definitive result.

    Returns (chosen, appearances, flags) where
      chosen      = {id: run}          the last run in `runs` order holding it
      appearances = {id: [run, ...]}   every run holding it, in `runs` order
      flags       = {(id, run): bool}  is_repeated set for that appearance

    See the module docstring: the last appearance is definitive because
    aluminion.sh prepends the accumulated prior reads before assembling an
    is_repeated sample, so each successive run assembles a larger read set.
    """
    chosen, appearances, flags = {}, {}, {}
    for run in runs:
        path = runs_dir / run / "list_seq.tsv"
        if not path.is_file():
            log.warning("%s/list_seq.tsv absent — its samples cannot be "
                        "collapsed and will be missing from the aggregate.", run)
            continue
        df = _read(path, "\t")
        key = next((k for k in ("ID", "ID único", "ID unico")
                    if k in df.columns), None)
        if key is None:
            log.warning("%s/list_seq.tsv has no ID column (found: %s) — skipped.",
                        run, ", ".join(map(str, df.columns[:6])))
            continue
        rep_col = next((c for c in ("is_repeated", "repetida", "Repetida")
                        if c in df.columns), None)
        for _, row in df.iterrows():
            sid = row[key]
            if pd.isna(sid) or str(sid).strip() == "":
                continue
            sid = str(sid).strip()
            chosen[sid] = run
            appearances.setdefault(sid, []).append(run)
            rep = row.get(rep_col) if rep_col else None
            flags[(sid, run)] = bool(rep is not None and not pd.isna(rep)
                                     and str(rep).strip() != "")

    # A strain whose definitive appearance is not flagged is_repeated was
    # assembled from that run's reads ALONE — collapsing to it drops the earlier
    # runs' reads from the picture. Report it rather than hiding it.
    not_merged = [(sid, apps) for sid, apps in appearances.items()
                  if len(apps) > 1 and not flags.get((sid, chosen[sid]), False)]
    if not_merged:
        log.warning("%d resequenced strain(s) whose definitive run is NOT "
                    "flagged is_repeated — their assembly is not cumulative, so "
                    "collapsing hides the earlier run(s). Mark is_repeated in "
                    "the definitive run's list_seq.tsv and re-run, or use "
                    "--keep-run-suffix to see every appearance:", len(not_merged))
        for sid, apps in sorted(not_merged):
            log.warning("    %-12s appears in %s; kept %s (not flagged)",
                        sid, ", ".join(apps), chosen[sid])

    repeats = {s: a for s, a in appearances.items() if len(a) > 1}
    if repeats:
        log.info("Collapsing %d resequenced strain(s) to their definitive run:",
                 len(repeats))
        for sid, apps in sorted(repeats.items()):
            log.info("    %-12s %s  ->  %s", sid, " ".join(apps), chosen[sid])
    return chosen, appearances, flags


def aggregate_table(runs_dir: Path, runs, fname: str, sep: str, keys, mode: str,
                    chosen=None):
    """Concatenate `fname` across `runs`.

    When `chosen` is None the key is rewritten to '<value>_<run>' and every
    appearance is kept (--keep-run-suffix). When `chosen` is a {id: run} map
    only rows belonging to each id's definitive run survive and the key is left
    untouched.

    Returns (df, n_runs_contributing, missing_runs).
    """
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

        path_col = PATH_KEYED.get(fname)
        if path_col and path_col in df.columns:
            # '08_Anotacion/143/abricate/143.tab' -> '143'
            df.insert(0, "Sample",
                      df[path_col].astype(str).map(lambda p: Path(p).stem))

        if mode == "suffix":
            key = next((k for k in keys if k in df.columns), None)
            if key is None:
                log.warning("%s/%s has none of the expected key columns %s "
                            "(found: %s) — concatenated as-is, rows may collide.",
                            run, fname, keys, ", ".join(map(str, df.columns[:6])))
            else:
                # Blank/NaN keys must neither become 'nan_Ia' nor be looked up.
                mask = df[key].notna() & (df[key].astype(str).str.strip() != "")
                df.loc[mask, key] = df.loc[mask, key].astype(str).str.strip()
                if chosen is None:
                    df.loc[mask, key] = df.loc[mask, key] + f"_{run}"
                else:
                    # Keep only rows whose strain has this run as its definitive
                    # one. An unknown key (present in a table but not in any
                    # list_seq.tsv) is kept, so a parser writing an unlisted
                    # sample stays visible instead of vanishing silently.
                    keep = df[key].map(lambda v: chosen.get(v, run) == run)
                    df = df[keep.fillna(True)]

        frames.append(df)

    if not frames:
        return None, 0, missing

    out = pd.concat(frames, ignore_index=True, sort=False)
    # Run first: it is the one column a reader of the aggregate always wants.
    cols = ["Run"] + [c for c in out.columns if c != "Run"]
    return out[cols], len(frames) - sum(1 for f in frames if f.empty), missing


def link_bandage_plots(runs_dir: Path, runs, out_dir: Path, chosen=None) -> int:
    """Mirror 03_assemblies/<sample>.png into the aggregate.

    aluminion_reporter.py embeds these Bandage graphs as base64 hover images
    keyed on the sample name (`inject_hover`), resolved relative to its working
    directory — so the mirrored filename must match the aggregated table key
    exactly, or the aggregate report loses every hover plot.

    With `chosen` (default mode) only the definitive run's plot is mirrored,
    under the bare sample name. Without it, every plot is mirrored as
    <sample>_<run>.png to match the suffixed keys.
    """
    dest = out_dir / "03_assemblies"
    dest.mkdir(parents=True, exist_ok=True)
    # Stale plots from a previous invocation in the other mode would be picked
    # up by inject_hover() and silently attached to the wrong row.
    for old in dest.glob("*.png"):
        old.unlink()
    n = 0
    for run in runs:
        src_dir = runs_dir / run / "03_assemblies"
        if not src_dir.is_dir():
            continue
        for png in src_dir.glob("*.png"):
            if chosen is None:
                target = dest / f"{png.stem}_{run}.png"
            else:
                if chosen.get(png.stem, run) != run:
                    continue
                target = dest / f"{png.stem}.png"
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
    ap.add_argument("--keep-run-suffix", action="store_true",
                    help="Keep every (sample, run) appearance as its own row, "
                         "with the key rewritten to '<ID>_<run>'. Default is "
                         "one row per strain, taken from the last run it "
                         "appears in (the assembly of all accumulated reads).")
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
    mode_label = ("per-appearance ('<ID>_<run>' keys)" if args.keep_run_suffix
                  else "definitive result per strain (last run, bare ID)")
    log.info("Aggregating %d run(s) into %s: %s", len(runs), out_dir,
             ", ".join(runs))
    log.info("Mode: %s", mode_label)

    if args.keep_run_suffix:
        chosen, appearances = None, {}
    else:
        chosen, appearances, _flags = build_definitive_map(runs_dir, runs)

    written, samples = [], set()
    for fname, sep, keys, mode in TABLES:
        df, n_found, missing = aggregate_table(runs_dir, runs, fname, sep,
                                               keys, mode, chosen=chosen)
        if df is None:
            log.warning("%-22s absent from every run — not written.", fname)
            continue

        key = next((k for k in keys if k in df.columns), None)
        if fname == "list_seq.tsv" and key:
            if appearances:
                # Provenance of the collapse, for whoever reads the TSV.
                df["N_runs"] = df[key].map(
                    lambda v: str(len(appearances.get(v, [v]))))
                df["Runs_all"] = df[key].map(
                    lambda v: ",".join(appearances.get(v, [])))
            samples = set(df[key].dropna())

        _write(df, out_dir / fname, sep)
        written.append(fname)
        note = f" (absent in: {', '.join(missing)})" if missing else ""
        cov = ""
        if samples and key and fname != "list_seq.tsv":
            n_cov = len(samples & set(df[key].dropna()))
            cov = f", {n_cov}/{len(samples)} strains"
        log.info("%-22s %4d rows from %d run(s)%s%s",
                 fname, len(df), n_found, cov, note)

    n_png = link_bandage_plots(runs_dir, runs, out_dir, chosen=chosen)
    log.info("Mirrored %d Bandage plot(s) into %s/03_assemblies/", n_png, out_dir)

    if samples:
        log.info("Aggregate covers %d %s across %d run(s).", len(samples),
                 "sample-appearance(s)" if args.keep_run_suffix else "strain(s)",
                 len(runs))

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
