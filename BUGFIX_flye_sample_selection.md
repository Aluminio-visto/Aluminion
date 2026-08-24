# Bug: a sample deleted from `list_seq.tsv` was still handed to Flye

Date: 2026-08-24 · Reported on run `ENTHERE_2026_JUL_07`, sample `SCT-HURS-74`
(Enthere batch, `aluminion_batch --runlist runs.tsv -d ./ -b ~/Databases -t 30 -- --repo ./repository --resume`)

## Symptom

`SCT-HURS-74` sequenced badly (likely contaminated / unassemblable). Flye spun
forever inside `Extending reads` (`0% 80%` and no further progress), stalling the
whole batch. The user removed the sample's row from `list_seq.tsv` and deleted
`03_assemblies/SCT-HURS-74/`, yet every `--resume` re-queued it for assembly.

## Root causes — three independent defects

### 1. `samples` was derived from disk state, not from the sample sheet

`aluminion.sh` built its per-run sample list with:

```bash
find 01_reads -type f -name "*.fastq.gz" -size "+${MIN_READ_MB}M" \
    -exec basename {} .fastq.gz \; | sort | uniq > samples
```

`list_seq.tsv` drove only the *read-concatenation* loop just above it, which
populates `01_reads/`. The list itself came from whatever FASTQs were on disk.
So the sheet was authoritative on a sample's **first** pass only: once
`01_reads/<id>.fastq.gz` existed, deleting the row changed nothing — the `find`
picked the leftover FASTQ back up on every subsequent `--resume`. Deleting
`03_assemblies/<id>/` did not help either, because that only cleared the
`resume_done` sentinel and made the sample look *unassembled*, i.e. due for a
fresh Flye attempt.

**Fix:** build `samples` as the intersection of the IDs declared in
`list_seq.tsv` (column 3, whitespace- and CRLF-normalised) and the on-disk FASTQs
passing `--min-read-mb`. The sheet is now authoritative on every pass. Reads
belonging to an ID no longer in the sheet are reported once in the log rather
than silently ignored, since ignoring a 400 MB FASTQ without comment looks like
a bug in its own right.

### 2. Flye hanging was never handled — only Flye *failing* was

The assembly-failure handler (skip / retry `--meta` / abort) is reached only when
`flye` **exits**. A sample that makes Flye loop forever never reaches it, so an
unattended batch stalls indefinitely instead of skipping the sample. This is the
reason the run had to be killed by hand.

**Fix:** every Flye invocation now runs through `run_flye()`, which wraps it in
`timeout --foreground --kill-after=60s "$FLYE_TIMEOUT"` (default `4h`, override
with `--flye-timeout` or `$ALUMINION_FLYE_TIMEOUT`, `0` to disable). Exit 124/137
is reported distinctly ("exceeded the limit", not "Flye could not assemble") and
then routed through the existing failure path, which in non-interactive runs
skips the sample and continues. `--foreground` is required so Flye's progress
output still reaches the tee'd log; the follow-up SIGKILL covers a Flye that
ignores SIGTERM.

### 3. Assembly failures were not remembered across resumes

A sample that failed assembly was dropped from `samples` for the remainder of
that run only. The next `--resume` re-queued it and paid the full assembly cost
again — now up to `FLYE_TIMEOUT` per pass — to fail identically.

**Fix:** failures are appended to `<run>/.failed_assemblies` and skipped on
resumed runs. `--retry-failed-assembly` clears the record, for use after
re-sequencing the isolate or raising the timeout.

## Latent defect found while testing (independent of the report)

`find ... > samples` wrote its result directly into `samples`, so an empty scan
**truncated the file**. That is the normal state of a completed run: the
end-of-run cleanup prunes `01_reads/`. Any `--resume` over such a run — exactly
what a `runs.tsv` loop does — truncated `samples` to zero rows, aborted with
"No samples passed the filter", and destroyed the record of which samples the run
contained. Verified against the real data: `2026_06_23` (28 samples) and
`2025_12_15` (12 samples) would both have been wiped.

**Fix:** stage into `.samples_new`, validate, and only then `mv` over `samples`.
When the scan is empty but a non-empty `samples` exists, the run aborts with a
message naming the likely cause (pruned reads) and pointing at
`--just-preprocessing` / `--skip-preprocessing`, leaving `samples` untouched.

## Validation

Logic was exercised in disposable sandboxes with a stub `flye`:

- sample deleted from the sheet but FASTQ still on disk → excluded (the reported bug)
- ID in the sheet with an undersized FASTQ → excluded (size gate still works)
- IDs with padded whitespace / CRLF and blank sheet rows → parsed correctly
- FASTQ with no sheet row → excluded and logged
- stub Flye that hangs → terminated at the timeout, `rc=124`, treated as failure
- stub Flye exiting 0 → unaffected; `--flye-timeout 0` restores legacy behaviour
- failure record → idempotent, honoured on resume, cleared by `--retry-failed-assembly`
- pruned `01_reads/` with a historical `samples` → `samples` preserved, clean abort

Checked: `bash -n` clean; `--help` renders (backticks escaped inside the
variable-expanding heredoc, per the project's heredoc convention); unknown-flag
rejection intact; `RETRY_FAILED_ASSEMBLY` initialised so `set -u` does not abort;
mode `100755` and LF endings preserved. All 13 `scripts/*.py` still import in the
production `aluminion_annot` env.

## Server state repaired for `ENTHERE_2026_JUL_07`

- `samples` rebuilt via the new intersection: 27 → 26 rows (`SCT-HURS-74` gone).
  Previous file kept as `samples.bak.20260824`.
- `SCT-HURS-74` written to `.failed_assemblies` so it stays skipped even if the
  row returns to the sheet.
- The incomplete `03_assemblies/SCT-HURS-74/` staging dir was removed.
- `SCT-HURS-70` is in the sheet but its FASTQ is 121 MB (< `--min-read-mb 135`),
  so it is correctly excluded — pre-existing behaviour, not a regression.
- Audited every run in `runs.tsv`: no other run had a sample in `samples` that
  was absent from its sheet.

## Follow-up for the user — NOT addressed here

**`aluminion_annot` now has pandas 3.0.3.** The project notes assume production
runs pandas 2.x and dev runs 1.5.3; the env has since moved to 3.0.3, a major
release with removals beyond the 2.x changes that produced the earlier
`LossySetitemError` fix. All 13 scripts still import, but importing is not
exercising: `pytest` is **not installed** in `aluminion_annot`, so the suite has
never run against the version production actually uses. Recommend
`mamba install -n aluminion_annot pytest` and a full run before the next batch.

---

# Follow-up (same day): pandas 3 silently disabled the merge-key strips in parser.py

Installing pytest into `aluminion_annot` — per the follow-up above — let the suite
run against the version production actually uses for the first time, and it
immediately found a live bug. **Not a test-harness artefact: `parser.py` produces
wrong output on the next run.**

## Symptom

Under pandas 3.0.3, `taxonomy.csv` came out with `Subspecies`, `MLST`,
`Serotype`, `KO_locus`, `Carbapenemase`, `ESBL` and the allele columns **entirely
empty** (19/19 rows), and `Aluminion_Report.html` lost the corresponding fields.
`parser.py` still exited 0 and logged "Process completed successfully" — the
failure is silent.

## Root cause

Two defensive whitespace-strip loops were gated on the column dtype:

```python
for c in kraken_df.columns:
    if kraken_df[c].dtype == 'object':      # <-- False on pandas 3
        kraken_df[c] = kraken_df[c].astype(str).str.strip()...
```

pandas 3 infers text columns as the dedicated `str` dtype rather than `object`,
so `dtype == 'object'` is False and **every strip in both loops was skipped**.
The `Sample` key coming out of the Kraken report (`Kraken report → awk →
tab-split`) carries a trailing space — `'Eclo_VC_600-1 '` — which the strip
existed precisely to remove. With the strip disabled, all five successive
`pd.merge(..., on='Sample')` joins matched **zero** rows:

```
merge kraken→gambit    claves comunes: 0    (19x19 rows)
merge      →mlst       claves comunes: 0
merge      →kleborate  claves comunes: 0
merge      →ectyper    claves comunes: 0
```

Because the joins are `how='left'`, the row count stayed at 19 and no error was
raised — the merged columns just came back all-NaN. This is the same class of bug
as the pandas-1.5.3 `DataFrame.apply(lambda x: x.str.strip())` no-op the comments
in that function already warn about: the *idiom* was fixed, but the *dtype guard
wrapping it* reintroduced the identical failure on a newer pandas.

## Fix

Drop the dtype guard in both loops (`scripts/parser.py`, the `gambit_df` and
`kraken_df` strips) and strip unconditionally. Both loops run on frames already
narrowed to their final text columns — `gambit_df` to `['Sample','Subspecies']`,
`kraken_df` to `kraken_cols` — so `astype(str)` cannot clobber a numeric column
that matters. Verified: all five merges recover 19/19 keys, no empty cells, and
real values land (`Enterobacter hormaechei subsp. xiangfangensis`, MLST `114`,
`KL53/O3/O3a`). Behaviour is identical under pandas 2.2.2.

Grepped the rest of `scripts/` for the same pattern: these two were the only
occurrences. `lab_db_updater.py` strips unconditionally and was never affected.

## Blast radius

**Existing tables are clean.** Every `taxonomy.csv` under
`/home/usuario/Seqs/Enthere` (runs `2025_08_27`, `2025_09_09`, `2025_09_16`,
`2025_12_15`, `2026_06_23`) shows 0 empty values in `Subspecies` / `MLST` /
`KO_locus` — they were generated when the env still had pandas 2 and the guard
still worked. The bug would have corrupted the *next* run, not past ones. No
historical repair needed.

## Note on running the tests

`pytest` must be run with `PYTHONSAFEPATH` unset. That variable disables Python's
automatic insertion of the script's own directory into `sys.path`, which is how
`scripts/*.py` resolve `from _log import get_logger`. It is not set in a normal
login shell — only in some sandboxed/tooling environments — but if it leaks in,
10 tests fail with a misleading `ModuleNotFoundError: No module named '_log'`
that has nothing to do with the code under test.
