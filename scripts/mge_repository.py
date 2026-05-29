#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Script:  mge_repository.py
# Purpose: Read/write API for the cross-run cumulative mobile genetic element
#          (MGE) repository. The repository accumulates assembled plasmids
#          and integrons across longitudinal Nanopore runs so that recurrent
#          elements can be detected by mge_alerts.py.
# Inputs:  - Plasmid FASTAs from MOB-suite (08_Anotacion/<sample>/mob_recon/)
#          - Integron coordinates from integron_summary.csv
#          - Per-sample host metadata from data_analysis.tsv
# Outputs: - repository/index_plasmids.tsv
#          - repository/index_integrons.tsv
#          - repository/hosts.tsv
#          - repository/plasmids/<uid>.fasta
#          - repository/integrons/<uid>.fasta
#          - repository/sketches/plasmids.sk (skani sketch index)
# Author:  Aluminion lab
# Date:    2026-05-20
# =============================================================================

"""Cumulative MGE repository — ingestion and lookup API.

This module is the single source of truth for the repository layout. Other
modules (``mge_alerts.py``, future ``seed_repo_from_runs.py``) must go
through this API rather than reading/writing the index files directly.

Repository layout
-----------------
::

    <repo>/
    ├── index_plasmids.tsv
    ├── index_integrons.tsv
    ├── hosts.tsv
    ├── plasmids/<uid>.fasta
    ├── integrons/<uid>.fasta
    └── sketches/plasmids.sk

UID conventions
---------------
- ``host_uid``     = ``<isolate_id>_<run_name>``
- ``plasmid_uid``  = ``<host_uid>__pl__<contig>``
- ``integron_uid`` = ``<host_uid>__int__<contig>__<start>-<end>``

Isolate ID is the leading token in ``host_uid`` so that alphabetical sorting
of the index files groups all entries of the same isolate together, which
matches the analyst's natural query pattern. The ``__`` (double underscore)
separator within MGE UIDs is filesystem-safe on every platform (avoiding
characters like ``:`` that are invalid on NTFS) while remaining visually
distinguishable from the single underscore inside isolate IDs and run
names. A separate ``sha1`` field on plasmid rows allows strict
deduplication when the same sequence is ingested from two different runs.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

# -----------------------------------------------------------------------------
# Tunable defaults — overridable by the mge_alerts CLI.
# -----------------------------------------------------------------------------
# Plasmids smaller than this (in bp) are treated as noise and skipped during
# ingestion. The threshold is generous on purpose: small Col-like plasmids
# can legitimately fall just above 1 kb in our datasets.
MIN_PLASMID_SIZE: int = 1000

# -----------------------------------------------------------------------------
# Index schemas — kept as module constants so that any future schema change
# happens in one place. Order matters: it defines column order on disk.
# -----------------------------------------------------------------------------
PLASMID_INDEX_COLUMNS: list[str] = [
    "uid",
    "sha1",
    "host_uid",
    "run_name",
    "sample_id",
    "contig",
    "ptu",
    "rep",
    "mob",
    "mpf",
    "size",
    "amr_genes",
    "vir_genes",
    "ingested_at",
]

INTEGRON_INDEX_COLUMNS: list[str] = [
    "uid",
    "host_uid",
    "run_name",
    "sample_id",
    "contig",
    "integron_type",
    "integrase",
    "gene_set_json",
    "amr_genes",
    "vir_genes",
    "start",
    "end",
    "size",
    "ingested_at",
]

HOSTS_COLUMNS: list[str] = [
    "host_uid",
    "run_name",
    "seq_date",
    "lab_id",
    "isolate_id",
    "strain",
    "genus",
    "species",
    "subspecies",
    "mlst",
    "serotype",
    "ko_locus",
    "amr_score",
    "vir_score",
    "ingested_at",
]


# -----------------------------------------------------------------------------
# UID helpers
# -----------------------------------------------------------------------------
def make_host_uid(run_name: str, isolate_id: str) -> str:
    """Build the canonical host UID.

    Format: ``<isolate_id>_<run_name>``. Isolate ID first so that sorting
    the index files alphabetically groups all runs of the same isolate
    together.
    """
    return f"{isolate_id}_{run_name}"


def make_plasmid_uid(host_uid: str, contig: str) -> str:
    """Build the canonical plasmid UID (filesystem-safe)."""
    return f"{host_uid}__pl__{contig}"


def make_integron_uid(host_uid: str, contig: str, start: int, end: int) -> str:
    """Build the canonical integron UID (filesystem-safe)."""
    return f"{host_uid}__int__{contig}__{start}-{end}"


def sha1_of_fasta(fasta_path: Path) -> str:
    """SHA1 of the FASTA sequence content (header lines excluded).

    Used to dedupe plasmids whose assembly is byte-identical despite being
    ingested from different samples or runs.
    """
    h = hashlib.sha1()
    with open(fasta_path, "rb") as fh:
        for line in fh:
            if not line.startswith(b">"):
                h.update(line.strip())
    return h.hexdigest()


# -----------------------------------------------------------------------------
# Repository class
# -----------------------------------------------------------------------------
@dataclass
class Repository:
    """Handle to an on-disk MGE repository.

    Use ``Repository.init(path)`` to create an empty one, or
    ``Repository.open(path)`` to attach to an existing one.
    """

    root: Path

    # ---------------------------------------------------------------------
    # Path properties
    # ---------------------------------------------------------------------
    @property
    def plasmid_index_path(self) -> Path:
        return self.root / "index_plasmids.tsv"

    @property
    def integron_index_path(self) -> Path:
        return self.root / "index_integrons.tsv"

    @property
    def hosts_path(self) -> Path:
        return self.root / "hosts.tsv"

    @property
    def plasmid_dir(self) -> Path:
        return self.root / "plasmids"

    @property
    def integron_dir(self) -> Path:
        return self.root / "integrons"

    @property
    def sketches_dir(self) -> Path:
        return self.root / "sketches"

    @property
    def plasmid_sketch(self) -> Path:
        return self.sketches_dir / "plasmids.sk"

    # ---------------------------------------------------------------------
    # Lifecycle
    # ---------------------------------------------------------------------
    @classmethod
    def init(cls, root: os.PathLike | str) -> "Repository":
        """Create an empty repository at ``root`` (idempotent)."""
        repo = cls(Path(root))
        repo.root.mkdir(parents=True, exist_ok=True)
        repo.plasmid_dir.mkdir(exist_ok=True)
        repo.integron_dir.mkdir(exist_ok=True)
        repo.sketches_dir.mkdir(exist_ok=True)
        # Initialize empty index files so downstream consumers can always
        # read them. Write headers only when the file does not exist yet
        # to keep --init-repo idempotent against an existing repo.
        if not repo.plasmid_index_path.exists():
            pd.DataFrame(columns=PLASMID_INDEX_COLUMNS).to_csv(
                repo.plasmid_index_path, sep="\t", index=False
            )
        if not repo.integron_index_path.exists():
            pd.DataFrame(columns=INTEGRON_INDEX_COLUMNS).to_csv(
                repo.integron_index_path, sep="\t", index=False
            )
        if not repo.hosts_path.exists():
            pd.DataFrame(columns=HOSTS_COLUMNS).to_csv(
                repo.hosts_path, sep="\t", index=False
            )
        return repo

    @classmethod
    def open(cls, root: os.PathLike | str) -> "Repository":
        """Attach to an existing repository. Raises if structure is missing."""
        repo = cls(Path(root))
        for required in (
            repo.plasmid_index_path,
            repo.integron_index_path,
            repo.hosts_path,
        ):
            if not required.exists():
                raise FileNotFoundError(
                    f"Repository at {repo.root} is missing {required.name}. "
                    "Run with --init-repo to create the structure."
                )
        return repo

    # ---------------------------------------------------------------------
    # Read API
    # ---------------------------------------------------------------------
    def load_plasmid_index(self) -> pd.DataFrame:
        return pd.read_csv(self.plasmid_index_path, sep="\t", dtype=str).fillna("")

    def load_integron_index(self) -> pd.DataFrame:
        return pd.read_csv(self.integron_index_path, sep="\t", dtype=str).fillna("")

    def load_hosts(self) -> pd.DataFrame:
        return pd.read_csv(self.hosts_path, sep="\t", dtype=str).fillna("")

    # ---------------------------------------------------------------------
    # Write API — ingestion
    # ---------------------------------------------------------------------
    def ingest_host(self, host_row: dict) -> str:
        """Upsert one host metadata row into ``hosts.tsv``.

        Idempotent: re-ingesting the same ``host_uid`` overwrites the prior
        row in place (handy when re-running ``--resume`` flows where a
        sample's metadata gained a new column between attempts).
        Returns the ``host_uid``.
        """
        host_uid = host_row["host_uid"]
        df = self.load_hosts()
        if (df["host_uid"] == host_uid).any():
            df = df[df["host_uid"] != host_uid]
        new_row = pd.DataFrame([host_row], columns=HOSTS_COLUMNS).fillna("")
        df = pd.concat([df, new_row], ignore_index=True)
        df.to_csv(self.hosts_path, sep="\t", index=False)
        return host_uid

    def ingest_plasmid(
        self,
        host_uid: str,
        contig: str,
        fasta_src: Path,
        metadata: dict,
    ) -> Optional[str]:
        """Copy ``fasta_src`` into the repo and add an index row.

        Returns the plasmid UID on success, ``None`` if the source FASTA is
        missing or the plasmid is below :data:`MIN_PLASMID_SIZE`. Idempotent:
        re-ingesting the same UID is a no-op (no re-copy, no duplicate row).
        The skani sketch is NOT updated here — call
        :meth:`update_plasmid_sketch` once per run after all plasmids have
        been ingested to avoid quadratic sketch rebuilds.
        """
        fasta_src = Path(fasta_src)
        if not fasta_src.exists():
            return None
        try:
            size = int(metadata.get("size", 0))
        except (TypeError, ValueError):
            size = 0
        if size and size < MIN_PLASMID_SIZE:
            return None

        uid = make_plasmid_uid(host_uid, contig)
        dest = self.plasmid_dir / f"{uid}.fasta"
        if not dest.exists():
            shutil.copy(fasta_src, dest)

        df = self.load_plasmid_index()
        if (df["uid"] == uid).any():
            return uid

        sha1 = sha1_of_fasta(dest)
        row = {
            "uid": uid,
            "sha1": sha1,
            "host_uid": host_uid,
            "run_name": metadata.get("run_name", ""),
            "sample_id": metadata.get("sample_id", ""),
            "contig": contig,
            "ptu": metadata.get("ptu", ""),
            "rep": metadata.get("rep", ""),
            "mob": metadata.get("mob", ""),
            "mpf": metadata.get("mpf", ""),
            "size": str(size),
            "amr_genes": metadata.get("amr_genes", ""),
            "vir_genes": metadata.get("vir_genes", ""),
            "ingested_at": now_iso(),
        }
        new = pd.DataFrame([row], columns=PLASMID_INDEX_COLUMNS).fillna("")
        df = pd.concat([df, new], ignore_index=True)
        df.to_csv(self.plasmid_index_path, sep="\t", index=False)
        return uid

    def ingest_integron(
        self,
        host_uid: str,
        contig: str,
        start: int,
        end: int,
        assembly_fasta: Path,
        metadata: dict,
    ) -> Optional[str]:
        """Extract an integron region from ``assembly_fasta`` and add an index row.

        Coordinates are 1-based inclusive (IntegronFinder convention). The
        extracted region is written as a single-record FASTA at
        ``integrons/<uid>.fasta``. Returns the integron UID on success,
        ``None`` if the assembly FASTA is missing or the contig cannot be
        located in it. Idempotent on UID.
        """
        assembly_fasta = Path(assembly_fasta)
        uid = make_integron_uid(host_uid, contig, start, end)
        dest = self.integron_dir / f"{uid}.fasta"

        if not dest.exists():
            if not assembly_fasta.exists():
                return None
            # Lazy biopython import: heavy dependency, only required when an
            # integron is actually being ingested.
            from Bio import SeqIO

            record = None
            for rec in SeqIO.parse(str(assembly_fasta), "fasta"):
                if rec.id == contig:
                    record = rec
                    break
            if record is None:
                return None
            # 1-based inclusive → Python half-open slicing.
            sub_seq = record.seq[start - 1 : end]
            with open(dest, "w") as fh:
                fh.write(f">{uid}\n{sub_seq}\n")

        df = self.load_integron_index()
        existing_mask = df["uid"] == uid
        if existing_mask.any():
            # In-place backfill of legacy entries: every integron ingested
            # before the parser dual-format fix carried gene_set_json="[]"
            # and amr_genes="" because _parse_cassette_gene_set silently
            # dropped every cassette under the current integron_summary.csv
            # format. Re-running aluminion --resume on a historical run
            # re-invokes mge_alerts.run_ingestion which lands here with the
            # *now correctly parsed* metadata; promote it onto the existing
            # row when the existing one is empty. Guarded against
            # overwriting good data: only touches a column if the existing
            # value is empty/"[]" AND the new value is non-empty/non-"[]".
            new_gs = metadata.get("gene_set_json", "[]")
            new_amr = metadata.get("amr_genes", "")
            new_vir = metadata.get("vir_genes", "")
            existing_gs = str(df.loc[existing_mask, "gene_set_json"].iloc[0] or "")
            existing_amr = str(df.loc[existing_mask, "amr_genes"].iloc[0] or "")
            existing_vir = str(df.loc[existing_mask, "vir_genes"].iloc[0] or "")
            changed = False
            if existing_gs in ("", "[]") and new_gs not in ("", "[]"):
                df.loc[existing_mask, "gene_set_json"] = new_gs
                changed = True
            if not existing_amr and new_amr:
                df.loc[existing_mask, "amr_genes"] = new_amr
                changed = True
            if not existing_vir and new_vir:
                df.loc[existing_mask, "vir_genes"] = new_vir
                changed = True
            if changed:
                df.to_csv(self.integron_index_path, sep="\t", index=False)
            return uid

        size = end - start + 1
        row = {
            "uid": uid,
            "host_uid": host_uid,
            "run_name": metadata.get("run_name", ""),
            "sample_id": metadata.get("sample_id", ""),
            "contig": contig,
            "integron_type": metadata.get("integron_type", ""),
            "integrase": metadata.get("integrase", ""),
            "gene_set_json": metadata.get("gene_set_json", "[]"),
            "amr_genes": metadata.get("amr_genes", ""),
            "vir_genes": metadata.get("vir_genes", ""),
            "start": str(start),
            "end": str(end),
            "size": str(size),
            "ingested_at": now_iso(),
        }
        new = pd.DataFrame([row], columns=INTEGRON_INDEX_COLUMNS).fillna("")
        df = pd.concat([df, new], ignore_index=True)
        df.to_csv(self.integron_index_path, sep="\t", index=False)
        return uid

    def update_plasmid_sketch(self) -> None:
        """Rebuild the skani sketch from all FASTAs in ``plasmids/``.

        Wraps ``skani sketch`` as a subprocess. The output is a sketch
        directory at ``sketches/plasmids.sk`` that ``skani search`` can
        query in P3. If skani is not on ``$PATH`` (e.g., the user has not
        run ``mamba env update`` after pulling the new aluminion_annot
        env), this logs a warning and returns; ANI-level matching will be
        unavailable until skani is installed. Tuple-level matching (P3) is
        unaffected.
        """
        if shutil.which("skani") is None:
            print(
                "[mge_repository] skani not on PATH; ANI matching will be skipped. "
                "Update the aluminion_annot env to install it.",
                file=sys.stderr,
            )
            return
        fastas = sorted(self.plasmid_dir.glob("*.fasta"))
        if not fastas:
            return
        # Recreate the sketch directory from scratch each time so it always
        # mirrors the on-disk plasmid set exactly (cheap: skani sketches are
        # tiny and indexing is parallel).
        if self.plasmid_sketch.exists():
            shutil.rmtree(self.plasmid_sketch)
        cmd = ["skani", "sketch"] + [str(f) for f in fastas] + ["-o", str(self.plasmid_sketch)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print(
                f"[mge_repository] skani sketch failed (rc={proc.returncode}): "
                f"{proc.stderr.strip()}",
                file=sys.stderr,
            )

    # ---------------------------------------------------------------------
    # Search API — matching
    # ---------------------------------------------------------------------
    def find_plasmid_matches_by_tuple(
        self,
        ptu: str,
        rep: str,
        mob: str,
        size: int,
        size_tolerance: float = 0.15,
    ) -> pd.DataFrame:
        """Coarse plasmid match by (PTU, Rep, MOB) tuple and size band.

        Returns repo plasmid rows with the same PTU, Rep and MOB strings and
        whose stored ``size`` falls within ``size_tolerance`` fractional
        difference of the query ``size``. The size band is skipped when
        ``size`` is non-positive (caller did not supply a size).
        """
        df = self.load_plasmid_index()
        if df.empty:
            return df
        mask = (df["ptu"] == ptu) & (df["rep"] == rep) & (df["mob"] == mob)
        df_t = df[mask].copy()
        if df_t.empty or size <= 0:
            return df_t
        df_t["size"] = pd.to_numeric(df_t["size"], errors="coerce").fillna(0).astype(int)
        lo = int(size * (1 - size_tolerance))
        hi = int(size * (1 + size_tolerance))
        return df_t[(df_t["size"] >= lo) & (df_t["size"] <= hi)].copy()

    def find_plasmid_matches_by_ani(
        self,
        query_fasta: Path,
        min_ani: float = 99.0,
    ) -> pd.DataFrame:
        """Fine plasmid match by average nucleotide identity (skani).

        Returns repo plasmid rows whose ANI to ``query_fasta`` is at least
        ``min_ani`` (percent), augmented with an ``ani`` column. Falls back
        to an empty frame — without crashing — when skani is not installed
        or the sketch has not been built yet; the caller is expected to
        treat this as "no identity-level matches available" and rely on
        tuple-level matches only.
        """
        empty = pd.DataFrame(columns=list(PLASMID_INDEX_COLUMNS) + ["ani"])
        if shutil.which("skani") is None or not self.plasmid_sketch.exists():
            return empty
        query_fasta = Path(query_fasta)
        if not query_fasta.exists():
            return empty

        cmd = [
            "skani", "search",
            "-d", str(self.plasmid_sketch),
            "-q", str(query_fasta),
            "--min-af", "50",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print(
                f"[mge_repository] skani search failed (rc={proc.returncode}): "
                f"{proc.stderr.strip()}",
                file=sys.stderr,
            )
            return empty

        records: list[tuple[str, float]] = []
        for ln in proc.stdout.splitlines():
            if not ln or ln.startswith("Ref_file"):
                continue
            parts = ln.split("\t")
            if len(parts) < 3:
                continue
            try:
                ani = float(parts[2])
            except ValueError:
                continue
            if ani < min_ani:
                continue
            # The ref path is the on-disk plasmid FASTA; the UID is its stem.
            ref_uid = Path(parts[0]).stem
            records.append((ref_uid, ani))

        if not records:
            return empty

        df = self.load_plasmid_index()
        matched_uids = [r[0] for r in records]
        df_match = df[df["uid"].isin(matched_uids)].copy()
        if df_match.empty:
            return empty
        ani_map = dict(records)
        df_match["ani"] = df_match["uid"].map(ani_map).astype(float)
        return df_match.sort_values("ani", ascending=False).reset_index(drop=True)

    def find_integron_matches_by_jaccard(
        self,
        query_gene_set: Iterable[str],
        min_jaccard: float = 0.8,
    ) -> pd.DataFrame:
        """Integron match by Jaccard index on the cassette gene set.

        Returns repo integron rows whose cassette gene set has a Jaccard
        index >= ``min_jaccard`` against ``query_gene_set``, sorted by
        descending Jaccard. An empty ``query_gene_set`` or empty repo index
        returns an empty frame.

        Both the query and every repo entry are run through
        :func:`normalize_gene_set` to strip Prokka per-hit ``_<digits>``
        suffixes and drop placeholder rows ("NA", "-"). Without this,
        biologically identical integrons frequently score Jaccard < 0.5
        (see the regression test for the 5x effect on a synthetic pair),
        which masks all real epidemiological recurrences.
        """
        df = self.load_integron_index()
        result_cols = list(INTEGRON_INDEX_COLUMNS) + ["jaccard"]
        if df.empty:
            return pd.DataFrame(columns=result_cols)
        query = normalize_gene_set(query_gene_set)
        if not query:
            return pd.DataFrame(columns=result_cols)

        scores: list[float] = []
        for _, row in df.iterrows():
            # Legacy index entries written before this normalization existed
            # still carry "aadA1_5" / "NA" tokens in gene_set_json — normalize
            # on-the-fly so they match cleanly without re-ingesting the repo.
            repo_set = normalize_gene_set(deserialize_gene_set(row["gene_set_json"]))
            if not repo_set:
                scores.append(0.0)
                continue
            union = query | repo_set
            intersection = query & repo_set
            scores.append(len(intersection) / len(union) if union else 0.0)

        df = df.copy()
        df["jaccard"] = scores
        return (
            df[df["jaccard"] >= min_jaccard]
            .sort_values("jaccard", ascending=False)
            .reset_index(drop=True)
        )


# -----------------------------------------------------------------------------
# Convenience helpers
# -----------------------------------------------------------------------------
def now_iso() -> str:
    """ISO-8601 timestamp in UTC for the ``ingested_at`` column."""
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def serialize_gene_set(genes: Iterable[str]) -> str:
    """Canonical JSON serialization of an integron cassette gene set.

    Inputs are normalized first (see :func:`normalize_gene_set`) so that
    Prokka-style per-hit suffixes and "NA" placeholder rows don't pollute
    the on-disk index and silently break Jaccard matching downstream.
    """
    return json.dumps(sorted(normalize_gene_set(genes)))


def deserialize_gene_set(s: str) -> list[str]:
    """Inverse of :func:`serialize_gene_set`. Tolerates empty / null cells."""
    if not s or s in {"-", "[]"}:
        return []
    try:
        loaded = json.loads(s)
        return list(loaded) if isinstance(loaded, list) else []
    except (json.JSONDecodeError, TypeError):
        return []


# Compile once at import time so the normalize_* helpers below stay hot
# in the per-row loops inside Repository.find_integron_matches_by_jaccard.
import re as _re
_GENE_HIT_SUFFIX_RE = _re.compile(r"_\d+$")

# Placeholder cells emitted by integron_parser when a Prokka CDS had no gene
# symbol assigned (only a product description like "hypothetical protein").
# These are NOT real gene names; treating them as such inflates the Jaccard
# union and silently masks real matches between epidemiologically identical
# integrons. Add new placeholders here if integron_parser ever introduces them.
_GENE_PLACEHOLDERS = {"NA", "-", "?", "Unknown", ""}


def _normalize_gene_name(name: str) -> str:
    """Strip Prokka's per-hit ``_<digits>`` suffix (e.g. ``aadA1_5`` -> ``aadA1``).

    The suffix is a hit-index from the Prokka run, NOT a stable allele
    identifier — IntegronFinder/Prokka can assign different numbers to the
    SAME gene across two integrons in the same run, which makes set-based
    Jaccard matching essentially useless for integrons until this is stripped.
    Returns the input unchanged when no suffix is present, or ``""`` when
    the input is empty / a known placeholder.
    """
    if not isinstance(name, str):
        return ""
    s = name.strip()
    if s in _GENE_PLACEHOLDERS:
        return ""
    return _GENE_HIT_SUFFIX_RE.sub("", s)


def normalize_gene_set(genes: Iterable[str]) -> set[str]:
    """Normalize a cassette gene set: strip per-hit suffixes, drop placeholders.

    Idempotent. Used on BOTH the ingest path (so future index entries are
    clean) and the matching path's repo side (so legacy entries written
    before this normalization existed still match correctly).
    """
    out: set[str] = set()
    for g in genes:
        n = _normalize_gene_name(g)
        if n:
            out.add(n)
    return out


# -----------------------------------------------------------------------------
# Self-test entry point
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # Minimal smoke test: create a tmp repo, open it, verify file presence.
    # Real tests live in tests/test_mge_repository.py (P6).
    import argparse
    import tempfile

    ap = argparse.ArgumentParser(
        description="Smoke test for the MGE repository scaffold."
    )
    ap.add_argument(
        "--path",
        default=None,
        help="Initialize at this path (default: a fresh tmp dir).",
    )
    args = ap.parse_args()

    target = Path(args.path) if args.path else Path(tempfile.mkdtemp(prefix="aluminion_repo_"))
    repo = Repository.init(target)
    Repository.open(target)
    print(f"OK — repository initialized at {target}")
    print(f"  plasmid index : {repo.plasmid_index_path}")
    print(f"  integron idx  : {repo.integron_index_path}")
    print(f"  hosts         : {repo.hosts_path}")
