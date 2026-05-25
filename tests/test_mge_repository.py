#!/usr/bin/env python3
"""Unit tests for ``scripts/mge_repository.py``.

Covers the on-disk repository API in isolation: lifecycle, host upsert,
plasmid ingestion (filesystem-safe UIDs, size cut, idempotency), integron
region extraction with biopython, and the three search functions
(tuple match, ANI match graceful fallback, jaccard match).

Usage:
    python -m pytest tests/test_mge_repository.py -v
"""
import json
import os
import sys
import shutil
import tempfile
import unittest
from pathlib import Path

REPO_ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(REPO_ROOT, "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from mge_repository import (  # noqa: E402
    MIN_PLASMID_SIZE,
    Repository,
    deserialize_gene_set,
    make_host_uid,
    make_integron_uid,
    make_plasmid_uid,
    serialize_gene_set,
    sha1_of_fasta,
)


class TestUidConventions(unittest.TestCase):
    """UID format is part of the public contract — sorting and on-disk
    filenames depend on it."""

    def test_host_uid_isolate_first(self):
        # Isolate ID must come before the run name so alphabetical sort
        # groups runs of the same isolate together.
        self.assertEqual(
            make_host_uid("run-2026-05", "Kpne_VC_175-1"),
            "Kpne_VC_175-1_run-2026-05",
        )

    def test_plasmid_uid_uses_double_underscore(self):
        # `__` separator stays valid on NTFS, unlike `::`.
        host = "Kpne_VC_175-1_run-2026-05"
        self.assertEqual(
            make_plasmid_uid(host, "AA002"),
            "Kpne_VC_175-1_run-2026-05__pl__AA002",
        )

    def test_integron_uid_carries_coords(self):
        host = "Kpne_VC_175-1_run-2026-05"
        self.assertEqual(
            make_integron_uid(host, "contig_2", 21441, 23662),
            "Kpne_VC_175-1_run-2026-05__int__contig_2__21441-23662",
        )


class TestRepositoryLifecycle(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="aluminion_repo_test_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_init_creates_indices_and_subdirs(self):
        repo = Repository.init(self.tmp / "repo")
        for required in (
            repo.plasmid_index_path,
            repo.integron_index_path,
            repo.hosts_path,
        ):
            self.assertTrue(required.exists(), f"{required} not created")
        for subdir in (repo.plasmid_dir, repo.integron_dir, repo.sketches_dir):
            self.assertTrue(subdir.is_dir(), f"{subdir} not a directory")

    def test_init_idempotent(self):
        # A second init on a populated repo must not wipe data.
        repo = Repository.init(self.tmp / "repo")
        (repo.plasmid_dir / "marker.fasta").write_text(">x\nACGT\n")
        Repository.init(self.tmp / "repo")
        self.assertTrue((repo.plasmid_dir / "marker.fasta").exists())

    def test_open_raises_on_missing_structure(self):
        with self.assertRaises(FileNotFoundError):
            Repository.open(self.tmp / "does-not-exist")


class TestIngestion(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="aluminion_repo_ingest_"))
        self.repo = Repository.init(self.tmp / "repo")
        self.host_uid = make_host_uid("run-X", "Kpne_001")
        self.repo.ingest_host({
            "host_uid": self.host_uid, "run_name": "run-X",
            "lab_id": "1", "isolate_id": "Kpne_001", "strain": "x",
            "genus": "Klebsiella", "species": "pneumoniae",
            "subspecies": "-", "mlst": "307", "serotype": "-",
            "ko_locus": "K19/O1", "amr_score": "3", "vir_score": "1",
            "ingested_at": "2026-05-01T00:00:00+00:00",
        })

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _plasmid_fasta(self, name: str, n: int) -> Path:
        path = self.tmp / f"{name}.fasta"
        seq = ("ACGT" * (n // 4 + 1))[:n]
        path.write_text(f">{name}\n{seq}\n")
        return path

    def test_ingest_host_upserts(self):
        # Re-ingest with a different MLST: the row must be replaced, not
        # duplicated.
        self.repo.ingest_host({
            "host_uid": self.host_uid, "run_name": "run-X",
            "lab_id": "1", "isolate_id": "Kpne_001", "strain": "x",
            "genus": "Klebsiella", "species": "pneumoniae",
            "subspecies": "-", "mlst": "147",  # changed
            "serotype": "-", "ko_locus": "K19/O1",
            "amr_score": "3", "vir_score": "1",
            "ingested_at": "2026-05-02T00:00:00+00:00",
        })
        df = self.repo.load_hosts()
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["mlst"], "147")

    def test_ingest_plasmid_writes_index_and_fasta(self):
        src = self._plasmid_fasta("plasmid_AA002", 65498)
        uid = self.repo.ingest_plasmid(
            host_uid=self.host_uid, contig="AA002", fasta_src=src,
            metadata={"size": 65498, "ptu": "PTU-L/M", "rep": "IncL/M",
                      "mob": "MOBP", "amr_genes": "OXA-48"},
        )
        self.assertEqual(uid, f"{self.host_uid}__pl__AA002")
        self.assertTrue((self.repo.plasmid_dir / f"{uid}.fasta").exists())
        df = self.repo.load_plasmid_index()
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["amr_genes"], "OXA-48")
        self.assertEqual(df.iloc[0]["sha1"], sha1_of_fasta(src))

    def test_ingest_plasmid_skips_below_min_size(self):
        # Build a plasmid one byte below the minimum size threshold.
        src = self._plasmid_fasta("tiny", MIN_PLASMID_SIZE - 1)
        uid = self.repo.ingest_plasmid(
            host_uid=self.host_uid, contig="TINY", fasta_src=src,
            metadata={"size": MIN_PLASMID_SIZE - 1, "ptu": "-", "rep": "-",
                      "mob": "-", "amr_genes": "-"},
        )
        self.assertIsNone(uid)
        self.assertEqual(len(self.repo.load_plasmid_index()), 0)

    def test_ingest_plasmid_idempotent(self):
        src = self._plasmid_fasta("plasmid_AA002", 65498)
        for _ in range(3):
            self.repo.ingest_plasmid(
                host_uid=self.host_uid, contig="AA002", fasta_src=src,
                metadata={"size": 65498, "ptu": "PTU-L/M", "rep": "IncL/M",
                          "mob": "MOBP", "amr_genes": "OXA-48"},
            )
        df = self.repo.load_plasmid_index()
        self.assertEqual(len(df), 1, "Re-ingesting the same UID duplicated rows")

    def test_ingest_integron_extracts_correct_region(self):
        asm = self.tmp / "assembly.fasta"
        contig = ("ACGT" * 7000)[:25000]
        asm.write_text(f">contig_2\n{contig}\n")
        uid = self.repo.ingest_integron(
            host_uid=self.host_uid, contig="contig_2",
            start=21441, end=23662, assembly_fasta=asm,
            metadata={"integron_type": "complete", "integrase": "tyr_intI",
                      "gene_set_json": '["dfrA7_1", "emrE"]',
                      "amr_genes": "dfrA7_1;emrE"},
        )
        out = self.repo.integron_dir / f"{uid}.fasta"
        seq_only = "".join(
            ln for ln in out.read_text().splitlines() if not ln.startswith(">")
        )
        # 1-based inclusive coordinates → length = end - start + 1.
        self.assertEqual(len(seq_only), 23662 - 21441 + 1)

    def test_ingest_integron_idempotent(self):
        asm = self.tmp / "assembly.fasta"
        asm.write_text(">contig_2\n" + ("ACGT" * 6250)[:24000] + "\n")
        for _ in range(3):
            self.repo.ingest_integron(
                host_uid=self.host_uid, contig="contig_2",
                start=1000, end=2222, assembly_fasta=asm,
                metadata={"integron_type": "complete", "integrase": "tyr_intI",
                          "gene_set_json": '["aadA2"]', "amr_genes": "aadA2"},
            )
        self.assertEqual(len(self.repo.load_integron_index()), 1)


class TestMatching(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="aluminion_repo_match_"))
        self.repo = Repository.init(self.tmp / "repo")
        host = make_host_uid("run-PRIOR", "Kpne_OLD")
        self.repo.ingest_host({
            "host_uid": host, "run_name": "run-PRIOR",
            "lab_id": "10", "isolate_id": "Kpne_OLD", "strain": "x",
            "genus": "Klebsiella", "species": "pneumoniae",
            "subspecies": "-", "mlst": "307", "serotype": "-",
            "ko_locus": "K19/O1", "amr_score": "3", "vir_score": "1",
            "ingested_at": "2026-03-01T00:00:00+00:00",
        })
        # Plasmid PTU-L/M, 65498 bp.
        p = self.tmp / "p.fasta"
        p.write_text(">p\n" + ("ACGT" * 16374)[:65498] + "\n")
        self.repo.ingest_plasmid(
            host_uid=host, contig="AA002", fasta_src=p,
            metadata={"run_name": "run-PRIOR", "sample_id": "Kpne_OLD",
                      "ptu": "PTU-L/M", "rep": "IncL/M(pOXA-48)",
                      "mob": "MOBP", "mpf": "-", "size": 65498,
                      "amr_genes": "OXA-48", "vir_genes": ""},
        )
        # Integron with gene set {aadA2, dfrA12, emrE}.
        asm = self.tmp / "asm.fasta"
        asm.write_text(">c\n" + ("ACGT" * 6250)[:24000] + "\n")
        self.repo.ingest_integron(
            host_uid=host, contig="c", start=1000, end=3000, assembly_fasta=asm,
            metadata={"run_name": "run-PRIOR", "sample_id": "Kpne_OLD",
                      "integron_type": "complete", "integrase": "tyr_intI",
                      "gene_set_json": serialize_gene_set({"aadA2", "dfrA12", "emrE"}),
                      "amr_genes": "aadA2;dfrA12;emrE"},
        )

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_plasmid_tuple_match_exact(self):
        df = self.repo.find_plasmid_matches_by_tuple(
            ptu="PTU-L/M", rep="IncL/M(pOXA-48)", mob="MOBP",
            size=65000, size_tolerance=0.15,
        )
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["amr_genes"], "OXA-48")

    def test_plasmid_tuple_match_outside_size_band(self):
        # Same tuple but query is 2x larger -> outside the size band.
        df = self.repo.find_plasmid_matches_by_tuple(
            ptu="PTU-L/M", rep="IncL/M(pOXA-48)", mob="MOBP",
            size=200000, size_tolerance=0.15,
        )
        self.assertEqual(len(df), 0)

    def test_plasmid_tuple_match_different_ptu(self):
        df = self.repo.find_plasmid_matches_by_tuple(
            ptu="PTU-XYZ", rep="IncL/M(pOXA-48)", mob="MOBP",
            size=65000, size_tolerance=0.15,
        )
        self.assertEqual(len(df), 0)

    def test_plasmid_ani_graceful_fallback(self):
        # In dev (no skani), the ANI method must return an empty frame
        # rather than crash. Production CI runs the aluminion_annot env.
        q = self.tmp / "q.fasta"
        q.write_text(">q\nACGT\n")
        df = self.repo.find_plasmid_matches_by_ani(q, min_ani=99.0)
        self.assertEqual(len(df), 0)

    def test_integron_jaccard_perfect_match(self):
        df = self.repo.find_integron_matches_by_jaccard(
            {"aadA2", "dfrA12", "emrE"}, min_jaccard=0.8,
        )
        self.assertEqual(len(df), 1)
        self.assertEqual(float(df.iloc[0]["jaccard"]), 1.0)

    def test_integron_jaccard_partial_above_threshold(self):
        # Two genes shared out of three, plus one new -> 2/4 = 0.5 (below).
        df = self.repo.find_integron_matches_by_jaccard(
            {"aadA2", "dfrA12", "newGene"}, min_jaccard=0.8,
        )
        self.assertEqual(len(df), 0)
        # But the same two with no new -> 2/3 = 0.66 (still below 0.8).
        df = self.repo.find_integron_matches_by_jaccard(
            {"aadA2", "dfrA12"}, min_jaccard=0.8,
        )
        self.assertEqual(len(df), 0)

    def test_integron_jaccard_empty_query(self):
        df = self.repo.find_integron_matches_by_jaccard(set(), min_jaccard=0.8)
        self.assertEqual(len(df), 0)


class TestGeneSetSerialization(unittest.TestCase):
    def test_round_trip(self):
        genes = {"dfrA12", "aadA2", "emrE"}
        s = serialize_gene_set(genes)
        # Canonical sorted order, deterministic for diffing.
        self.assertEqual(json.loads(s), sorted(genes))
        self.assertEqual(set(deserialize_gene_set(s)), genes)

    def test_deserialize_empty_inputs(self):
        for empty in ("", "-", "[]", None):
            self.assertEqual(deserialize_gene_set(empty), [])


if __name__ == "__main__":
    unittest.main()
