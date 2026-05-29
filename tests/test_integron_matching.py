#!/usr/bin/env python3
"""Regression tests for the integron alerting path.

Three classes of bug led to ZERO integron alerts on the 2026-05-13
production batch despite known epidemiological recurrences:

  Bug A (gene set pollution): "NA" placeholder rows survived parsing and
        inflated the union, lowering Jaccard below the 0.8 threshold.
  Bug B (Prokka per-hit suffix): "aadA1_1" vs "aadA1_5" treated as
        different genes; the trailing "_<digits>" is a hit counter, not
        a stable allele identifier.
  Bug D (parser format mismatch): _parse_cassette_gene_set expected a
        Python list literal ("['blaOXA-1_1']") but integron_parser.py
        had switched to writing bare comma-separated strings
        ("blaOXA-1, NA"). ast.literal_eval silently failed; every
        cassette was dropped; every repo entry was ingested with
        gene_set_json="[]". 80 repo entries on the user's production
        server were biologically useless until both this fix and the
        ingest_integron backfill landed.

Usage:
    python -m pytest tests/test_integron_matching.py -v
"""
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(REPO_ROOT, "scripts")
sys.path.insert(0, SCRIPTS_DIR)

import pandas as pd  # noqa: E402

from mge_alerts import _parse_cassette_gene_set  # noqa: E402
from mge_repository import (  # noqa: E402
    Repository, make_host_uid, normalize_gene_set, serialize_gene_set,
)


class TestCassetteParser(unittest.TestCase):
    """The parser must handle BOTH formats produced by integron_parser.py
    across its lifetime, applying the normalization defined in
    normalize_gene_set."""

    def test_current_format_single_bare_token(self):
        """Single-gene cell as a bare string, no brackets, no quotes."""
        row = pd.Series({"Cassette 1": "blaGES-13"})
        self.assertEqual(_parse_cassette_gene_set(row), {"blaGES-13"})

    def test_current_format_comma_separated_with_NA(self):
        """Multi-gene cell as 'gene1, NA, NA, NA' — NA placeholders dropped."""
        row = pd.Series({
            "Cassette 1": "aac(6')-Ib-D181Y",
            "Cassette 2": "blaOXA-1, NA, NA, NA",
        })
        self.assertEqual(
            _parse_cassette_gene_set(row),
            {"aac(6')-Ib-D181Y", "blaOXA-1"},
        )

    def test_legacy_format_list_literal(self):
        """Backward compat: '['blaOXA-1_1']' list literal still parses."""
        row = pd.Series({"Cassette 1": "['blaOXA-1_1']"})
        # _1 hit suffix stripped by normalize.
        self.assertEqual(_parse_cassette_gene_set(row), {"blaOXA-1"})

    def test_legacy_format_gene_semicolon_description(self):
        """Legacy 'NA;hypothetical protein' rows are entirely filtered."""
        row = pd.Series({"Cassette 1": "['NA;hypothetical protein']"})
        self.assertEqual(_parse_cassette_gene_set(row), set())

    def test_prokka_hit_suffix_stripped(self):
        """Bug B: 'aadA1_1' and 'aadA1_5' must both normalize to 'aadA1'."""
        row1 = pd.Series({"Cassette 1": "aadA1_1"})
        row5 = pd.Series({"Cassette 1": "aadA1_5"})
        self.assertEqual(_parse_cassette_gene_set(row1),
                         _parse_cassette_gene_set(row5))
        self.assertEqual(_parse_cassette_gene_set(row1), {"aadA1"})

    def test_empty_and_placeholder_cells_drop(self):
        row = pd.Series({"Cassette 1": "", "Cassette 2": "-", "Cassette 3": "NA"})
        self.assertEqual(_parse_cassette_gene_set(row), set())

    def test_two_integrons_identical_biology_match_jaccard_one(self):
        """The headline regression: two integrons whose cassettes differ
        only by the per-hit suffix or NA pollution must score Jaccard 1.0
        after parsing."""
        old = pd.Series({"Cassette 1": "['blaOXA-1_1', 'aadA1_5']",
                         "Cassette 2": "['NA;hypothetical protein']"})
        new = pd.Series({"Cassette 1": "blaOXA-1, aadA1, NA"})
        gs_old = _parse_cassette_gene_set(old)
        gs_new = _parse_cassette_gene_set(new)
        self.assertEqual(gs_old, gs_new)
        union = gs_old | gs_new
        inter = gs_old & gs_new
        self.assertEqual(len(inter) / len(union), 1.0)


class TestNormalizeGeneSet(unittest.TestCase):
    """normalize_gene_set is the single chokepoint shared by ingest and match."""

    def test_strip_hit_suffix(self):
        self.assertEqual(normalize_gene_set(["aadA1_1", "aadA1_5"]),
                         {"aadA1"})

    def test_drop_placeholders(self):
        self.assertEqual(normalize_gene_set(["NA", "-", "?", "blaOXA-1"]),
                         {"blaOXA-1"})

    def test_idempotent(self):
        once = normalize_gene_set(["aadA1_1", "NA", "sul1"])
        twice = normalize_gene_set(once)
        self.assertEqual(once, twice)

    def test_preserves_gene_internal_digits(self):
        """dfrA12 ends in '2' but has no underscore — must stay intact."""
        self.assertEqual(normalize_gene_set(["dfrA12"]), {"dfrA12"})
        # Only the trailing _digits is stripped.
        self.assertEqual(normalize_gene_set(["dfrA12_1"]), {"dfrA12"})

    def test_serialize_runs_through_normalize(self):
        """serialize_gene_set must produce a clean JSON so the on-disk
        index never carries the noisy tokens."""
        out = serialize_gene_set(["aadA1_1", "NA", "blaOXA-1_3"])
        self.assertEqual(json.loads(out), ["aadA1", "blaOXA-1"])


class TestIngestIntegronBackfill(unittest.TestCase):
    """When the repository already holds an integron entry with empty
    gene_set_json (the broken-parser legacy state on production), a
    subsequent ingest of the SAME UID with a populated gene_set must
    UPDATE the existing row instead of silently returning early.

    This is what makes `aluminion --resume` on a historical run fix the
    biologically empty entries without requiring a separate backfill script.
    """

    def setUp(self):
        # mkdtemp + rmtree(ignore_errors=True) mirrors test_mge_repository.py;
        # TemporaryDirectory.cleanup() trips on a Windows + biopython interaction
        # where SeqIO.parse() leaves the FASTA handle open after `break` and
        # the tempdir teardown then can't delete it.
        self.tmp = tempfile.mkdtemp()
        self.repo_dir = Path(self.tmp) / "repo"
        self.repo = Repository.init(self.repo_dir)
        self.assembly = Path(self.tmp) / "asm.fasta"
        # 2 kb synthetic contig that covers our integron coordinates.
        self.assembly.write_text(">contig_1\n" + "ACGT" * 500 + "\n")
        self.host_uid = make_host_uid("run-2026_03_05", "Kpne_OLD")
        self.repo.ingest_host({
            "host_uid": self.host_uid, "run_name": "run-2026_03_05",
            "lab_id": "1", "isolate_id": "Kpne_OLD", "strain": "x",
            "genus": "Klebsiella", "species": "K. pneumoniae",
            "subspecies": "-", "mlst": "11", "serotype": "-",
            "ko_locus": "-", "amr_score": "0", "vir_score": "0",
            "ingested_at": "2026-03-05T08:00:00+00:00",
        })

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _ingest(self, gene_set_json: str, amr_genes: str = ""):
        return self.repo.ingest_integron(
            self.host_uid, "contig_1", 100, 1100, self.assembly,
            metadata={
                "run_name": "run-2026_03_05",
                "sample_id": "Kpne_OLD",
                "integron_type": "complete",
                "integrase": "intI",
                "gene_set_json": gene_set_json,
                "amr_genes": amr_genes,
                "vir_genes": "",
            },
        )

    def _read_index(self) -> pd.DataFrame:
        return pd.read_csv(self.repo.integron_index_path, sep="\t").fillna("")

    def test_empty_then_populated_backfills(self):
        """Initial ingest with [] gene_set; re-ingest with real genes
        updates the existing row in place."""
        uid1 = self._ingest("[]", amr_genes="")
        self.assertIsNotNone(uid1)
        df = self._read_index()
        self.assertEqual(df.loc[df["uid"] == uid1, "gene_set_json"].iloc[0], "[]")

        # Now simulate the parser-fix re-run with the correctly parsed cassettes.
        uid2 = self._ingest(json.dumps(["aadA1", "sul1"]),
                            amr_genes="aadA1;sul1")
        self.assertEqual(uid1, uid2)
        df = self._read_index()
        self.assertEqual(df.loc[df["uid"] == uid2, "gene_set_json"].iloc[0],
                         '["aadA1", "sul1"]')
        self.assertEqual(df.loc[df["uid"] == uid2, "amr_genes"].iloc[0],
                         "aadA1;sul1")

    def test_populated_then_empty_does_not_overwrite(self):
        """A re-ingest carrying an empty gene_set must NOT clobber the
        already-populated value (guard against bad inputs damaging the repo)."""
        uid1 = self._ingest(json.dumps(["aadA1"]), amr_genes="aadA1")
        # Re-ingest with empty payload (e.g. a broken re-parse).
        self._ingest("[]", amr_genes="")
        df = self._read_index()
        self.assertEqual(df.loc[df["uid"] == uid1, "gene_set_json"].iloc[0],
                         '["aadA1"]')
        self.assertEqual(df.loc[df["uid"] == uid1, "amr_genes"].iloc[0],
                         "aadA1")

    def test_populated_stable_on_second_ingest_with_same_data(self):
        """Idempotency: re-ingesting the same row doesn't churn the file."""
        gs = json.dumps(["aadA1", "sul1"])
        self._ingest(gs, amr_genes="aadA1;sul1")
        mtime1 = self.repo.integron_index_path.stat().st_mtime_ns
        self._ingest(gs, amr_genes="aadA1;sul1")
        mtime2 = self.repo.integron_index_path.stat().st_mtime_ns
        # No write happens when nothing would change.
        self.assertEqual(mtime1, mtime2)


class TestJaccardMatchesRealRecurrence(unittest.TestCase):
    """End-to-end: ingest one historical integron, then look up a current
    integron that carries the same biology under both the legacy and
    current formats. Must score Jaccard 1.0 and surface as a match."""

    def setUp(self):
        # Same Windows teardown workaround as TestIngestIntegronBackfill above.
        self.tmp = tempfile.mkdtemp()
        self.repo_dir = Path(self.tmp) / "repo"
        self.repo = Repository.init(self.repo_dir)
        self.assembly = Path(self.tmp) / "asm.fasta"
        self.assembly.write_text(">contig_1\n" + "ACGT" * 500 + "\n")
        host_uid = make_host_uid("run-2026_03_05", "Kpne_OLD")
        self.repo.ingest_host({
            "host_uid": host_uid, "run_name": "run-2026_03_05",
            "lab_id": "1", "isolate_id": "Kpne_OLD", "strain": "x",
            "genus": "Klebsiella", "species": "K. pneumoniae",
            "subspecies": "-", "mlst": "11", "serotype": "-",
            "ko_locus": "-", "amr_score": "0", "vir_score": "0",
            "ingested_at": "2026-03-05T08:00:00+00:00",
        })
        # Prior integron written under the LEGACY format (list literal with
        # "_N" suffixes), serialized through the new normalize-aware
        # serialize_gene_set so it lands clean in the index.
        self.repo.ingest_integron(
            host_uid, "contig_1", 100, 2000, self.assembly,
            metadata={
                "run_name": "run-2026_03_05",
                "sample_id": "Kpne_OLD",
                "integron_type": "complete",
                "integrase": "intI",
                "gene_set_json": serialize_gene_set(["aadA1_1", "sul1_1", "NA"]),
                "amr_genes": "aadA1;sul1",
                "vir_genes": "",
            },
        )

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_matches_current_format_with_NA_pollution(self):
        """Today's current-format cassette with 'NA' noise matches the
        repo's older entry once both go through normalize_gene_set."""
        row = pd.Series({"Cassette 1": "aadA1, sul1, NA"})
        query = _parse_cassette_gene_set(row)
        self.assertEqual(query, {"aadA1", "sul1"})
        df = self.repo.find_integron_matches_by_jaccard(query, min_jaccard=0.8)
        self.assertEqual(len(df), 1, "expected exactly one Jaccard match")
        self.assertEqual(df["jaccard"].iloc[0], 1.0)

    def test_legacy_repo_entry_with_unstripped_genes_still_matches(self):
        """Defensive: even if a legacy entry slipped into the repo with
        un-normalized gene_set_json (e.g. ingested before this fix),
        find_integron_matches_by_jaccard normalizes the repo side at
        query time so the match still fires."""
        # Patch the index entry to simulate a legacy un-normalized record.
        df = pd.read_csv(self.repo.integron_index_path, sep="\t").fillna("")
        df["gene_set_json"] = json.dumps(["aadA1_5", "sul1_3", "NA"])
        df.to_csv(self.repo.integron_index_path, sep="\t", index=False)

        row = pd.Series({"Cassette 1": "aadA1, sul1"})
        query = _parse_cassette_gene_set(row)
        out = self.repo.find_integron_matches_by_jaccard(query, min_jaccard=0.8)
        self.assertEqual(len(out), 1)
        self.assertEqual(out["jaccard"].iloc[0], 1.0)


if __name__ == "__main__":
    unittest.main()
