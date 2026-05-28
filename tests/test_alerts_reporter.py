#!/usr/bin/env python3
"""Tests for the pattern-grouped HTML rendering in ``alerts_reporter.py``.

These tests do NOT exercise mge_alerts.emit_alerts (covered in
test_mge_alerts.py) — they consume a synthetic alerts.tsv directly and
verify the rendered HTML structure. Stream A of the 2026-05-28 UX rewrite:
one collapsible card per plasmid/integron pattern instead of one per
individual MGE, with carriers and prior hits nested inside.

Usage:
    python -m pytest tests/test_alerts_reporter.py -v
"""
import json
import os
import re
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(REPO_ROOT, "scripts")
sys.path.insert(0, SCRIPTS_DIR)

import pandas as pd  # noqa: E402

from alerts_reporter import render_alerts_html  # noqa: E402
from mge_alerts import ALERTS_COLUMNS  # noqa: E402


def _row(**overrides) -> dict:
    """Build a single alerts.tsv row with sensible defaults."""
    base = {col: "" for col in ALERTS_COLUMNS}
    base.update({
        "alert_category": "RECURRENT_RELEVANT",
        "priority": "high",
        "mge_type": "plasmid",
        "current_host_uid": "Kpne_X_2026_04_28",
        "current_isolate_id": "Kpne_X",
        "current_species": "Klebsiella pneumoniae (95%)",
        "current_mlst": "307",
        "n_hits": "1",
        "match_level": "identity",
        "cross_species": "no",
        "ptu": "PTU-L/M",
        "mob": "MOBP",
        "mpf": "typeI",
        "rep": "IncL/M(pOXA-48)",
        "size": "65530",
        "gene_set": "",
        "amr_genes": "OXA-48",
        "vir_genes": "",
        "priority_genes_detected": "OXA-48",
        "priority_categories": "carbapenemase",
        "match_hits_json": "[]",
    })
    base.update(overrides)
    return base


def _hit(**overrides) -> dict:
    """Build one prior-hit dict for embedding in match_hits_json."""
    base = {
        "match_uid": "Kpne_PRIOR_run-X__pl__AA001",
        "match_level": "identity",
        "previous_host_uid": "Kpne_PRIOR_run-X",
        "previous_isolate_id": "Kpne_PRIOR",
        "previous_seq_date": "2026_03_05",
        "previous_species": "Klebsiella pneumoniae (90%)",
        "previous_mlst": "307",
        "previous_ptu": "PTU-L/M",
        "previous_rep": "IncL/M(pOXA-48)",
        "previous_mob": "MOBP",
        "previous_mpf": "typeI",
        "previous_size": "65500",
        "previous_amr_genes": "OXA-48",
        "cross_species": "no",
        "ingested_at": "2026-03-05T08:00:00+00:00",
    }
    base.update(overrides)
    return base


def _write_alerts_tsv(rows: list, path: Path) -> None:
    df = pd.DataFrame(rows, columns=ALERTS_COLUMNS)
    df.to_csv(path, sep="\t", index=False)


class TestPatternGrouping(unittest.TestCase):
    """The collapsed card count and titles must match the user's mental model."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmp_path = Path(self.tmp.name)
        self.alerts_tsv = self.tmp_path / "alerts.tsv"
        self.html_path = self.tmp_path / "Alerts_Report.html"

    def tearDown(self):
        self.tmp.cleanup()

    def test_three_plasmid_carriers_one_pattern(self):
        """Three carriers sharing (PTU, primary Rep, top priority gene) collapse
        into ONE pattern card, not three."""
        rows = [
            _row(current_isolate_id="Kpne_A", n_hits="2",
                 match_hits_json=json.dumps([_hit(), _hit()])),
            _row(current_isolate_id="Ecol_B", current_species="Escherichia coli (88%)",
                 cross_species="yes", n_hits="1",
                 match_hits_json=json.dumps([_hit(cross_species="yes")])),
            _row(current_isolate_id="Eclo_C", current_species="Enterobacter hormaechei (80%)",
                 cross_species="yes", n_hits="3",
                 match_hits_json=json.dumps([_hit(), _hit(), _hit(cross_species="yes")])),
        ]
        _write_alerts_tsv(rows, self.alerts_tsv)

        render_alerts_html(self.alerts_tsv, self.html_path, run_name="test")
        html = self.html_path.read_text(encoding="utf-8")

        # Exactly ONE pattern <details> block at depth 1.
        pattern_details = re.findall(r"<details[^>]*class=\"pattern[^\"]*\"", html)
        self.assertEqual(len(pattern_details), 1,
                         f"expected 1 pattern card, got {len(pattern_details)}")

        # The collapsed summary should announce 3 carriers and the total prior count.
        self.assertIn("Current-run carriers", html)
        self.assertIn(">3<", html, "expected the bold carrier count")
        # 2+1+3 = 6 hits total
        self.assertIn(">6<", html, "expected the bold total prior-occurrences count")

        # Each carrier's isolate ID appears inside the body.
        for iso in ("Kpne_A", "Ecol_B", "Eclo_C"):
            self.assertIn(iso, html)

    def test_distinct_patterns_get_separate_cards(self):
        """Different (PTU, Rep, priority_gene) tuples must NOT merge."""
        rows = [
            # pOXA-48 pattern
            _row(current_isolate_id="Kpne_A"),
            # pCTX-M-15 / IncF pattern
            _row(current_isolate_id="Ecol_X",
                 current_species="Escherichia coli (90%)",
                 ptu="PTU-F", rep="IncFIA;IncFIB", mob="MOBF", mpf="typeF",
                 amr_genes="CTX-M-15",
                 priority_genes_detected="CTX-M-15",
                 priority_categories="esbl",
                 priority="medium"),
        ]
        _write_alerts_tsv(rows, self.alerts_tsv)
        render_alerts_html(self.alerts_tsv, self.html_path, run_name="test")
        html = self.html_path.read_text(encoding="utf-8")

        pattern_details = re.findall(r"<details[^>]*class=\"pattern[^\"]*\"", html)
        self.assertEqual(len(pattern_details), 2)
        # Both priority gene tokens appear as pattern titles.
        self.assertIn("OXA-48", html)
        self.assertIn("CTX-M-15", html)

    def test_first_occurrence_badge(self):
        """NEW_PRIORITY carriers inside a recurrent pattern keep a visible badge."""
        rows = [
            _row(current_isolate_id="Kpne_recurrent"),
            _row(current_isolate_id="Kpne_brandnew",
                 alert_category="NEW_PRIORITY", n_hits="0",
                 match_hits_json="[]"),
        ]
        _write_alerts_tsv(rows, self.alerts_tsv)
        render_alerts_html(self.alerts_tsv, self.html_path, run_name="test")
        html = self.html_path.read_text(encoding="utf-8")

        # One pattern card containing both carriers.
        pattern_details = re.findall(r"<details[^>]*class=\"pattern[^\"]*\"", html)
        self.assertEqual(len(pattern_details), 1)
        # The FIRST OCCURRENCE badge appears at least once.
        self.assertIn("FIRST OCCURRENCE", html)
        # The pattern-level badge counts the new occurrences.
        self.assertRegex(html, r"\b1\s+NEW\b")

    def test_high_priority_ordered_before_medium(self):
        """HIGH section must render before MEDIUM."""
        rows = [
            # MEDIUM-priority pattern with 1 carrier
            _row(current_isolate_id="Ecol_med",
                 priority="medium", amr_genes="CTX-M-15",
                 priority_genes_detected="CTX-M-15",
                 priority_categories="esbl",
                 ptu="PTU-F", rep="IncFII"),
            # HIGH-priority pattern with 2 carriers (pOXA-48)
            _row(current_isolate_id="Kpne_h1"),
            _row(current_isolate_id="Ecol_h2",
                 current_species="Escherichia coli (88%)"),
        ]
        _write_alerts_tsv(rows, self.alerts_tsv)
        render_alerts_html(self.alerts_tsv, self.html_path, run_name="test")
        html = self.html_path.read_text(encoding="utf-8")

        high_idx = html.find("HIGH priority")
        med_idx = html.find("MEDIUM priority")
        self.assertGreater(high_idx, -1, "HIGH section heading missing")
        self.assertGreater(med_idx, -1, "MEDIUM section heading missing")
        self.assertLess(high_idx, med_idx, "HIGH must render before MEDIUM")

    def test_integron_pattern_uses_cassette_gene_set(self):
        """Integrons group by the canonical (sorted) cassette gene set."""
        gs1 = json.dumps(sorted(["aadA1", "dfrA1", "sul1"]))
        gs2 = json.dumps(sorted(["dfrA1", "aadA1", "sul1"]))  # same set, diff order
        gs_other = json.dumps(sorted(["aac6", "aadA2"]))      # different set
        rows = [
            _row(mge_type="integron", current_isolate_id="Kpne_I1",
                 ptu="", mob="", mpf="", rep="",
                 gene_set=gs1, amr_genes="aadA1;dfrA1;sul1",
                 priority_genes_detected="", priority_categories="",
                 priority="medium"),
            _row(mge_type="integron", current_isolate_id="Ecol_I2",
                 current_species="Escherichia coli (88%)",
                 ptu="", mob="", mpf="", rep="",
                 gene_set=gs2, amr_genes="aadA1;dfrA1;sul1",
                 priority_genes_detected="", priority_categories="",
                 priority="medium"),
            _row(mge_type="integron", current_isolate_id="Kpne_I3",
                 ptu="", mob="", mpf="", rep="",
                 gene_set=gs_other, amr_genes="aac6;aadA2",
                 priority_genes_detected="", priority_categories="",
                 priority="medium"),
        ]
        _write_alerts_tsv(rows, self.alerts_tsv)
        render_alerts_html(self.alerts_tsv, self.html_path, run_name="test")
        html = self.html_path.read_text(encoding="utf-8")

        # gs1 and gs2 collapse, gs_other stays separate → 2 patterns.
        pattern_details = re.findall(r"<details[^>]*class=\"pattern[^\"]*\"", html)
        self.assertEqual(len(pattern_details), 2)
        # Pattern title contains the cassette gene names.
        self.assertIn("aadA1", html)
        self.assertIn("aac6", html)

    def test_plasmid_without_priority_gene_falls_into_other_bucket(self):
        """A plasmid match with no HIGH/MEDIUM gene still gets grouped by
        (PTU, primary Rep) so the user can scan it."""
        rows = [
            _row(current_isolate_id="Kpne_q1",
                 priority="medium",
                 amr_genes="", priority_genes_detected="", priority_categories="",
                 ptu="PTU-Q1", rep="ColRNAI"),
            _row(current_isolate_id="Ecol_q2",
                 current_species="Escherichia coli (88%)",
                 priority="medium",
                 amr_genes="", priority_genes_detected="", priority_categories="",
                 ptu="PTU-Q1", rep="ColRNAI"),
        ]
        _write_alerts_tsv(rows, self.alerts_tsv)
        render_alerts_html(self.alerts_tsv, self.html_path, run_name="test")
        html = self.html_path.read_text(encoding="utf-8")

        pattern_details = re.findall(r"<details[^>]*class=\"pattern[^\"]*\"", html)
        self.assertEqual(len(pattern_details), 1)
        # The pattern title should still reflect the Rep / PTU even without a gene.
        self.assertTrue(("ColRNAI" in html) or ("PTU-Q1" in html))

    def test_empty_alerts_tsv_renders_empty_banner(self):
        """No alerts → friendly banner, not a crash."""
        _write_alerts_tsv([], self.alerts_tsv)
        render_alerts_html(self.alerts_tsv, self.html_path, run_name="test")
        html = self.html_path.read_text(encoding="utf-8")
        self.assertIn("No alerts triggered", html)

    def test_match_hits_with_mixed_int_string_size(self):
        """Regression: the JSON round-trip used to leak ints for previous_size
        when the match came from the tuple branch. The reporter must render
        both happily without 'int has no attribute replace'."""
        rows = [
            _row(current_isolate_id="Kpne_A",
                 match_hits_json=json.dumps([
                     _hit(previous_size="65530"),  # quoted as string
                     {**_hit(), "previous_size": 60569},  # bare int (the bug)
                 ])),
        ]
        _write_alerts_tsv(rows, self.alerts_tsv)
        render_alerts_html(self.alerts_tsv, self.html_path, run_name="test")
        html = self.html_path.read_text(encoding="utf-8")
        self.assertIn("65530", html)
        self.assertIn("60569", html)


if __name__ == "__main__":
    unittest.main()
