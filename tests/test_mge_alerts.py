#!/usr/bin/env python3
"""Integration tests for ``scripts/mge_alerts.py``.

These exercise the end-to-end matching + relevance filter + alert
emission flow against synthetic run directories. Each test builds its
own temp run and temp repository, so the suite is fully isolated and
does NOT require skani, blast, biopython databases, or any tool from the
aluminion_annot environment beyond pandas + biopython.

Usage:
    python -m pytest tests/test_mge_alerts.py -v
"""
import argparse
import os
import sys
import shutil
import tempfile
import unittest
from pathlib import Path

REPO_ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(REPO_ROOT, "scripts")
sys.path.insert(0, SCRIPTS_DIR)

import pandas as pd  # noqa: E402

from mge_repository import Repository, make_host_uid, serialize_gene_set  # noqa: E402
from mge_alerts import _load_run_inputs, emit_alerts, run_matching  # noqa: E402


def _default_args(run_dir: Path, run_name: str, alert_new_priority: bool = False):
    """Build an argparse.Namespace mimicking what the CLI would produce."""
    return argparse.Namespace(
        run_dir=run_dir, run_name=run_name,
        ani_threshold=99.0, jaccard_threshold=0.8,
        size_tolerance=0.15, min_plasmid_size=1000,
        alert_new_priority=alert_new_priority,
        no_ingest=False, repo=None,
    )


class _SyntheticRunBuilder:
    """Helper that builds a synthetic run directory under a tmp root."""

    def __init__(self, root: Path, run_name: str):
        self.root = root
        self.run_dir = root / run_name
        self.run_dir.mkdir()
        self.run_name = run_name

    def write_data_analysis(self, rows: list[dict]) -> None:
        cols = ["Lab_id", "ID", "Barcode", "Strain", "Majority_genus",
                "Majority_species", "Subspecies", "MLST", "Serotype",
                "KO_locus", "AMRscore", "VIRscore"]
        out = "\t".join(cols) + "\n"
        for r in rows:
            out += "\t".join(str(r.get(c, "")) for c in cols) + "\n"
        (self.run_dir / "data_analysis.tsv").write_text(out)

    def write_taxonomy(self, rows: list[dict]) -> None:
        cols = ["Sample", "Majority_genus", "Majority_species", "Subspecies",
                "MLST", "Serotype", "KO_locus", "Contaminants",
                "Carbapenemase", "ESBL", "Other_resistance", "N_AMR_genes",
                "AMRscore", "VIRscore"]
        out = ",".join(cols) + "\n"
        for r in rows:
            out += ",".join(str(r.get(c, "")) for c in cols) + "\n"
        (self.run_dir / "taxonomy.csv").write_text(out)

    def write_copla(self, rows: list[dict]) -> None:
        cols = ["Sample", "Contig", "PTU", "Size", "MOB", "MPF", "Rep", "AbR"]
        out = ",".join(cols) + "\n"
        for r in rows:
            out += ",".join(str(r.get(c, "-")) for c in cols) + "\n"
        (self.run_dir / "copla_modif.csv").write_text(out)

    def write_integrons(self, rows: list[dict]) -> None:
        cols = ["Sample", "Pl/Chr", "Name", "Size", "Start", "End", "Type",
                "Integrase", "Cassette 1", "Cassette 2"]
        out = ",".join(cols) + "\n"
        for r in rows:
            out += ",".join(str(r.get(c, "")) for c in cols) + "\n"
        (self.run_dir / "integron_summary.csv").write_text(out)

    def add_plasmid_fasta(self, sample: str, contig: str, n_bp: int) -> None:
        d = self.run_dir / "08_Anotacion" / sample / "mob_recon"
        d.mkdir(parents=True, exist_ok=True)
        seq = ("ACGT" * (n_bp // 4 + 1))[:n_bp]
        (d / f"plasmid_{contig}.fasta").write_text(f">plasmid_{contig}\n{seq}\n")

    def add_assembly(self, sample: str, contig: str, n_bp: int) -> None:
        d = self.run_dir / "03_assemblies"
        d.mkdir(exist_ok=True)
        seq = ("ACGT" * (n_bp // 4 + 1))[:n_bp]
        (d / f"{sample}.fasta").write_text(f">{contig}\n{seq}\n")


def _seed_repo_with_kpne_oxa48(repo_dir: Path) -> Repository:
    """Standard prior state: one K. pneumoniae with PTU-L/M OXA-48 plasmid."""
    repo = Repository.init(repo_dir)
    host_uid = make_host_uid("run-PRIOR", "Kpne_VC_OLD")
    repo.ingest_host({
        "host_uid": host_uid, "run_name": "run-PRIOR",
        "lab_id": "10", "isolate_id": "Kpne_VC_OLD", "strain": "x",
        "genus": "Klebsiella (96%)",
        "species": "Klebsiella pneumoniae (95%)",
        "subspecies": "-", "mlst": "307", "serotype": "-",
        "ko_locus": "K19/O1", "amr_score": "3", "vir_score": "1",
        "ingested_at": "2026-03-01T00:00:00+00:00",
    })
    p = repo_dir.parent / "prior.fasta"
    p.write_text(">p\n" + ("ACGT" * 16374)[:65498] + "\n")
    repo.ingest_plasmid(
        host_uid=host_uid, contig="AA002", fasta_src=p,
        metadata={"run_name": "run-PRIOR", "sample_id": "Kpne_VC_OLD",
                  "ptu": "PTU-L/M", "rep": "IncL/M(pOXA-48)",
                  "mob": "MOBP", "mpf": "-", "size": 65498,
                  "amr_genes": "OXA-48", "vir_genes": ""},
    )
    return repo


class TestRecurrentCrossSpecies(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="aluminion_alerts_cs_"))
        _seed_repo_with_kpne_oxa48(self.tmp / "repo")
        builder = _SyntheticRunBuilder(self.tmp, "run-CURR")
        # E. coli carrying the same PTU-L/M OXA-48 plasmid: cross-species.
        builder.write_data_analysis([{"Lab_id": "1", "ID": "Eclo_VC_NEW"}])
        builder.write_taxonomy([{
            "Sample": "Eclo_VC_NEW",
            "Majority_genus": "Enterobacter (95%)",
            "Majority_species": "Enterobacter cloacae (94%)",
            "MLST": "114",
        }])
        builder.write_copla([{
            "Sample": "Eclo_VC_NEW", "Contig": "AA002",
            "PTU": "PTU-L/M", "Size": "65498", "MOB": "MOBP",
            "Rep": "IncL/M(pOXA-48)", "AbR": "OXA-48",
        }])
        builder.add_plasmid_fasta("Eclo_VC_NEW", "AA002", 65498)
        self.run_dir = builder.run_dir

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_high_priority_cross_species_alert(self):
        args = _default_args(self.run_dir, "run-CURR")
        repo = Repository.open(self.tmp / "repo")
        matches = run_matching(args, repo)
        emit_alerts(matches, args, repo)
        df = pd.read_csv(self.run_dir / "alerts.tsv", sep="\t",
                         dtype=str).fillna("")
        self.assertEqual(len(df), 1)
        row = df.iloc[0]
        self.assertEqual(row["alert_category"], "RECURRENT_RELEVANT")
        self.assertEqual(row["priority"], "high")
        self.assertEqual(row["cross_species"], "yes")
        self.assertIn("OXA-48", row["priority_genes_detected"])
        self.assertIn("Enterobacter", row["current_species"])
        self.assertIn("Klebsiella", row["previous_species"])


class TestRecurrentSameSpeciesMedium(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="aluminion_alerts_med_"))
        repo = Repository.init(self.tmp / "repo")
        prior = make_host_uid("run-PRIOR", "Kpne_VC_OLD")
        repo.ingest_host({
            "host_uid": prior, "run_name": "run-PRIOR",
            "lab_id": "10", "isolate_id": "Kpne_VC_OLD", "strain": "x",
            "genus": "Klebsiella", "species": "Klebsiella pneumoniae",
            "subspecies": "-", "mlst": "307", "serotype": "-",
            "ko_locus": "K19/O1", "amr_score": "3", "vir_score": "1",
            "ingested_at": "2026-03-01T00:00:00+00:00",
        })
        asm = self.tmp / "prior_asm.fasta"
        asm.write_text(">contig_X\n" + ("ACGT" * 7000)[:25000] + "\n")
        repo.ingest_integron(
            host_uid=prior, contig="contig_X", start=1000, end=3000,
            assembly_fasta=asm,
            metadata={"run_name": "run-PRIOR", "sample_id": "Kpne_VC_OLD",
                      "integron_type": "complete", "integrase": "tyr_intI",
                      "gene_set_json": serialize_gene_set({"aadA2", "dfrA12"}),
                      "amr_genes": "aadA2;dfrA12"},
        )
        builder = _SyntheticRunBuilder(self.tmp, "run-CURR")
        # Another K. pneumoniae with the same integron.
        builder.write_data_analysis([{"Lab_id": "1", "ID": "Kpne_VC_NEW"}])
        builder.write_taxonomy([{
            "Sample": "Kpne_VC_NEW",
            "Majority_genus": "Klebsiella",
            "Majority_species": "Klebsiella pneumoniae",
            "MLST": "15",
        }])
        builder.write_integrons([{
            "Sample": "Kpne_VC_NEW", "Pl/Chr": "contig_5",
            "Name": "int", "Size": "2000",
            "Start": "1000", "End": "2999",
            "Type": "complete", "Integrase": "tyr_intI",
            "Cassette 1": "['aadA2']", "Cassette 2": "['dfrA12']",
        }])
        builder.add_assembly("Kpne_VC_NEW", "contig_5", 4000)
        self.run_dir = builder.run_dir

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_medium_priority_no_escalation(self):
        args = _default_args(self.run_dir, "run-CURR")
        repo = Repository.open(self.tmp / "repo")
        matches = run_matching(args, repo)
        emit_alerts(matches, args, repo)
        df = pd.read_csv(self.run_dir / "alerts.tsv", sep="\t",
                         dtype=str).fillna("")
        self.assertEqual(len(df), 1)
        row = df.iloc[0]
        self.assertEqual(row["alert_category"], "RECURRENT_RELEVANT")
        self.assertEqual(row["mge_type"], "integron")
        # Same species + no priority gene -> medium.
        self.assertEqual(row["priority"], "medium")
        self.assertEqual(row["cross_species"], "no")


class TestNewPriorityFlag(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="aluminion_alerts_np_"))
        Repository.init(self.tmp / "repo")  # empty repo
        builder = _SyntheticRunBuilder(self.tmp, "run-CURR")
        builder.write_data_analysis([{"Lab_id": "1", "ID": "Eclo_X"}])
        builder.write_taxonomy([{
            "Sample": "Eclo_X", "Majority_genus": "Enterobacter",
            "Majority_species": "Enterobacter cloacae", "MLST": "114",
        }])
        # Novel plasmid carrying NDM-1 (HIGH-priority gene).
        builder.write_copla([{
            "Sample": "Eclo_X", "Contig": "P1",
            "PTU": "PTU-XYZ", "Size": "15000", "MOB": "MOBQ",
            "Rep": "IncN", "AbR": "NDM-1",
        }])
        builder.add_plasmid_fasta("Eclo_X", "P1", 15000)
        self.run_dir = builder.run_dir

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_silenced_by_default(self):
        args = _default_args(self.run_dir, "run-X", alert_new_priority=False)
        repo = Repository.open(self.tmp / "repo")
        matches = run_matching(args, repo)
        emit_alerts(matches, args, repo)
        df = pd.read_csv(self.run_dir / "alerts.tsv", sep="\t",
                         dtype=str).fillna("")
        self.assertEqual(len(df), 0,
                         "NEW_PRIORITY leaked without --alert-new-priority")

    def test_fires_with_flag(self):
        args = _default_args(self.run_dir, "run-X", alert_new_priority=True)
        repo = Repository.open(self.tmp / "repo")
        matches = run_matching(args, repo)
        emit_alerts(matches, args, repo)
        df = pd.read_csv(self.run_dir / "alerts.tsv", sep="\t",
                         dtype=str).fillna("")
        self.assertEqual(len(df), 1)
        row = df.iloc[0]
        self.assertEqual(row["alert_category"], "NEW_PRIORITY")
        self.assertEqual(row["priority"], "high")
        self.assertEqual(row["match_uid"], "")
        self.assertIn("NDM-1", row["priority_genes_detected"])


class TestRecurrenceWithoutRelevanceSilenced(unittest.TestCase):
    """A matched plasmid that carries no AMR or virulence must NOT alert."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="aluminion_alerts_silent_"))
        repo = Repository.init(self.tmp / "repo")
        prior = make_host_uid("run-PRIOR", "Kpne_VC_OLD")
        repo.ingest_host({
            "host_uid": prior, "run_name": "run-PRIOR",
            "lab_id": "10", "isolate_id": "Kpne_VC_OLD", "strain": "x",
            "genus": "Klebsiella", "species": "Klebsiella pneumoniae",
            "subspecies": "-", "mlst": "307", "serotype": "-",
            "ko_locus": "K19/O1", "amr_score": "0", "vir_score": "0",
            "ingested_at": "2026-03-01T00:00:00+00:00",
        })
        p = self.tmp / "p.fasta"
        p.write_text(">p\n" + ("ACGT" * 1214)[:4853] + "\n")
        repo.ingest_plasmid(
            host_uid=prior, contig="AA517", fasta_src=p,
            metadata={"run_name": "run-PRIOR", "sample_id": "Kpne_VC_OLD",
                      "ptu": "-", "rep": "Col440II", "mob": "MOBP",
                      "mpf": "-", "size": 4853, "amr_genes": "-",
                      "vir_genes": ""},
        )
        builder = _SyntheticRunBuilder(self.tmp, "run-CURR")
        builder.write_data_analysis([{"Lab_id": "1", "ID": "Eclo_X"}])
        builder.write_taxonomy([{
            "Sample": "Eclo_X", "Majority_genus": "Enterobacter",
            "Majority_species": "Enterobacter cloacae", "MLST": "114",
        }])
        # Same Col440II plasmid -> recurrent, but no AMR or virulence.
        builder.write_copla([{
            "Sample": "Eclo_X", "Contig": "AA517",
            "PTU": "-", "Size": "4853", "MOB": "MOBP",
            "Rep": "Col440II", "AbR": "-",
        }])
        builder.add_plasmid_fasta("Eclo_X", "AA517", 4853)
        self.run_dir = builder.run_dir

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_no_alert_emitted(self):
        args = _default_args(self.run_dir, "run-CURR")
        repo = Repository.open(self.tmp / "repo")
        matches = run_matching(args, repo)
        emit_alerts(matches, args, repo)
        df = pd.read_csv(self.run_dir / "alerts.tsv", sep="\t",
                         dtype=str).fillna("")
        self.assertEqual(len(df), 0,
                         "Matched plasmid without AMR/VIR should not alert")


class TestTaxonomyMerge(unittest.TestCase):
    """_load_run_inputs must overwrite empty data_analysis placeholders with
    populated taxonomy.csv values."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="aluminion_alerts_tx_"))
        builder = _SyntheticRunBuilder(self.tmp, "run")
        builder.write_data_analysis([{
            "Lab_id": "1", "ID": "Kpne", "Barcode": "01",
            # All taxonomy fields empty (mimics the real pre-merge state).
        }])
        builder.write_taxonomy([{
            "Sample": "Kpne",
            "Majority_genus": "Klebsiella (96%)",
            "Majority_species": "Klebsiella pneumoniae (95%)",
            "MLST": "307", "KO_locus": "K19/O1",
            "AMRscore": "3", "VIRscore": "1",
        }])
        self.run_dir = builder.run_dir

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_taxonomy_overrides_empty_data_analysis(self):
        inputs = _load_run_inputs(self.run_dir)
        da = inputs["data_analysis"]
        row = da[da["ID"] == "Kpne"].iloc[0]
        self.assertEqual(row["Majority_genus"], "Klebsiella (96%)")
        self.assertEqual(row["MLST"], "307")
        self.assertEqual(row["KO_locus"], "K19/O1")
        # data_analysis-only column preserved.
        self.assertEqual(row["Lab_id"], "1")


if __name__ == "__main__":
    unittest.main()
