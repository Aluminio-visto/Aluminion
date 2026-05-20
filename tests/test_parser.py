#!/usr/bin/env python3
"""
Aluminion parsing pipeline tests.

Uses the example files in examples/ to verify that parser.py and
aluminion_reporter.py produce correct outputs without requiring any
databases or bioinformatics tools to be installed.

Usage (from repo root):
    python -m pytest tests/test_parser.py -v
    # or without pytest:
    python tests/test_parser.py
"""

import os
import sys
import shutil
import tempfile
import subprocess
import unittest

# Resolve repo root regardless of where the test is invoked from
REPO_ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(REPO_ROOT, "scripts")
EXAMPLES    = os.path.join(REPO_ROOT, "examples")


def build_test_dir(tmp_dir):
    """Populate a temporary run directory with the example input files."""
    # Sub-directories expected by parser.py
    os.makedirs(os.path.join(tmp_dir, "04_taxonomies", "kraken2"),  exist_ok=True)
    os.makedirs(os.path.join(tmp_dir, "04_taxonomies", "kleborate"), exist_ok=True)
    os.makedirs(os.path.join(tmp_dir, "04_taxonomies", "ectyper"),   exist_ok=True)

    copies = {
        "mlst.csv":             "mlst.csv",
        "AbR_report.csv":       "AbR_report.csv",
        "species.csv":          "04_taxonomies/kraken2/species.csv",
        "genus.csv":            "04_taxonomies/kraken2/genus.csv",
        "gambit.csv":           "04_taxonomies/gambit.csv",
        "output.tsv":           "04_taxonomies/ectyper/output.tsv",
        "enterobacterales__species_output.txt":
                                "04_taxonomies/kleborate/enterobacterales__species_output.txt",
        "phage_summary.csv":    "phage_summary.csv",
        "integron_summary.csv": "integron_summary.csv",
        "copla_modif.csv":      "copla_modif.csv",
        "kleborate.tsv":        "kleborate.tsv",
        "list_seq.tsv":         "list_seq.tsv",
        "QC_reads.csv":         "QC_reads.csv",
        "QC_assembly.csv":      "QC_assembly.csv",
    }
    for src_name, dst_rel in copies.items():
        src = os.path.join(EXAMPLES, src_name)
        dst = os.path.join(tmp_dir, dst_rel)
        if os.path.exists(src):
            shutil.copy(src, dst)


def run_script(script_name, args, cwd=None):
    """Run a Python script as a subprocess and return (returncode, stdout, stderr)."""
    cmd = [sys.executable, os.path.join(SCRIPTS_DIR, script_name)] + args
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    return result.returncode, result.stdout, result.stderr


class TestParserPy(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp(prefix="aluminion_test_")
        build_test_dir(self.tmp_dir)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def test_parser_exits_cleanly(self):
        rc, out, err = run_script("parser.py", [
            "-i", self.tmp_dir,
            "--skip-phages", "--skip-integrons", "--skip-plasmids",
        ])
        self.assertEqual(rc, 0, msg=f"parser.py exited with code {rc}.\nSTDOUT:\n{out}\nSTDERR:\n{err}")

    def test_taxonomy_csv_row_count(self):
        run_script("parser.py", [
            "-i", self.tmp_dir,
            "--skip-phages", "--skip-integrons", "--skip-plasmids",
        ])
        import pandas as pd
        tax = pd.read_csv(os.path.join(self.tmp_dir, "taxonomy.csv"))
        self.assertEqual(len(tax), 19, msg=f"taxonomy.csv has {len(tax)} rows, expected 19")

    def test_taxonomy_no_duplicates(self):
        run_script("parser.py", [
            "-i", self.tmp_dir,
            "--skip-phages", "--skip-integrons", "--skip-plasmids",
        ])
        import pandas as pd
        tax = pd.read_csv(os.path.join(self.tmp_dir, "taxonomy.csv"))
        dups = tax["Sample"].duplicated().sum()
        self.assertEqual(dups, 0, msg=f"taxonomy.csv has {dups} duplicated Sample rows")

    def test_taxonomy_key_columns_populated(self):
        run_script("parser.py", [
            "-i", self.tmp_dir,
            "--skip-phages", "--skip-integrons", "--skip-plasmids",
        ])
        import pandas as pd
        tax = pd.read_csv(os.path.join(self.tmp_dir, "taxonomy.csv"))
        for col in ["Sample", "Majority_genus", "Majority_species",
                    "Subspecies", "MLST", "KO_locus"]:
            empty = tax[col].isna().sum() + (tax[col].astype(str).str.strip() == "").sum()
            self.assertEqual(empty, 0, msg=f"Column '{col}' has {empty} empty values")

    def test_abr_file_column_clean(self):
        """#FILE in AbR_modif.xlsx must be a bare sample name without file extension."""
        run_script("parser.py", [
            "-i", self.tmp_dir,
            "--skip-phages", "--skip-integrons", "--skip-plasmids",
        ])
        import pandas as pd
        abr = pd.read_excel(os.path.join(self.tmp_dir, "AbR_modif.xlsx"))
        has_ext = abr["#FILE"].astype(str).str.contains(r"\.", na=False).any()
        self.assertFalse(has_ext, msg="#FILE column still contains file extensions")

    def test_mlst_modif_row_count(self):
        run_script("parser.py", [
            "-i", self.tmp_dir,
            "--skip-phages", "--skip-integrons", "--skip-plasmids",
        ])
        import pandas as pd
        mlst = pd.read_csv(os.path.join(self.tmp_dir, "mlst_modif.csv"), sep="\t")
        self.assertEqual(len(mlst), 19, msg=f"mlst_modif.csv has {len(mlst)} rows, expected 19")

    def test_known_resistance_genes(self):
        """OXA-48 and VIM-1 must appear in the correct samples in AbR_modif.xlsx."""
        run_script("parser.py", [
            "-i", self.tmp_dir,
            "--skip-phages", "--skip-integrons", "--skip-plasmids",
        ])
        import pandas as pd
        abr = pd.read_excel(os.path.join(self.tmp_dir, "AbR_modif.xlsx"))
        abr_idx = abr.set_index("#FILE")["Resistance_genes"].to_dict()
        self.assertIn("OXA-48",  abr_idx.get("Eclo_VC_600-1", ""), msg="OXA-48 missing from Eclo_VC_600-1")
        self.assertIn("VIM-1",   abr_idx.get("Eclo_VC_79371", ""), msg="VIM-1 missing from Eclo_VC_79371")
        self.assertIn("KPC-2",   abr_idx.get("Kpne_VC_175-1", ""), msg="KPC-2 missing from Kpne_VC_175-1")


class TestReporterPy(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp(prefix="aluminion_rep_test_")
        build_test_dir(self.tmp_dir)
        # Run parser first so reporter has all its inputs
        run_script("parser.py", [
            "-i", self.tmp_dir,
            "--skip-phages", "--skip-integrons", "--skip-plasmids",
        ])

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def test_reporter_exits_cleanly(self):
        rc, out, err = run_script("aluminion_reporter.py", [self.tmp_dir])
        self.assertEqual(rc, 0, msg=f"aluminion_reporter.py exited with code {rc}.\nSTDOUT:\n{out}\nSTDERR:\n{err}")

    def test_html_report_created(self):
        run_script("aluminion_reporter.py", [self.tmp_dir])
        html_path = os.path.join(self.tmp_dir, "Aluminion_Report.html")
        self.assertTrue(os.path.exists(html_path), msg="Aluminion_Report.html was not created")

    def test_html_contains_key_data(self):
        run_script("aluminion_reporter.py", [self.tmp_dir])
        with open(os.path.join(self.tmp_dir, "Aluminion_Report.html")) as f:
            html = f.read()
        checks = {
            "Eclo_VC_600-1": "sample name",
            "OXA-48":        "resistance gene",
            "KPC-2":         "resistance gene",
            "ecloacae":      "MLST scheme",
        }
        for value, label in checks.items():
            self.assertIn(value, html, msg=f"Expected {label} '{value}' not found in HTML report")

    def test_html_no_warnings_in_output(self):
        rc, out, err = run_script("aluminion_reporter.py", [self.tmp_dir])
        self.assertNotIn("Aviso:", out + err,
                         msg=f"Reporter emitted unexpected warnings:\n{out}\n{err}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
