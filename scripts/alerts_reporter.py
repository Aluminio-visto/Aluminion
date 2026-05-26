#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Script:  alerts_reporter.py
# Purpose: Render the cross-run MGE alerts table (alerts.tsv) into a
#          standalone, human-readable HTML report (Alerts_Report.html).
# Input:   <run>/alerts.tsv   (produced by mge_alerts.py)
# Output:  <run>/Alerts_Report.html
# Author:  Aluminion lab
# Date:    2026-05-21
# =============================================================================

"""HTML rendering for the MGE alerts table.

The report has two sections — HIGH priority first, then MEDIUM — each with
the alerts listed as cards (one card per row). Cards display the current
host, the matched repository host (or "NEW_PRIORITY" badge for first
occurrences), the MGE metadata, and a comma-separated list of priority
genes detected with the firing category. Pure server-side rendering: no
JavaScript dependencies, inline CSS only, so the report opens cleanly on
any browser without a network connection.
"""

from __future__ import annotations

import argparse
import datetime
import html
from pathlib import Path

import pandas as pd


# -----------------------------------------------------------------------------
# Styling — kept inline so the HTML file is fully self-contained.
# -----------------------------------------------------------------------------
_CSS = """
:root {
  --c-bg:        #fafafa;
  --c-card:      #ffffff;
  --c-border:    #e0e0e0;
  --c-text:      #1a1a1a;
  --c-muted:     #6b7280;
  --c-high:      #b91c1c;
  --c-high-bg:   #fee2e2;
  --c-medium:    #b45309;
  --c-medium-bg: #fef3c7;
  --c-cross:     #1e3a8a;
}
* { box-sizing: border-box; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
               "Helvetica Neue", Arial, sans-serif;
  background: var(--c-bg);
  color: var(--c-text);
  margin: 0;
  padding: 24px;
  line-height: 1.45;
}
h1, h2, h3 { margin-top: 0; }
.header {
  display: flex; justify-content: space-between; align-items: baseline;
  border-bottom: 2px solid var(--c-text);
  padding-bottom: 8px;
  margin-bottom: 24px;
}
.header .meta { color: var(--c-muted); font-size: 0.9em; }
.empty-banner {
  background: var(--c-card);
  border: 1px solid var(--c-border);
  border-radius: 8px;
  padding: 24px;
  text-align: center;
  color: var(--c-muted);
}
.section { margin-bottom: 32px; }
.section h2 {
  display: inline-block;
  padding: 4px 12px;
  border-radius: 4px;
  font-size: 1.1em;
}
.section.high h2   { background: var(--c-high-bg);   color: var(--c-high); }
.section.medium h2 { background: var(--c-medium-bg); color: var(--c-medium); }
.alert {
  background: var(--c-card);
  border: 1px solid var(--c-border);
  border-left: 4px solid var(--c-muted);
  border-radius: 6px;
  padding: 16px;
  margin-bottom: 12px;
}
.alert.priority-high   { border-left-color: var(--c-high); }
.alert.priority-medium { border-left-color: var(--c-medium); }
.alert .row {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 16px;
  margin-bottom: 10px;
}
.alert .label {
  color: var(--c-muted);
  font-size: 0.78em;
  text-transform: uppercase;
  letter-spacing: 0.05em;
  margin-bottom: 2px;
}
.alert .value { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
.alert .badge {
  display: inline-block;
  padding: 2px 8px;
  border-radius: 999px;
  font-size: 0.75em;
  font-weight: 600;
  margin-left: 6px;
  vertical-align: middle;
}
.badge.cat-new      { background: #dbeafe; color: #1e3a8a; }
.badge.cat-rec      { background: #ede9fe; color: #5b21b6; }
.badge.cross        { background: #dbeafe; color: var(--c-cross); }
.badge.mge-plasmid  { background: #d1fae5; color: #065f46; }
.badge.mge-integron { background: #fce7f3; color: #9d174d; }
.alert .priority-genes {
  margin-top: 6px;
  font-size: 0.9em;
}
.alert .priority-genes .pg {
  display: inline-block;
  background: var(--c-high-bg);
  color: var(--c-high);
  padding: 1px 6px;
  border-radius: 3px;
  margin-right: 4px;
  font-family: ui-monospace, monospace;
  font-size: 0.85em;
}
"""


def _badge(label: str, css_class: str) -> str:
    return f'<span class="badge {css_class}">{html.escape(label)}</span>'


def _field(label: str, value: str) -> str:
    return (
        f'<div><div class="label">{html.escape(label)}</div>'
        f'<div class="value">{html.escape(value) if value else "&mdash;"}</div></div>'
    )


def _render_alert(row: pd.Series) -> str:
    """Render one alerts.tsv row as an HTML card."""
    category = row["alert_category"]
    priority = row["priority"]
    mge_type = row["mge_type"]
    cross = row["cross_species"] == "yes"

    head_badges: list[str] = []
    head_badges.append(_badge(mge_type.upper(), f"mge-{mge_type}"))
    if category == "NEW_PRIORITY":
        head_badges.append(_badge("NEW PRIORITY", "cat-new"))
    else:
        head_badges.append(_badge("RECURRENT", "cat-rec"))
    if cross:
        head_badges.append(_badge("CROSS-SPECIES", "cross"))

    title = (
        f'{html.escape(row["current_isolate_id"])} '
        f'<span class="value">({html.escape(_species_short(row["current_species"]))})</span>'
    )

    # Match information block.
    if row["match_uid"]:
        prev_species = _species_short(row["previous_species"])
        prev_mlst = row["previous_mlst"]
        if prev_mlst and prev_mlst not in ("-", "ST-"):
            prev_species_mlst = f"{prev_species}  ST{prev_mlst}"
        else:
            prev_species_mlst = prev_species
        match_block = (
            "<div class='row'>"
            + _field("Match UID", row["match_uid"])
            + _field("Match level", row["match_level"])
            + "</div>"
            + "<div class='row'>"
            + _field("Prior isolate", _isolate_from_uid(row["previous_host_uid"]))
            + _field("Prior species / MLST", prev_species_mlst)
            + "</div>"
            + "<div class='row'>"
            + _field("Prior occurrences", row["n_prior_occurrences"])
            + _field("Prior AMR genes", row["previous_amr_genes"])
            + "</div>"
        )
    else:
        match_block = (
            "<div class='row'>"
            + _field("Match", "(new — no prior occurrence)")
            + _field("Prior occurrences", "0")
            + "</div>"
        )

    # MGE-specific metadata.
    if mge_type == "plasmid":
        rep_mob_mpf = f"{row['rep'] or '-'} / {row['mob'] or '-'} / {row['mpf'] or '-'}"
        mge_block = (
            "<div class='row'>"
            + _field("PTU", row["ptu"])
            + _field("Rep / MOB / MPF", rep_mob_mpf)
            + "</div>"
            + "<div class='row'>"
            + _field("Size (bp)", row["size"])
            + _field("AMR genes", row["amr_genes"])
            + "</div>"
        )
        if row["vir_genes"]:
            mge_block += (
                "<div class='row'>"
                f"{_field('Virulence genes', row['vir_genes'])}"
                "<div></div></div>"
            )
    else:  # integron
        mge_block = (
            "<div class='row'>"
            f"{_field('Cassette gene set', row['gene_set'])}"
            f"{_field('Size (bp)', row['size'])}"
            "</div>"
        )

    # Priority gene callouts.
    priority_genes_html = ""
    if row["priority_genes_detected"]:
        chips = "".join(
            f'<span class="pg">{html.escape(g)}</span>'
            for g in row["priority_genes_detected"].split(";")
            if g
        )
        cats = row["priority_categories"]
        priority_genes_html = (
            f"<div class='priority-genes'>"
            f"<strong>Priority genes:</strong> {chips}"
            + (f" <span class='label'>[{html.escape(cats)}]</span>" if cats else "")
            + "</div>"
        )

    return (
        f'<div class="alert priority-{priority}">'
        f'  <h3 style="margin-bottom:8px">{title} {" ".join(head_badges)}</h3>'
        f'  {match_block}'
        f'  {mge_block}'
        f'  {priority_genes_html}'
        f'</div>'
    )


def _species_short(s: str) -> str:
    """Trim the trailing ``" (NN.N%)"`` confidence suffix from data_analysis."""
    if not s:
        return ""
    return s.split(" (")[0]


def _isolate_from_uid(host_uid: str) -> str:
    """Recover the isolate ID from a host UID of the form
    ``<isolate_id>_<run_name>``. Best-effort: returns the full UID if the
    expected suffix is not present.
    """
    return host_uid or ""


def render_alerts_html(alerts_tsv: Path, output_html: Path, run_name: str = "") -> None:
    """Read ``alerts_tsv`` and write a fully-rendered HTML report to ``output_html``."""
    df = pd.read_csv(alerts_tsv, sep="\t", dtype=str).fillna("")

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n_high = (df["priority"] == "high").sum() if not df.empty else 0
    n_medium = (df["priority"] == "medium").sum() if not df.empty else 0

    header = (
        f'<div class="header">'
        f'  <div><h1>Cross-run MGE Alerts</h1>'
        f'  <div class="meta">Run: <b>{html.escape(run_name or "(unspecified)")}</b> '
        f' &middot; Generated: {timestamp}</div></div>'
        f'  <div class="meta">'
        f'    <b style="color:var(--c-high)">{n_high}</b> high &nbsp; '
        f'    <b style="color:var(--c-medium)">{n_medium}</b> medium'
        f'  </div>'
        f'</div>'
    )

    if df.empty:
        body = (
            '<div class="empty-banner">No alerts triggered for this run. '
            "Either no MGE recurrences with AMR/virulence carriage were "
            "detected, or the repository is empty (first run after init).</div>"
        )
    else:
        sections: list[str] = []
        for level, css_class, title in [
            ("high", "high", "HIGH priority"),
            ("medium", "medium", "MEDIUM priority"),
        ]:
            sub = df[df["priority"] == level]
            if sub.empty:
                continue
            cards = "\n".join(_render_alert(row) for _, row in sub.iterrows())
            sections.append(
                f'<div class="section {css_class}">'
                f'  <h2>{html.escape(title)} ({len(sub)})</h2>'
                f'  {cards}'
                f'</div>'
            )
        body = "\n".join(sections)

    full = (
        "<!DOCTYPE html>\n<html lang='en'>\n<head>\n"
        f"<meta charset='utf-8'><title>MGE Alerts — {html.escape(run_name)}</title>\n"
        f"<style>{_CSS}</style>\n</head>\n<body>\n"
        f"{header}\n{body}\n</body>\n</html>\n"
    )
    output_html.write_text(full, encoding="utf-8")


# -----------------------------------------------------------------------------
# Stand-alone invocation: useful when re-rendering the HTML without re-running
# matching / ingestion (e.g., after tweaking the report styling).
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Render alerts.tsv to Alerts_Report.html.")
    ap.add_argument("--alerts-tsv", required=True, type=Path)
    ap.add_argument("--output-html", required=True, type=Path)
    ap.add_argument("--run-name", default="")
    args = ap.parse_args()
    render_alerts_html(args.alerts_tsv, args.output_html, run_name=args.run_name)
    print(f"Wrote {args.output_html}")
