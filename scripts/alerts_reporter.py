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
import json
from pathlib import Path

import pandas as pd

# At most this many prior repository hits are rendered per alert card; older
# occurrences are summarized as a count. Prevalent plasmids (pOXA-48 et al.)
# can match dozens of priors, so showing the most recent few keeps the card
# readable while still flagging recurrence depth.
_MAX_HITS_SHOWN = 5


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
.hits-header {
  margin-top: 14px;
  font-size: 0.82em;
  text-transform: uppercase;
  letter-spacing: 0.05em;
  color: var(--c-muted);
  border-top: 1px solid var(--c-border);
  padding-top: 10px;
}
.hit {
  background: #f8fafc;
  border: 1px solid var(--c-border);
  border-radius: 5px;
  padding: 10px 12px;
  margin-top: 8px;
}
.hit .row { margin-bottom: 6px; }
.hit .row:last-child { margin-bottom: 0; }
.hits-more { color: var(--c-muted); font-size: 0.85em; margin-top: 8px; }
"""


def _badge(label: str, css_class: str) -> str:
    return f'<span class="badge {css_class}">{html.escape(label)}</span>'


def _field(label: str, value: str) -> str:
    return (
        f'<div><div class="label">{html.escape(label)}</div>'
        f'<div class="value">{html.escape(value) if value else "&mdash;"}</div></div>'
    )


def _species_mlst(species: str, mlst: str) -> str:
    """Format a "species ST<n>" label, omitting MLST when absent."""
    sp = _species_short(species)
    if mlst and mlst not in ("-", "ST-"):
        return f"{sp}  ST{mlst}"
    return sp


def _render_hit(hit: dict, mge_type: str) -> str:
    """Render one prior repository occurrence as a sub-box inside an alert."""
    # Prior isolate, annotated with its sequencing date (or run name) in parens.
    isolate = hit.get("previous_isolate_id", "") or _isolate_from_uid(hit.get("previous_host_uid", ""))
    when = hit.get("previous_seq_date", "")
    isolate_label = f"{isolate} ({when})" if (isolate and when) else isolate

    fields = [
        _field("Prior isolate", isolate_label),
        _field("Prior species / MLST",
               _species_mlst(hit.get("previous_species", ""), hit.get("previous_mlst", ""))),
    ]
    if mge_type == "plasmid":
        rep_mob_mpf = (
            f"{hit.get('previous_rep') or '-'} / "
            f"{hit.get('previous_mob') or '-'} / "
            f"{hit.get('previous_mpf') or '-'}"
        )
        fields += [
            _field("Prior PTU", hit.get("previous_ptu", "")),
            _field("Prior Rep / MOB / MPF", rep_mob_mpf),
        ]
    fields += [
        _field("Prior size (bp)", hit.get("previous_size", "")),
        _field("Match level", hit.get("match_level", "")),
    ]
    if hit.get("previous_amr_genes"):
        fields.append(_field("Prior AMR genes", hit["previous_amr_genes"]))

    # Lay the fields out two-per-row, padding the last row when odd.
    rows_html = ""
    for i in range(0, len(fields), 2):
        pair = fields[i:i + 2]
        if len(pair) == 1:
            pair.append("<div></div>")
        rows_html += "<div class='row'>" + "".join(pair) + "</div>"

    cross_badge = _badge("CROSS-SPECIES", "cross") if hit.get("cross_species") == "yes" else ""
    return f"<div class='hit'>{rows_html}{cross_badge}</div>"


def _render_alert(row: pd.Series) -> str:
    """Render one alerts.tsv row as an HTML card (one alert == one MGE)."""
    category = row["alert_category"]
    priority = row["priority"]
    mge_type = row["mge_type"]
    cross = row["cross_species"] == "yes"

    head_badges: list[str] = [_badge(mge_type.upper(), f"mge-{mge_type}")]
    if category == "NEW_PRIORITY":
        head_badges.append(_badge("NEW PRIORITY", "cat-new"))
    else:
        head_badges.append(_badge("RECURRENT", "cat-rec"))
    if cross:
        head_badges.append(_badge("CROSS-SPECIES", "cross"))

    title = (
        f'{html.escape(row["current_isolate_id"])} '
        f'<span class="value">({html.escape(_species_mlst(row["current_species"], row["current_mlst"]))})</span>'
    )

    # MGE-specific metadata for the CURRENT element.
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

    # Prior repository occurrences — stacked inside this single card.
    hits_html = ""
    try:
        hits = json.loads(row["match_hits_json"]) if row["match_hits_json"] else []
    except (json.JSONDecodeError, TypeError):
        hits = []
    if hits:
        shown = hits[:_MAX_HITS_SHOWN]
        cards = "".join(_render_hit(h, mge_type) for h in shown)
        more = ""
        if len(hits) > _MAX_HITS_SHOWN:
            more = (
                f"<div class='hits-more'>&hellip; and {len(hits) - _MAX_HITS_SHOWN} "
                f"older occurrence(s) not shown.</div>"
            )
        hits_html = (
            f"<div class='hits-header'>Prior repository occurrences "
            f"({len(hits)})</div>{cards}{more}"
        )
    elif category != "NEW_PRIORITY":
        hits_html = "<div class='hits-header'>No prior repository occurrence.</div>"

    return (
        f'<div class="alert priority-{priority}">'
        f'  <h3 style="margin-bottom:8px">{title} {" ".join(head_badges)}</h3>'
        f'  {mge_block}'
        f'  {priority_genes_html}'
        f'  {hits_html}'
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
