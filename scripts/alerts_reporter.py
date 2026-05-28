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

# Priority tier ordering for the pattern grouping logic below. The user-facing
# pattern title prefers the gene from the highest-tier category present on the
# carriers — so a plasmid carrying OXA-48 (carbapenemase) is named after the
# carbapenemase even if the alert row also lists an ESBL.
_PRIORITY_TIER_ORDER = ("carbapenemase", "colistin", "hv_klebsiella", "esbl")

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


def _field(label, value) -> str:
    """Render one label/value pair as a small HTML block.

    Coerces every argument to ``str`` before passing it to ``html.escape``
    because ``html.escape`` calls ``s.replace('&', '&amp;')`` and dies on
    non-string scalars with ``'int' object has no attribute 'replace'``.
    The mge_alerts JSON round-trip used to leak ints when ``repo_row['size']``
    came from the tuple-match branch (cast via ``.astype(int)``); fixing that
    at the source removed the immediate trigger, but keeping this defensive
    here means a future column whose dtype drifts cannot crash the report.
    """
    label_s = str(label) if label is not None else ""
    if value is None or value == "":
        value_html = "&mdash;"
    else:
        value_html = html.escape(str(value))
    return (
        f'<div><div class="label">{html.escape(label_s)}</div>'
        f'<div class="value">{value_html}</div></div>'
    )


def _species_mlst(species: str, mlst: str) -> str:
    """Format a "species ST<n>" label, omitting MLST when absent."""
    sp = _species_short(species)
    if mlst and mlst not in ("-", "ST-"):
        return f"{sp}  ST{mlst}"
    return sp


# =============================================================================
# Pattern grouping (Stream A of the 2026-05-28 UX rewrite)
# =============================================================================
# The alerts.tsv schema stays one-row-per-MGE so external consumers don't
# break, but the HTML report collapses N alert rows sharing a biological
# "pattern" (e.g. pOXA-48-like = PTU-L/M + IncL/M + OXA-48) into a single
# expandable card. Carriers and their prior hits live inside that card.
# Pattern key:
#   plasmid  -> (PTU, primary Rep token, top priority gene)
#   integron -> ('integron', canonical sorted cassette gene set)
# Plasmids without any priority gene fall into an "Other" bucket keyed only
# by (PTU, primary Rep) so they still cluster by tuple instead of fragmenting.
# =============================================================================

def _primary_rep(rep: str) -> str:
    """First ';'-separated Rep token; '' if none."""
    if not rep:
        return ""
    return rep.split(";")[0].strip()


def _top_priority_gene(row: pd.Series) -> tuple[str, str]:
    """Pick the (gene, category) used to title the pattern.

    Returns the first gene of the highest-tier category present in
    ``priority_categories``; falls back to the first gene of
    ``priority_genes_detected`` if neither maps to a known tier. Empty
    string when no priority gene was detected.
    """
    detected = (row.get("priority_genes_detected", "") or "").strip()
    if not detected:
        return ("", "")
    categories = {c.strip() for c in (row.get("priority_categories", "") or "").split(",") if c.strip()}
    genes = [g.strip() for g in detected.split(";") if g.strip()]
    if not genes:
        return ("", "")
    # Pick the tier that's actually present on this row, by global tier order.
    for tier in _PRIORITY_TIER_ORDER:
        if tier in categories:
            # Match the first gene whose pattern matches this tier. We avoid
            # re-importing _priority_genes' compiled patterns here; the row's
            # `priority_genes_detected` is already the matched set, so the
            # tier label + a gene from that set is good enough as a title.
            return (genes[0], tier)
    return (genes[0], next(iter(categories), ""))


def _pattern_key(row: pd.Series) -> tuple:
    """Compute the hashable pattern key for one alerts.tsv row."""
    mge_type = row.get("mge_type", "")
    if mge_type == "integron":
        gs_json = row.get("gene_set", "") or ""
        try:
            genes = json.loads(gs_json) if gs_json else []
        except (json.JSONDecodeError, TypeError):
            genes = []
        canonical = tuple(sorted({str(g).strip() for g in genes if str(g).strip()}))
        return ("integron", canonical)
    # plasmid (and any unknown mge_type) keyed by tuple.
    ptu = (row.get("ptu", "") or "").strip()
    rep_primary = _primary_rep(row.get("rep", "") or "")
    top_gene, _ = _top_priority_gene(row)
    return ("plasmid", ptu, rep_primary, top_gene)


def _pattern_title(key: tuple, sample_row: pd.Series) -> str:
    """Build the human-readable title shown on the pattern card."""
    if key[0] == "integron":
        _, genes = key
        if not genes:
            return "Integron (empty cassette)"
        shown = list(genes)[:3]
        more = f" + {len(genes) - 3} more" if len(genes) > 3 else ""
        return f"Integron: {', '.join(shown)}{more}"
    # plasmid
    _, ptu, rep, gene = key
    # Treat the bare dash as "not determined" — appending it as "(-)" is just
    # noise in the title. The metadata line under the title carries the raw
    # PTU value for users who care about the difference.
    ptu_show = ptu if (ptu and ptu != "-") else ""
    rep_show = rep if (rep and rep != "-") else ""
    if gene:
        # "OXA-48 / IncL/M (PTU-L/M)" — the gene comes first because that's
        # what the clinician is searching the page for; the PTU goes in
        # parens so two patterns sharing gene+Rep but differing by PTU don't
        # render with identical titles (real case on the production run:
        # CTX-M-15 / IncFIB(K) split across PTU-F and PTU-FE).
        title = f"{gene} / {rep_show}" if rep_show else gene
        if ptu_show:
            title += f" ({ptu_show})"
        return title
    if rep_show:
        return f"Plasmid: {rep_show}" + (f" ({ptu_show})" if ptu_show else "")
    if ptu_show:
        return f"Plasmid: {ptu_show}"
    return "Plasmid: unclassified"


def _pattern_summary_meta(block: dict) -> str:
    """One-line "metadata strip" rendered under the pattern title.

    Aggregates the carriers' shared metadata into a compact summary. Per-
    carrier variation (different MOB / MPF strings detected across samples)
    is exposed in each carrier sub-block when the card is expanded.
    """
    parts: list[str] = []
    if block["ptu"]:
        parts.append(f"PTU: {block['ptu']}")
    if block["rep"] or block["mob"] or block["mpf"]:
        parts.append(
            f"Rep / MOB / MPF: {block['rep'] or '-'} / "
            f"{block['mob'] or '-'} / {block['mpf'] or '-'}"
        )
    if block["size_repr"]:
        parts.append(f"~{block['size_repr']} bp")
    if block["top_gene"]:
        cat = f" [{block['top_category']}]" if block["top_category"] else ""
        parts.append(f"Priority: {block['top_gene']}{cat}")
    return " · ".join(parts)


def _representative_value(values: list[str]) -> str:
    """Most common non-empty value across carriers; '' when all blank."""
    cleaned = [v for v in values if v and v != "-"]
    if not cleaned:
        return ""
    # Counter-by-hand to avoid a stdlib import cost; the list is small.
    counts: dict[str, int] = {}
    for v in cleaned:
        counts[v] = counts.get(v, 0) + 1
    return max(counts.items(), key=lambda x: x[1])[0]


def _group_by_pattern(df: pd.DataFrame) -> list[dict]:
    """Group an alerts.tsv DataFrame into pattern blocks.

    Returns a list of dicts shaped for ``_render_pattern_card``. The caller
    splits the list by priority and orders within each section.
    """
    blocks: dict[tuple, dict] = {}
    for _, row in df.iterrows():
        key = _pattern_key(row)
        blk = blocks.get(key)
        if blk is None:
            top_gene, top_cat = _top_priority_gene(row)
            blk = {
                "key": key,
                "mge_type": row.get("mge_type", ""),
                "carriers": [],
                "priorities": set(),
                "cross_species_carriers": 0,
                "first_occurrence_carriers": 0,
                "total_prior_hits": 0,
                "top_gene": top_gene,
                "top_category": top_cat,
                # Representative metadata collected after the loop.
                "_ptu_vals": [],
                "_rep_vals": [],
                "_mob_vals": [],
                "_mpf_vals": [],
                "_size_vals": [],
            }
            blocks[key] = blk
        blk["carriers"].append(row)
        blk["priorities"].add(row.get("priority", "") or "")
        if (row.get("cross_species", "") or "") == "yes":
            blk["cross_species_carriers"] += 1
        if (row.get("alert_category", "") or "") == "NEW_PRIORITY":
            blk["first_occurrence_carriers"] += 1
        try:
            blk["total_prior_hits"] += int(row.get("n_hits", "0") or 0)
        except (TypeError, ValueError):
            pass
        blk["_ptu_vals"].append(row.get("ptu", "") or "")
        blk["_rep_vals"].append(_primary_rep(row.get("rep", "") or ""))
        blk["_mob_vals"].append(row.get("mob", "") or "")
        blk["_mpf_vals"].append(row.get("mpf", "") or "")
        blk["_size_vals"].append(row.get("size", "") or "")
        # Prefer to remember a top gene/category if the first carrier didn't
        # have one but a later carrier in the same key does (rare edge case;
        # the key would force them to share, so a non-empty top_gene on any
        # carrier is the right value to display).
        if not blk["top_gene"]:
            g, c = _top_priority_gene(row)
            if g:
                blk["top_gene"], blk["top_category"] = g, c

    # Finalise representative fields and pattern-level metadata.
    final: list[dict] = []
    for key, blk in blocks.items():
        sample_row = blk["carriers"][0]
        blk["ptu"] = _representative_value(blk["_ptu_vals"])
        blk["rep"] = _representative_value(blk["_rep_vals"])
        blk["mob"] = _representative_value(blk["_mob_vals"])
        blk["mpf"] = _representative_value(blk["_mpf_vals"])
        # Size: numeric median of the carriers; falls back to representative
        # string if nothing parses.
        numeric_sizes: list[int] = []
        for s in blk["_size_vals"]:
            try:
                if s:
                    numeric_sizes.append(int(float(s)))
            except (TypeError, ValueError):
                pass
        if numeric_sizes:
            numeric_sizes.sort()
            mid = numeric_sizes[len(numeric_sizes) // 2]
            blk["size_repr"] = f"{mid:,}"
        else:
            blk["size_repr"] = _representative_value(blk["_size_vals"])
        # Pattern priority = the highest priority any carrier reaches.
        blk["priority"] = "high" if "high" in blk["priorities"] else (
            "medium" if "medium" in blk["priorities"] else ""
        )
        blk["any_cross_species"] = blk["cross_species_carriers"] > 0
        blk["n_carriers"] = len(blk["carriers"])
        blk["title"] = _pattern_title(key, sample_row)
        blk["meta_line"] = _pattern_summary_meta(blk)
        # Drop the working buffers so the dict is clean for rendering.
        for k in list(blk.keys()):
            if k.startswith("_"):
                del blk[k]
        final.append(blk)
    return final


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


def _render_carrier(row: pd.Series) -> str:
    """Render one alerts.tsv row as a carrier sub-card inside a pattern.

    Was previously the top-level ``_render_alert``; renamed for Stream A
    of the 2026-05-28 UX rewrite, where the top level is now the pattern
    card and individual MGE rows render here. Other than the FIRST
    OCCURRENCE badge (replacing the old NEW PRIORITY one with semantics
    closer to what the user reads in a grouped report), the body is
    unchanged on purpose so the per-carrier detail view stays familiar.
    """
    category = row["alert_category"]
    priority = row["priority"]
    mge_type = row["mge_type"]
    cross = row["cross_species"] == "yes"

    head_badges: list[str] = [_badge(mge_type.upper(), f"mge-{mge_type}")]
    if category == "NEW_PRIORITY":
        # Inside a recurrent pattern, this carrier is the first time the MGE
        # appears in the repository. Tag it so the user spots novel hosts
        # at a glance without having to scan match_hits_json counts.
        head_badges.append(_badge("FIRST OCCURRENCE", "first-occ"))
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


def _render_pattern_card(block: dict) -> str:
    """Render one pattern as a collapsible <details> card with carriers inside."""
    head_badges: list[str] = [_badge(block["mge_type"].upper(), f"mge-{block['mge_type']}")]
    if block["priority"]:
        # Reuse the cat-rec colors for the priority badge to keep the existing
        # palette small; the priority tier is already visible via the left
        # border color and the section heading.
        head_badges.append(_badge(f"{block['priority'].upper()} PRIORITY", "cat-rec"))
    if block["any_cross_species"]:
        head_badges.append(_badge("CROSS-SPECIES", "cross"))
    if block["first_occurrence_carriers"] > 0:
        head_badges.append(_badge(f"{block['first_occurrence_carriers']} NEW", "first-occ"))

    counts = (
        f"Current-run carriers: <b>{block['n_carriers']}</b> &middot; "
        f"Prior occurrences (sum across carriers): <b>{block['total_prior_hits']}</b>"
    )

    carriers_html = "\n".join(_render_carrier(r) for r in block["carriers"])

    summary = (
        f'<summary>'
        f'<h3>{html.escape(block["title"])} {" ".join(head_badges)}</h3>'
        f'<div class="pattern-meta">{html.escape(block["meta_line"])}</div>'
        f'<div class="pattern-counts">{counts}</div>'
        f'</summary>'
    )
    return (
        f'<details class="pattern priority-{block["priority"]}">'
        f'{summary}'
        f'<div class="pattern-body">{carriers_html}</div>'
        f'</details>'
    )


def render_alerts_html(alerts_tsv: Path, output_html: Path, run_name: str = "") -> None:
    """Read ``alerts_tsv`` and write a fully-rendered HTML report to ``output_html``.

    The report groups alerts by biological pattern (Stream A of the
    2026-05-28 UX rewrite). One ``<details>`` card per pattern; carriers
    and their prior hits render inside when expanded. Pattern key:
    ``(PTU, primary Rep, top priority gene)`` for plasmids, the canonical
    sorted cassette gene set for integrons.
    """
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
        # Build pattern blocks once; sections filter by priority. Within a
        # section, order patterns by (number of current-run carriers DESC,
        # number of prior hits DESC) so the most prevalent show up first.
        all_blocks = _group_by_pattern(df)
        all_blocks.sort(
            key=lambda b: (b["n_carriers"], b["total_prior_hits"]),
            reverse=True,
        )

        sections: list[str] = []
        for level, css_class, title in [
            ("high", "high", "HIGH priority"),
            ("medium", "medium", "MEDIUM priority"),
        ]:
            blocks = [b for b in all_blocks if b["priority"] == level]
            if not blocks:
                continue
            cards = "\n".join(_render_pattern_card(b) for b in blocks)
            n_carriers = sum(b["n_carriers"] for b in blocks)
            sections.append(
                f'<div class="section {css_class}">'
                f'  <h2>{html.escape(title)} '
                f'({len(blocks)} pattern{"s" if len(blocks) != 1 else ""}, '
                f'{n_carriers} carrier{"s" if n_carriers != 1 else ""})</h2>'
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
