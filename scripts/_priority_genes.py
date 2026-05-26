#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Script:  _priority_genes.py
# Purpose: Clinically prioritized AMR and virulence gene catalog plus a
#          classifier used by the MGE alert system. Genes are matched by
#          case-insensitive regular expressions so common naming variants
#          (OXA-48, blaOXA-48, OXA_48) are all recognized.
# Used by: mge_alerts.py
# Author:  Aluminion lab
# Date:    2026-05-20
# =============================================================================

"""Priority gene catalog for MGE alert classification.

Categories
----------
HIGH priority (clinical sentinels that almost always warrant an alert):
  - carbapenemase : OXA-48-like, KPC, NDM, VIM, IMP, GES (carbapenemase variants)
  - colistin      : mcr family
  - hv_klebsiella : hypervirulence markers of Klebsiella pneumoniae
                    (rmpA / rmpA2, aerobactin and yersiniabactin operons,
                    salmochelin)

MEDIUM priority (ESBL and other AMR carried by MGEs — relevant but less
acute than carbapenemase / colistin / hypervirulence):
  - esbl          : CTX-M, SHV-ESBL variants, TEM-ESBL variants

Any AMR or virulence gene that is not listed here counts as medium priority
when paired with a recurrent MGE; new (non-recurrent) MGEs only trigger an
alert when they carry at least one HIGH-priority gene (NEW_PRIORITY category).
"""

import re
from typing import Iterable

# -----------------------------------------------------------------------------
# Catalog: category -> list of regexes (case-insensitive, anchored to a word
# boundary on the left to avoid matching e.g. "ant(3'')-IIa" when looking for
# "TEM-1"). Patterns intentionally accept an optional "bla" prefix and either
# hyphen or underscore as separator before the variant number.
# -----------------------------------------------------------------------------
HIGH_PRIORITY_PATTERNS: dict[str, list[str]] = {
    "carbapenemase": [
        r"\b(bla)?OXA[-_]?(48|181|232|244|204|162|23|24|40|58|143|72)\b",
        r"\b(bla)?KPC[-_]?\d+\b",
        r"\b(bla)?NDM[-_]?\d+\b",
        r"\b(bla)?VIM[-_]?\d+\b",
        r"\b(bla)?IMP[-_]?\d+\b",
        # GES variants 5, 6, 11, 14, 16, 18, 20, 21, 24 are carbapenemase
        r"\b(bla)?GES[-_]?(5|6|11|14|16|18|20|21|24)\b",
    ],
    "colistin": [
        r"\bmcr[-_]?\d+(\.\d+)?\b",
    ],
    "hv_klebsiella": [
        # Hypermucoid regulators
        r"\brmpA2?\b",
        # Aerobactin operon (any of iucA/B/C/D, iutA)
        r"\biuc[ABCD]\b",
        r"\biutA\b",
        # Salmochelin operon
        r"\biro[BCDN]\b",
        # Yersiniabactin high-pathogenicity island marker
        r"\bybtS\b",
    ],
}

MEDIUM_PRIORITY_PATTERNS: dict[str, list[str]] = {
    "esbl": [
        r"\b(bla)?CTX[-_]?M[-_]?\d+\b",
        # SHV ESBL variants (the common ones; non-ESBL SHV-1/11/etc. excluded)
        r"\b(bla)?SHV[-_]?(2|2a|5|7|8|12|18|27|28|30|31|55|106)\b",
        # TEM ESBL variants
        r"\b(bla)?TEM[-_]?(3|10|12|26|52|158)\b",
    ],
}


def _compile(patterns: dict[str, list[str]]) -> dict[str, list[re.Pattern]]:
    """Compile all patterns once at import time."""
    return {
        category: [re.compile(p, re.IGNORECASE) for p in plist]
        for category, plist in patterns.items()
    }


_HIGH = _compile(HIGH_PRIORITY_PATTERNS)
_MEDIUM = _compile(MEDIUM_PRIORITY_PATTERNS)


def classify_priority(genes: Iterable[str]) -> dict:
    """Classify a list of gene names against the priority catalog.

    Parameters
    ----------
    genes : iterable of str
        Gene identifiers as emitted by Copla (AbR column), ABRicate (gene
        column) or VFDB. Empty strings, ``None`` and the literal ``"-"``
        placeholder are ignored.

    Returns
    -------
    dict
        ``{"priority": "high" | "medium" | "none",
           "matched": [genes that matched any pattern],
           "categories": [unique category names that fired]}``

        ``priority == "high"`` if any HIGH category fires.
        ``priority == "medium"`` if at least one gene is present but only
        MEDIUM (or no) categories fire — i.e., AMR/VIR carriage exists but
        no sentinel hit.
        ``priority == "none"`` if no genes are supplied (or only blanks).
    """
    matched: list[str] = []
    categories: list[str] = []
    any_gene_seen = False

    for gene in genes:
        if not gene or gene == "-":
            continue
        any_gene_seen = True
        for category, patterns in _HIGH.items():
            if any(p.search(gene) for p in patterns):
                matched.append(gene)
                if category not in categories:
                    categories.append(category)
        for category, patterns in _MEDIUM.items():
            if any(p.search(gene) for p in patterns):
                matched.append(gene)
                if category not in categories:
                    categories.append(category)

    if any(c in categories for c in HIGH_PRIORITY_PATTERNS):
        priority = "high"
    elif any_gene_seen:
        priority = "medium"
    else:
        priority = "none"

    return {"priority": priority, "matched": matched, "categories": categories}


def split_amr_field(amr_string: str) -> list[str]:
    """Split a Copla ``AbR`` cell or ABRicate gene list into individual genes.

    Handles both semicolon-separated (Copla) and comma-separated inputs and
    returns a list with empty strings and the ``-`` placeholder removed.
    """
    if not amr_string or amr_string == "-":
        return []
    # Copla uses ';', some downstream tables use ',' — accept both.
    raw = re.split(r"[;,]", amr_string)
    return [g.strip() for g in raw if g.strip() and g.strip() != "-"]
