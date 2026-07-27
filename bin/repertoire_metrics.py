#!/usr/bin/env python3
"""Clonal abundance metrics for nanobody repertoires.

The functions in this module take a vector of clone sizes -- the ``Count``
column of ``clusterbig.csv``, i.e. the number of sequences assigned to each
CD-HIT cluster of one sample -- and return the diversity statistics used in the
AIRR literature. They are kept separate from the report generator so that they
can be unit tested independently of the HTML rendering.

Every metric is undefined for some input; ``None`` is returned in those cases
and rendered as ``NA``, rather than substituting a number that would read as a
measurement.
"""

import math

__all__ = [
    "shannon_entropy",
    "shannon_norm",
    "clonality",
    "gini_simpson",
    "gini",
    "d50",
    "top_fraction",
    "abundance_metrics",
]


def _clean(counts):
    """Keep the strictly positive clone sizes."""
    return [float(c) for c in counts if c is not None and float(c) > 0]


def shannon_entropy(counts):
    """Shannon index H' = -sum(p_i * ln p_i) over clone frequencies.

    Returns ``None`` when no sequence has been observed.
    """
    counts = _clean(counts)
    total = sum(counts)
    if total <= 0:
        return None
    return -sum((c / total) * math.log(c / total) for c in counts)


def shannon_norm(counts):
    """Pielou's evenness, H' / ln(S).

    Undefined for a single clone, where ln(S) is 0.
    """
    counts = _clean(counts)
    h = shannon_entropy(counts)
    if h is None or len(counts) < 2:
        return None
    return h / math.log(len(counts))


def clonality(counts):
    """1 - Pielou's evenness. 0 for a perfectly even repertoire."""
    evenness = shannon_norm(counts)
    return None if evenness is None else 1.0 - evenness


def gini_simpson(counts):
    """Gini-Simpson index, 1 - sum(p_i^2): probability that two sequences
    drawn at random belong to different clones."""
    counts = _clean(counts)
    total = sum(counts)
    if total <= 0:
        return None
    return 1.0 - sum((c / total) ** 2 for c in counts)


def gini(counts):
    """Gini coefficient of the clone size distribution.

    0 means every clone has the same size, and the value approaches 1 as the
    sequences concentrate into a single clone.
    """
    counts = sorted(_clean(counts))
    n = len(counts)
    total = sum(counts)
    if n == 0 or total <= 0:
        return None
    if n == 1:
        return 0.0
    weighted = sum((i + 1) * c for i, c in enumerate(counts))
    return (2.0 * weighted) / (n * total) - (n + 1.0) / n


def d50(counts):
    """Smallest number of clones accounting for at least 50% of the sequences."""
    counts = sorted(_clean(counts), reverse=True)
    total = sum(counts)
    if total <= 0:
        return None
    cumulative = 0.0
    for i, c in enumerate(counts, start=1):
        cumulative += c
        if cumulative >= total / 2.0:
            return i
    return len(counts)


def top_fraction(counts, n=1):
    """Fraction of the sequences carried by the ``n`` largest clones."""
    counts = sorted(_clean(counts), reverse=True)
    total = sum(counts)
    if total <= 0:
        return None
    return sum(counts[:n]) / total


def abundance_metrics(counts):
    """Return every clonal abundance metric for one sample.

    ``counts`` is the list of cluster sizes of that sample. Fractions are
    returned as percentages where the key ends in ``_pct``.
    """
    counts = _clean(counts)
    total = sum(counts)
    richness = len(counts)

    metrics = {
        "Total_sequences": int(total),
        "Clusters": richness,
        "Shannon_abundance": shannon_entropy(counts),
        "Shannon_norm": shannon_norm(counts),
        "Clonality": clonality(counts),
        "Simpson": gini_simpson(counts),
        "Gini": gini(counts),
        "D50": d50(counts),
        "Top1_pct": None,
        "Top10_pct": None,
        "Clusters_per_1k_seqs": None,
    }

    top1 = top_fraction(counts, 1)
    top10 = top_fraction(counts, 10)
    metrics["Top1_pct"] = None if top1 is None else top1 * 100.0
    metrics["Top10_pct"] = None if top10 is None else top10 * 100.0
    if total > 0:
        metrics["Clusters_per_1k_seqs"] = 1000.0 * richness / total

    return metrics
