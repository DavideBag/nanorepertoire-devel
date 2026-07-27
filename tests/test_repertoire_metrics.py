"""Unit tests for ``bin/repertoire_metrics.py``.

Run with ``pytest tests/test_repertoire_metrics.py`` from the repository root.
"""

import importlib.util
import math
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_metrics():
    spec = importlib.util.spec_from_file_location(
        "repertoire_metrics", REPO_ROOT / "bin" / "repertoire_metrics.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


rm = _load_metrics()

UNIFORM = [10] * 8
DOMINATED = [10000] + [1] * 9


def test_uniform_repertoire_has_zero_clonality():
    assert rm.clonality(UNIFORM) == pytest.approx(0.0)
    assert rm.shannon_norm(UNIFORM) == pytest.approx(1.0)


def test_uniform_shannon_equals_log_richness():
    assert rm.shannon_entropy(UNIFORM) == pytest.approx(math.log(8))


def test_dominated_repertoire_approaches_clonality_one():
    assert rm.clonality(DOMINATED) > 0.95
    assert rm.top_fraction(DOMINATED, 1) > 0.99
    assert rm.d50(DOMINATED) == 1


def test_single_clone_leaves_evenness_undefined():
    # ln(S) is 0, so Pielou evenness and clonality cannot be computed.
    assert rm.shannon_norm([500]) is None
    assert rm.clonality([500]) is None
    # The entropy itself is still defined, and is 0.
    assert rm.shannon_entropy([500]) == pytest.approx(0.0)


def test_empty_repertoire_returns_na_everywhere():
    m = rm.abundance_metrics([])
    assert m["Total_sequences"] == 0
    assert m["Clusters"] == 0
    for key in ("Shannon_abundance", "Shannon_norm", "Clonality", "Simpson",
                "Gini", "D50", "Top1_pct", "Top10_pct", "Clusters_per_1k_seqs"):
        assert m[key] is None, key


def test_zero_counts_are_ignored():
    assert rm.abundance_metrics([10, 10, 0])["Clusters"] == 2


def test_gini_simpson_bounds():
    assert rm.gini_simpson(UNIFORM) == pytest.approx(1 - 8 * (1 / 8) ** 2)
    assert rm.gini_simpson([100]) == pytest.approx(0.0)


def test_gini_coefficient():
    assert rm.gini(UNIFORM) == pytest.approx(0.0)
    assert rm.gini([100]) == pytest.approx(0.0)
    assert rm.gini(DOMINATED) > 0.85


def test_d50_on_a_known_vector():
    # 50 out of 100 sequences are reached by the two largest clones.
    assert rm.d50([30, 20, 20, 20, 10]) == 2


def test_top_fractions():
    counts = [50, 30, 20]
    assert rm.top_fraction(counts, 1) == pytest.approx(0.5)
    assert rm.top_fraction(counts, 10) == pytest.approx(1.0)


def test_depth_normalised_richness():
    m = rm.abundance_metrics([10] * 5)
    assert m["Total_sequences"] == 50
    assert m["Clusters_per_1k_seqs"] == pytest.approx(100.0)


def test_metrics_are_scale_invariant():
    small = rm.abundance_metrics([1, 2, 3, 4])
    large = rm.abundance_metrics([100, 200, 300, 400])
    for key in ("Shannon_abundance", "Clonality", "Simpson", "Gini", "Top1_pct"):
        assert small[key] == pytest.approx(large[key]), key
