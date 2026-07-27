"""Unit tests for ``bin/readcdout.py``.

Run with ``pytest tests/test_readcdout.py`` from the repository root.
"""

import csv
import importlib.util
import io
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
CLSTR = REPO_ROOT / "tests" / "data" / "test_clusters.clstr"


def _load_readcdout():
    spec = importlib.util.spec_from_file_location("readcdout", REPO_ROOT / "bin" / "readcdout.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


readcdout = _load_readcdout()


@pytest.fixture(scope="module")
def summary_rows(tmp_path_factory):
    out = tmp_path_factory.mktemp("readcdout") / "test_clusters.summary"
    readcdout.main([str(CLSTR), str(out)])
    with open(out, newline="") as handle:
        return list(csv.DictReader(handle))


def test_header_schema(summary_rows):
    assert list(summary_rows[0].keys()) == ["Representative", "Count", "Identity", "Identity_n"]


def test_all_clusters_reported(summary_rows):
    reps = [r["Representative"] for r in summary_rows]
    assert reps == ["singleton_rep", "pair_rep", "five_rep", "local_rep", "local_rev_rep"]


def test_singleton_identity_is_na(summary_rows):
    singleton = summary_rows[0]
    assert singleton["Count"] == "1"
    assert singleton["Identity"] in ("NA", "")
    assert singleton["Identity_n"] == "0"


def test_two_member_cluster_uses_the_single_member_identity(summary_rows):
    pair = summary_rows[1]
    assert pair["Count"] == "2"
    assert float(pair["Identity"]) == pytest.approx(90.44)
    assert pair["Identity_n"] == "1"


def test_five_member_cluster_is_the_mean_of_non_representatives(summary_rows):
    five = summary_rows[2]
    assert five["Count"] == "5"
    assert five["Identity_n"] == "4"
    expected = (92.65 + 93.38 + 96.32 + 97.79) / 4
    assert float(five["Identity"]) == pytest.approx(expected)


def test_local_alignment_format_is_parsed(summary_rows):
    local = summary_rows[3]
    assert float(local["Identity"]) == pytest.approx(95.83)
    assert local["Identity_n"] == "1"


def test_strand_annotated_format_is_parsed(summary_rows):
    local_rev = summary_rows[4]
    assert float(local_rev["Identity"]) == pytest.approx(99.15)


def test_no_row_carries_the_legacy_placeholder_identity(summary_rows):
    # The old implementation wrote the fraction ``1`` for singletons, mixing
    # units within the same column.
    assert all(r["Identity"] != "1" for r in summary_rows)
    assert all(r["Identity"] not in ("1.0", "1") for r in summary_rows)


def test_empty_input_writes_only_the_header(tmp_path):
    empty = tmp_path / "empty.clstr"
    empty.write_text("")
    out = tmp_path / "empty.summary"
    assert readcdout.main([str(empty), str(out)]) == 0
    assert out.read_text() == "Representative,Count,Identity,Identity_n\n"


def test_unparsable_identity_raises_with_line_number():
    broken = io.StringIO(">Cluster 0\n0\t120aa, >a... *\n1\t120aa, >b... at unknown\n")
    with pytest.raises(ValueError) as excinfo:
        list(readcdout.parse_clstr(broken))
    assert "line 3" in str(excinfo.value)


def test_member_before_header_raises():
    orphan = io.StringIO("0\t120aa, >a... *\n")
    with pytest.raises(ValueError) as excinfo:
        list(readcdout.parse_clstr(orphan))
    assert "before any" in str(excinfo.value)


def test_cli_help_is_available(capsys):
    with pytest.raises(SystemExit) as excinfo:
        readcdout.main(["--help"])
    assert excinfo.value.code == 0
    assert "clstr" in capsys.readouterr().out
