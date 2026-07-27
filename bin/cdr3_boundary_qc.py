#!/usr/bin/env python3
"""Quality control of the CDR3 boundaries called by nanoCDR-X.

The workflow annotates the CDR3 with a deep-learning caller (nanoCDR-X). This
script compares those calls, on the very same translated sequences, with the
motif-based definition documented in ``assets/analysis_report.qmd`` and
implemented in ``bin/getcdr3.py``:

* the CDR3 **starts** immediately after the cysteine of the framework-3 anchor
  ``T.{2}Y.{1}C`` (the ``YYC`` motif);
* the CDR3 **ends** immediately before the tryptophan of the framework-4 anchor
  ``WG.G`` (the ``WGQ`` motif).

Offsets are expressed in residues, relative to that reference:

``n_offset`` = start(nanoCDR-X) − start(reference)
    positive when the caller starts further into the loop, negative when it
    reaches back into framework 3.
``n_offset`` = end(nanoCDR-X) − end(reference)
    positive when the caller extends past the reference end, i.e. when
    framework-4 residues leak into the annotated CDR3.

This script only measures; it changes no annotation.

Input is one or more ``fastaSeq.csv``-style tables carrying at least the
columns ``CDR3`` and ``sequence`` (``Sample``/``sample`` optional).
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

# Framework-3 anchor, as used by bin/getcdr3.py. The CDR3 reference start is the
# end of this match, i.e. the residue following the anchor cysteine.
FR3_ANCHOR = re.compile(r"T.{2}Y.{1}C")
# Framework-4 anchor. The CDR3 reference end is the start of this match.
FR4_ANCHOR = re.compile(r"WG.G")
# Reads in which the conserved framework-4 tryptophan is substituted still carry
# the C-terminal ``TVSS`` motif used by bin/getcdr3.py. The residue occupying the
# tryptophan position is the sixth one before ``TVSS``, so the reference end is
# one residue before that match, which is exactly what ``WG.G`` would have given.
FR4_FALLBACK = re.compile(r".{6}TVSS")

NULL_VALUES = {"", "NA", "NaN", "nan", "None"}


def reference_boundaries(sequence):
    """Return ``(start, end, anchor)`` of the motif-based CDR3, or ``None``.

    ``anchor`` records which framework-4 rule was used, ``WG.G`` or the
    ``TVSS`` fallback. ``None`` is returned when either anchor is missing or
    they are out of order, the same failure mode ``bin/getcdr3.py`` reports as
    ``no-cdr3``.
    """
    fr3 = FR3_ANCHOR.search(sequence)
    if fr3 is None:
        return None

    fr4 = FR4_ANCHOR.search(sequence, fr3.end())
    if fr4 is not None:
        end, anchor = fr4.start(), "WG.G"
    else:
        fallback = FR4_FALLBACK.search(sequence, fr3.end())
        if fallback is None:
            return None
        end, anchor = fallback.start() - 1, "TVSS"

    if end <= fr3.end():
        return None
    return fr3.end(), end, anchor


def compare_row(cdr3, sequence):
    """Compare one nanoCDR-X call with the motif-based reference.

    Returns a dict with ``n_offset``, ``c_offset`` and ``exact``, or a dict with
    a ``skipped`` reason when the comparison is not possible.
    """
    cdr3 = (cdr3 or "").strip().upper()
    sequence = (sequence or "").strip().upper()
    if cdr3 in NULL_VALUES or not sequence:
        return {"skipped": "no_call"}

    reference = reference_boundaries(sequence)
    if reference is None:
        return {"skipped": "no_reference_motif"}

    start = sequence.find(cdr3)
    if start < 0:
        return {"skipped": "call_not_found_in_sequence"}

    ref_start, ref_end, anchor = reference
    return {
        "n_offset": start - ref_start,
        "c_offset": (start + len(cdr3)) - ref_end,
        "exact": start == ref_start and start + len(cdr3) == ref_end,
        "called_length": len(cdr3),
        "reference_length": ref_end - ref_start,
        "anchor": anchor,
    }


def _median(values):
    if not values:
        return None
    values = sorted(values)
    mid = len(values) // 2
    if len(values) % 2:
        return float(values[mid])
    return (values[mid - 1] + values[mid]) / 2.0


def _quantile(values, q):
    if not values:
        return None
    values = sorted(values)
    idx = min(len(values) - 1, max(0, int(round(q * (len(values) - 1)))))
    return float(values[idx])


def summarise(rows):
    """Aggregate per-sample QC statistics from ``(sample, cdr3, sequence)`` rows."""
    per_sample = defaultdict(lambda: {
        "total": 0, "compared": 0, "exact": 0,
        "n_offsets": [], "c_offsets": [],
        "called_lengths": [], "reference_lengths": [],
        "anchors": Counter(), "skipped": Counter(),
    })

    for sample, cdr3, sequence in rows:
        acc = per_sample[sample]
        acc["total"] += 1
        result = compare_row(cdr3, sequence)
        if "skipped" in result:
            acc["skipped"][result["skipped"]] += 1
            continue
        acc["compared"] += 1
        acc["exact"] += 1 if result["exact"] else 0
        acc["n_offsets"].append(result["n_offset"])
        acc["c_offsets"].append(result["c_offset"])
        acc["called_lengths"].append(result["called_length"])
        acc["reference_lengths"].append(result["reference_length"])
        acc["anchors"][result["anchor"]] += 1

    summary = {}
    for sample, acc in per_sample.items():
        compared = acc["compared"]
        summary[sample] = {
            "sequences": acc["total"],
            "compared": compared,
            "exact_concordance_pct": (100.0 * acc["exact"] / compared) if compared else None,
            "median_n_offset": _median(acc["n_offsets"]),
            "median_c_offset": _median(acc["c_offsets"]),
            "q1_n_offset": _quantile(acc["n_offsets"], 0.25),
            "q3_n_offset": _quantile(acc["n_offsets"], 0.75),
            "q1_c_offset": _quantile(acc["c_offsets"], 0.25),
            "q3_c_offset": _quantile(acc["c_offsets"], 0.75),
            "median_called_length": _median(acc["called_lengths"]),
            "median_reference_length": _median(acc["reference_lengths"]),
            "n_offset_distribution": Counter(acc["n_offsets"]),
            "c_offset_distribution": Counter(acc["c_offsets"]),
            "reference_from_tvss_fallback": acc["anchors"].get("TVSS", 0),
            "skipped": dict(acc["skipped"]),
        }
    return summary


SUMMARY_COLUMNS = [
    "Sample", "Sequences", "Compared", "Exact_concordance_pct",
    "Median_N_offset", "Q1_N_offset", "Q3_N_offset",
    "Median_C_offset", "Q1_C_offset", "Q3_C_offset",
    "Median_called_length", "Median_reference_length", "Reference_from_TVSS_fallback",
    "Skipped_no_call", "Skipped_no_reference_motif", "Skipped_call_not_found",
]


def _fmt(value, digits=2):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def write_summary_tsv(summary, path):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(SUMMARY_COLUMNS)
        for sample in sorted(summary):
            s = summary[sample]
            writer.writerow([
                sample, s["sequences"], s["compared"],
                _fmt(s["exact_concordance_pct"]),
                _fmt(s["median_n_offset"], 1), _fmt(s["q1_n_offset"], 1), _fmt(s["q3_n_offset"], 1),
                _fmt(s["median_c_offset"], 1), _fmt(s["q1_c_offset"], 1), _fmt(s["q3_c_offset"], 1),
                _fmt(s["median_called_length"], 1), _fmt(s["median_reference_length"], 1),
                s["reference_from_tvss_fallback"],
                s["skipped"].get("no_call", 0),
                s["skipped"].get("no_reference_motif", 0),
                s["skipped"].get("call_not_found_in_sequence", 0),
            ])


def write_distribution_tsv(summary, path):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", "Terminus", "Offset", "Sequences", "Pct"])
        for sample in sorted(summary):
            s = summary[sample]
            for terminus, dist in (("N", s["n_offset_distribution"]),
                                   ("C", s["c_offset_distribution"])):
                total = sum(dist.values())
                for offset in sorted(dist):
                    count = dist[offset]
                    pct = (100.0 * count / total) if total else 0.0
                    writer.writerow([sample, terminus, offset, count, f"{pct:.2f}"])


def read_tables(paths):
    """Yield ``(sample, cdr3, sequence)`` from fastaSeq-style CSV/TSV files."""
    for path in paths:
        ext = os.path.splitext(path)[1].lower()
        delimiter = "\t" if ext in (".tsv", ".tab") else ","
        default_sample = os.path.basename(path).split(".")[0]
        with open(path, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter=delimiter):
                sample = row.get("Sample") or row.get("sample") or default_sample
                yield sample, row.get("CDR3"), row.get("sequence") or row.get("Sequence")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fastaseq", nargs="+", required=True,
                        help="One or more fastaSeq.csv-style tables with CDR3 and sequence columns")
    parser.add_argument("--summary", default="cdr3_boundary_qc.tsv",
                        help="Output TSV with the per-sample summary")
    parser.add_argument("--distribution", default="cdr3_boundary_offsets.tsv",
                        help="Output TSV with the full offset distribution")
    args = parser.parse_args(argv)

    summary = summarise(read_tables(args.fastaseq))
    write_summary_tsv(summary, args.summary)
    write_distribution_tsv(summary, args.distribution)

    for sample in sorted(summary):
        s = summary[sample]
        print(
            f"{sample}: compared {s['compared']}/{s['sequences']} sequences, "
            f"exact concordance {_fmt(s['exact_concordance_pct'])}%, "
            f"median N offset {_fmt(s['median_n_offset'], 1)}, "
            f"median C offset {_fmt(s['median_c_offset'], 1)}",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
