#!/usr/bin/env python3
"""Annotate CDR3s with nanoCDR-X (predict_cdrs) and write the per-sample CDR3 outputs.

Every input sequence appears exactly once in <prefix>_cdr3.tsv, with one status:
  unique                CDR3 called, first occurrence in the sample
  non-unique            CDR3 called, already seen in the sample
  cdr3-too-long         CDR3 called, longer than MAX_CDR3 residues (left out of the histogram)
  no-cdr3               sent to the model, no CDR3 called
  excluded-length       not sent to the model: length outside MIN_LENGTH-MAX_LENGTH
  excluded-nonstandard  not sent to the model: contains a residue the model does not know (e.g. X)
<prefix>_cdr3_summary.tsv counts the sequences in each group. The script exits
with an error if predict_cdrs fails or does not return every sequence it was given.
"""
import argparse
import csv
import subprocess
import sys
from collections import Counter

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LENGTH = 70
MAX_LENGTH = 150  # input length of the nanoCDR-X model; longer sequences make predict_cdrs fail
MAX_CDR3 = 50


def read_sequences(path):
    """Yield (identifier, sequence); headers may start with '>' or '@' (translate.py keeps the FASTQ header)."""
    identifier, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line[0] in ">@":
                if identifier is not None:
                    yield identifier, "".join(chunks).upper()
                # first word only, without commas, which would break the CSV given to predict_cdrs
                identifier, chunks = line[1:].split()[0].replace(",", "_"), []
            else:
                chunks.append(line)
    if identifier is not None:
        yield identifier, "".join(chunks).upper()


def predict(sequences, extra_args):
    """Run predict_cdrs once on all sequences and return {identifier: predicted CDR3}."""
    with open("nanocdrx_input.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["identifier", "input"])
        writer.writerows(sequences)
    cmd = ["predict_cdrs", "-i", "nanocdrx_input.csv", "-o", "nanocdrx_output.csv"] + extra_args
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stdout + result.stderr)
        sys.exit(f"predict_cdrs failed (exit code {result.returncode}) on {len(sequences)} sequences")

    with open("nanocdrx_output.csv", newline="") as handle:
        calls = {row["identifier"]: row.get("predicted_cdr3") or "" for row in csv.DictReader(handle)}
    missing = [identifier for identifier, _ in sequences if identifier not in calls]
    if missing:
        sys.exit(f"predict_cdrs returned no result for {len(missing)} of {len(sequences)} sequences, "
                 f"e.g. {', '.join(missing[:5])}")
    return calls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-p", "--prefix", required=True)
    args, extra_args = parser.parse_known_args()  # unrecognised arguments are passed on to predict_cdrs

    rows = []  # [identifier, CDR3, sequence, status], in input order
    for identifier, sequence in read_sequences(args.input):
        if not MIN_LENGTH <= len(sequence) <= MAX_LENGTH:
            status = "excluded-length"
        elif not set(sequence) <= STANDARD_AA:
            status = "excluded-nonstandard"
        else:
            status = None  # decided after the prediction
        rows.append([identifier, "NA", sequence, status])

    to_model = [(row[0], row[2]) for row in rows if row[3] is None]
    calls = predict(to_model, extra_args) if to_model else {}

    seen = set()
    lengths = Counter()
    for row in rows:
        if row[3] is not None:
            continue
        cdr3 = calls[row[0]]
        if cdr3 in ("", "nan"):
            row[3] = "no-cdr3"
            continue
        row[1] = cdr3
        if len(cdr3) > MAX_CDR3:
            row[3] = "cdr3-too-long"
        elif cdr3 in seen:
            row[3] = "non-unique"
        else:
            row[3] = "unique"
            seen.add(cdr3)
            lengths[len(cdr3)] += 1

    with open(f"{args.prefix}_cdr3.tsv", "w") as tsv, open(f"{args.prefix}_cdr3.fasta", "w") as fasta:
        tsv.write("ID\tCDR3\tsequence\tunique\n")
        for identifier, cdr3, sequence, status in rows:
            tsv.write(f"{identifier}\t{cdr3}\t{sequence}\t{status}\n")
            if status == "unique":
                fasta.write(f">{identifier}\n{cdr3}\n")

    with open(f"{args.prefix}_cdr3.hist", "w") as hist:
        hist.write("Size,Count\n")
        for size in range(MAX_CDR3 + 1):
            hist.write(f"{size},{lengths[size]}\n")

    status = Counter(row[3] for row in rows)
    summary = {
        "sample": args.prefix,
        "input_sequences": len(rows),
        "excluded_length": status["excluded-length"],
        "excluded_nonstandard": status["excluded-nonstandard"],
        "sent_to_model": len(to_model),
        "with_cdr3": status["unique"] + status["non-unique"] + status["cdr3-too-long"],
        "no_cdr3": status["no-cdr3"],
        "unique": status["unique"],
        "non_unique": status["non-unique"],
        "cdr3_too_long": status["cdr3-too-long"],
    }
    with open(f"{args.prefix}_cdr3_summary.tsv", "w") as handle:
        handle.write("\t".join(summary) + "\n")
        handle.write("\t".join(str(value) for value in summary.values()) + "\n")


if __name__ == "__main__":
    main()
