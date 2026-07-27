#!/usr/bin/env python
"""Summarise a CD-HIT ``.clstr`` file into a per-cluster CSV table.

For every cluster the script reports the identifier of the cluster
representative, the number of member sequences, and the mean percent identity
of the non-representative members to that representative.

Output schema (CSV)::

    Representative,Count,Identity,Identity_n

``Identity``   mean of the percent identities reported by CD-HIT for the
               non-representative members of the cluster. It is ``NA`` for
               singleton clusters: with a single member there is no sequence to
               compare against the representative, so the mean is undefined
               (it is *not* 100%).
``Identity_n`` number of members contributing to that mean, i.e. ``Count - 1``
               (0 for singletons).

Both CD-HIT identity notations are accepted:

* global identity (``-G 1``, the pipeline default): ``... at 95.83%``
* local alignment (``-G 0``):                       ``... at 1:120:1:118/95.83%``
* nucleotide input, with strand:                    ``... at 1:118:1:118/-/99.15%``
"""

import argparse
import re
import sys

# Matches the trailing identity field of a CD-HIT member line. Everything that
# may precede the percentage (alignment coordinates in local-alignment mode,
# strand for nucleotide input) is consumed by the optional prefix group.
IDENTITY_RE = re.compile(r"at\s+(?:\S*/)?([\d.]+)%\s*$")


def parse_clstr(handle):
    """Yield ``(representative, count, identity, identity_n)`` per cluster.

    ``identity`` is ``None`` for singleton clusters.
    """
    clusterrep = None
    clustermembers = 0
    clusteridentitysum = 0.0
    seen_header = False

    def finish():
        if clustermembers == 0:
            return None
        if clusterrep is None:
            raise ValueError("cluster with {0} member(s) has no representative ('*')".format(clustermembers))
        if clustermembers > 1:
            identity = clusteridentitysum / (clustermembers - 1)
        else:
            identity = None
        return (clusterrep, clustermembers, identity, clustermembers - 1)

    for lineno, line in enumerate(handle, start=1):
        line = line.rstrip("\n")
        if not line.strip():
            continue
        if line.startswith(">Cluster"):
            record = finish()
            if record is not None:
                yield record
            seen_header = True
            clusterrep = None
            clustermembers = 0
            clusteridentitysum = 0.0
            continue

        if not seen_header:
            raise ValueError(
                "line {0}: member line found before any '>Cluster' header: {1!r}".format(lineno, line)
            )

        if ">" not in line:
            raise ValueError(
                "line {0} is neither a cluster header nor a member line: {1!r}".format(lineno, line)
            )

        readname = line.split(">", 1)[1].split("...")[0]
        clustermembers += 1
        if line.rstrip().endswith("*"):
            if clusterrep is not None:
                raise ValueError(
                    "line {0}: two representatives found in the same cluster".format(lineno)
                )
            clusterrep = readname
        else:
            match = IDENTITY_RE.search(line)
            if match is None:
                raise ValueError(
                    "line {0}: cannot parse the identity field of a CD-HIT member line: "
                    "{1!r}".format(lineno, line)
                )
            clusteridentitysum += float(match.group(1))

    record = finish()
    if record is not None:
        yield record


def format_identity(identity):
    return "NA" if identity is None else str(identity)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Summarise a CD-HIT .clstr file into a per-cluster CSV table."
    )
    parser.add_argument("clstr", help="CD-HIT cluster file (*.clstr)")
    parser.add_argument("summary", help="output CSV file")
    args = parser.parse_args(argv)

    with open(args.clstr) as filein, open(args.summary, "w") as fileout:
        fileout.write("Representative,Count,Identity,Identity_n\n")
        for rep, count, identity, identity_n in parse_clstr(filein):
            fileout.write(
                "{0},{1},{2},{3}\n".format(rep, count, format_identity(identity), identity_n)
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
