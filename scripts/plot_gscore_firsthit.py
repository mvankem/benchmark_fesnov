#!/usr/bin/env python3
"""Cumulative plot of gscore for the top hit per query.

Queries with no hit are counted as gscore=0 (universe taken from qDB .lookup).
"""

import argparse
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# columns in gscore_*.tsv: query,target,evalue,bits,qstart,tstart,cigar,gscore,alntmscore


def first_hit_scores(path: str, col: int) -> dict[str, float]:
    out: dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= col:
                continue  # would be better to fail here
            q = parts[0]
            if q in out:
                continue
            v = float(parts[col])
            out[q] = 0.0 if math.isnan(v) else v
    return out


def queries_from_lookup(path: str) -> list[str]:
    qs = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                qs.append(parts[1])
    return qs


def label_for(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--queries", required=True, help="path to qDB .lookup")
    ap.add_argument("--col", type=int, default=7, help="score column index (0-based)")
    ap.add_argument("--ylabel", default="gscore of top hit")
    ap.add_argument("output")
    ap.add_argument("inputs", nargs="+")
    args = ap.parse_args()

    queries = queries_from_lookup(args.queries)
    n = len(queries)

    fig, ax = plt.subplots(figsize=(7, 4))
    for p in args.inputs:
        hits = first_hit_scores(p, args.col)
        vals = sorted((hits.get(q, 0.0) for q in queries), reverse=True)
        xs = [(i + 1) / n for i in range(n)]
        ax.plot(xs, vals, label=label_for(p))
    ax.set_xlabel("fraction of queries")
    ax.set_ylabel(args.ylabel)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.output, dpi=300)


if __name__ == "__main__":
    main()
