#!/usr/bin/env python3
"""Placeholder comparison plot.

Reads one search-result TSV per tool and emits a bar chart of total hit counts.
The real comparison logic (sensitivity vs rank, etc.) will replace this later.
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def tool_name(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def hit_count(path: str) -> int:
    with open(path) as fh:
        return sum(1 for line in fh if line.strip())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("results", nargs="*")
    args = ap.parse_args()

    names  = [tool_name(p) for p in args.results]
    counts = [hit_count(p) for p in args.results]

    fig, ax = plt.subplots(figsize=(6, 4))
    if names:
        ax.bar(names, counts)
        ax.set_ylabel("hits")
    else:
        ax.text(0.5, 0.5, "no results yet", ha="center", va="center")
        ax.set_axis_off()
    ax.set_title("FESNov search results (placeholder)")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)


if __name__ == "__main__":
    main()
