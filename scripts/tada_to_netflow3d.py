#!/usr/bin/env python3
"""Convert gene-level TADA Bayes factors into NetFlow3D heat scores."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tada-results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--gene-column", default="gene_id")
    parser.add_argument("--bf-column", default="bf")
    parser.add_argument(
        "--transform",
        choices=("log", "log10", "none"),
        default="log",
        help="Transformation applied to the Bayes factor (default: log).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    table = pd.read_csv(args.tada_results, sep="\t")
    required = {args.gene_column, args.bf_column}
    missing = required.difference(table.columns)
    if missing:
        raise SystemExit(f"Missing TADA column(s): {', '.join(sorted(missing))}")

    bf = pd.to_numeric(table[args.bf_column], errors="coerce")
    if bf.isna().any() or not np.isfinite(bf).all() or (bf <= 0).any():
        raise SystemExit("All Bayes factors must be finite and greater than zero.")

    if args.transform == "log":
        score = np.log(bf)
    elif args.transform == "log10":
        score = np.log10(bf)
    else:
        score = bf

    output = pd.DataFrame({"Gene": table[args.gene_column], "Score": score})
    if output["Gene"].isna().any() or output["Gene"].duplicated().any():
        raise SystemExit("Gene identifiers must be non-missing and unique.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
