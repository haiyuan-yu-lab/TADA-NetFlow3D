#!/usr/bin/env python3
"""Validate public TADA-NetFlow3D input schemas before a full run."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


SCHEMAS = {
    "info": {"gene_id", "LOEUF", "pLI"},
    "coding": {
        "gene_id", "mut.ptv", "mut.misb", "mut.misa",
        "prior.dn.ptv", "prior.dn.misb", "prior.dn.misa",
        "prior.in.ptv", "prior.in.misb", "prior.in.misa",
        "dn.ptv", "dn.ptv.sib", "dn.misb", "dn.misb.sib",
        "dn.misa", "dn.misa.sib", "in.t.ptv", "in.u.ptv",
        "in.t.misb", "in.u.misb", "in.t.misa", "in.u.misa",
    },
    "noncoding": {"gene_id", "mut.em", "dn.em", "dn.em.sib", "in.t.em", "in.u.em"},
    "variables": {"variable", "value"},
    "gene_score": {"Gene", "Score"},
    "disrupted": {"P1", "P2"},
}


def validate(path: Path, schema: str, header: bool = True, mut_col: str = "mut.em") -> None:
    table = pd.read_csv(path, sep="\t", nrows=10, header=0 if header else None)
    required = SCHEMAS[schema]
    if schema == "noncoding":
        required = (required - {"mut.em"}) | {mut_col}
    missing = required.difference(table.columns)
    if schema in {"info", "coding"} and not {"Hugo_Symbol", "gene", "gene_gencodeV33"}.intersection(table.columns):
        missing.add("gene symbol (Hugo_Symbol, gene, or gene_gencodeV33)")
    if missing:
        raise ValueError(f"{path}: missing {', '.join(sorted(missing))}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--info", type=Path)
    parser.add_argument("--coding", type=Path)
    parser.add_argument("--noncoding", type=Path)
    parser.add_argument("--mut-col", default="mut.em", help="Noncoding mutation-rate column")
    parser.add_argument("--priors", type=Path)
    parser.add_argument("--sample-sizes", type=Path)
    parser.add_argument("--gene-score", type=Path)
    parser.add_argument("--disrupted-interactions", type=Path)
    args = parser.parse_args()
    mapping = {
        "info": "info", "coding": "coding", "noncoding": "noncoding",
        "priors": "variables", "sample_sizes": "variables",
        "gene_score": "gene_score", "disrupted_interactions": "disrupted",
    }
    selected = 0
    for argument, schema in mapping.items():
        path = getattr(args, argument)
        if path:
            validate(path, schema, mut_col=args.mut_col)
            print(f"OK: {path}")
            selected += 1
    if not selected:
        parser.error("provide at least one input file")


if __name__ == "__main__":
    main()
