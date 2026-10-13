#!/usr/bin/env python3
"""Combine gene-level TADA evidence from the study datasets."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


META = ["gene_id", "Hugo_Symbol", "LOEUF", "pLI"]
COHORTS = ("SSC_KOR_STAR", "SPARK")


def bayesian_fdr(bf, pi0=0.95):
    bf = np.asarray(bf, dtype=float)
    if not np.isfinite(bf).all() or (bf <= 0).any():
        raise ValueError("Combined BF values must be finite and positive")
    order = np.argsort(-bf)
    posterior_null = 1 - ((1 - pi0) * bf[order] / (pi0 + (1 - pi0) * bf[order]))
    ranked = np.cumsum(posterior_null) / np.arange(1, len(bf) + 1)
    output = np.empty_like(ranked)
    output[order] = ranked
    return output


def q2p(fdr, gene_ids):
    series = pd.Series(np.asarray(fdr, dtype=float), index=pd.Index(gene_ids, dtype=str))
    ordered = series.sort_values()
    scaled = ordered * np.arange(1, len(ordered) + 1) / len(ordered) / ordered.max()
    return scaled.reindex(series.index).to_numpy()


def validate_input(frame, label):
    required = set(META + ["bf.dn.em", "bf.in.em", "bf.em", "fdr", "pval"])
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(sorted(missing))}")
    if frame["gene_id"].isna().any() or frame["gene_id"].astype(str).eq("").any():
        raise ValueError(f"{label} contains missing gene IDs")
    if frame["gene_id"].duplicated().any():
        raise ValueError(f"{label} contains duplicate gene IDs")
    for column in ("bf.dn.em", "bf.in.em", "bf.em"):
        values = pd.to_numeric(frame[column], errors="raise")
        if not np.isfinite(values).all() or (values <= 0).any():
            raise ValueError(f"{label} column {column} must be finite and positive")


def combine_tables(ssc, spark):
    validate_input(ssc, COHORTS[0])
    validate_input(spark, COHORTS[1])
    if set(ssc["gene_id"].astype(str)) != set(spark["gene_id"].astype(str)):
        raise ValueError("Cohort gene sets do not match")
    merged = ssc.merge(spark, on=META, how="inner", suffixes=(f"_{COHORTS[0]}", f"_{COHORTS[1]}"))
    if len(merged) != len(ssc):
        raise ValueError("Cohort metadata do not match for shared gene IDs")
    for component in ("dn", "in"):
        columns = [f"bf.{component}.em_{cohort}" for cohort in COHORTS]
        merged[f"bf.{component}.em_cohortBF"] = merged[columns].clip(lower=1).prod(axis=1)
    merged["bf_floor_separate"] = merged["bf.dn.em_cohortBF"] * merged["bf.in.em_cohortBF"]
    merged["bf"] = merged["bf_floor_separate"]
    merged["fdr"] = bayesian_fdr(merged["bf"])
    merged["pval"] = q2p(merged["fdr"], merged["gene_id"])
    return merged


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ssc", type=Path, required=True)
    parser.add_argument("--spark", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    output = combine_tables(pd.read_csv(args.ssc, sep="\t"), pd.read_csv(args.spark, sep="\t"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
