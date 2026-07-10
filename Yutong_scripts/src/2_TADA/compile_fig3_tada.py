"""
Compile Fig3 downsampled coding TADA outputs.

1. summarize coding-only downsampled TADA replicates by mean log component BF,
2. add cohort-specific coding evidence across cohorts,
3. add existing cohort-specific noncoding EM evidence,
4. write unclipped TADA summaries plus clipped NetFlow3D heat-score inputs.
"""

import argparse
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd


COMPONENT_BF_COLUMNS = (
    "bf.dn.ptv",
    "bf.in.ptv",
    "bf.dn.misb",
    "bf.in.misb",
    "bf.dn.misa",
    "bf.in.misa",
)
VARIANT_TYPES = ("ptv", "misb", "misa")
DEFAULT_COHORTS = ("SSC_KOR_STAR", "SPARK")


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compile Fig3 coding-only downsampled TADA replicates and combine "
            "them with existing noncoding EM evidence."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--downsample-dir",
        type=Path,
        required=True,
        help=(
            "Directory containing <cohort>/coding_only/TADA_results; compiled "
            "summary tables are also written here."
        ),
    )
    parser.add_argument(
        "--noncoding-file",
        type=Path,
        default=(
            "data/TADA/InhAll_DNVAll/ALL/TADA_results/"
            "TADA_NoncodingOnly_Combined_Simple_Roulette_MR_rd_control.txt"
        ),
        help="Existing cohort-specific noncoding EM summary to add to coding evidence.",
    )
    parser.add_argument(
        "--netflow3d-input-dir",
        type=Path,
        default="data/netflow3d/input",
        help="Directory for NetFlow3D input tables.",
    )
    parser.add_argument(
        "--cohorts",
        nargs="+",
        default=list(DEFAULT_COHORTS),
        help="Cohorts to compile and sum as separated coding evidence.",
    )
    parser.add_argument(
        "--n-rep",
        type=int,
        default=1000,
        help="Number of TADA replicate outputs expected per cohort.",
    )
    parser.add_argument(
        "--output-prefix",
        default="Fig3_downsampled",
        help="Prefix for compiled summary and NetFlow3D input filenames.",
    )
    return parser.parse_args(argv)


def bayesian_fdr(bayes_factors: Union[pd.Series, np.ndarray], pi0: float = 0.95) -> np.ndarray:
    """
    Calculate a Bayesian false discovery rate for each Bayes factor.

    Bayes factors are ranked from strongest to weakest evidence. The FDR for
    each rank is the cumulative mean posterior probability of the null, using
    'pi0' as the prior probability that a gene is not associated.

    Returns:
        An array of Bayesian FDR values in the same order as the input Bayes
        factors.
    """
    bf = np.asarray(bayes_factors, dtype=float)
    order = np.argsort(-bf)
    bf_sorted = bf[order]
    pi_alt = 1.0 - pi0
    posterior_alt = pi_alt * bf_sorted / (1.0 - pi_alt + pi_alt * bf_sorted)
    posterior_null = 1.0 - posterior_alt
    fdr_sorted = np.cumsum(posterior_null) / np.arange(1, len(bf_sorted) + 1)
    fdr = np.empty_like(fdr_sorted)
    fdr[order] = fdr_sorted
    return fdr


def compile_cohort(result_dir: Path, n_rep: int, cohort: str) -> pd.DataFrame:
    """
    Summarize coding Bayes factors across one cohort's TADA replicates.

    For each gene and coding evidence component, this function calculates the
    mean log Bayes factor across independently downsampled replicates generated
    by binomial thinning. Every requested replicate is expected to contain the
    same genes and components, and files are processed one at a time to limit
    memory use.

    Returns:
        A gene-level table containing available gene metadata, the cohort
        label, and the mean log Bayes factor for each coding component.
    """
    for replicate in range(n_rep):
        df = pd.read_csv(result_dir / f"Downsampled_rep{replicate:04d}_TADA_output.txt", sep="\t")
        df["gene_id"] = df["gene_id"].astype(str)
        log_bf = np.log(df.set_index("gene_id")[list(COMPONENT_BF_COLUMNS)])
        if replicate == 0:
            output = df[["gene_id", "Hugo_Symbol", "LOEUF", "pLI"]].set_index("gene_id")
            log_bf_sum = log_bf
        else:
            log_bf_sum = log_bf_sum.add(log_bf)

    mean_log_bf = (log_bf_sum / n_rep).add_prefix("mean_log_")
    output = output.join(mean_log_bf).reset_index()
    output.insert(1, "cohort", cohort)
    return output


def combine_cohort_summaries(cohort_summaries: list[pd.DataFrame]) -> pd.DataFrame:
    """
    Merge and combine two cohort-specific coding summaries by gene.

    De novo mean log Bayes factors are added directly; inherited means are
    floored at zero within each cohort before addition.

    Returns:
        A gene-level coding table containing cohort-specific component values,
        cross-cohort component totals, total coding evidence, Bayesian FDR, and
        the coding heat score.

    """
    first, second = cohort_summaries
    first_cohort = str(first["cohort"].iloc[0])
    second_cohort = str(second["cohort"].iloc[0])
    mean_columns = [f"mean_log_{component}" for component in COMPONENT_BF_COLUMNS]
    combined = first[["gene_id", "Hugo_Symbol", "LOEUF", "pLI"] + mean_columns].merge(
        second[["gene_id"] + mean_columns],
        on="gene_id",
        suffixes=(f"_{first_cohort}", f"_{second_cohort}"),
    )

    for component in COMPONENT_BF_COLUMNS:
        first_column = f"mean_log_{component}_{first_cohort}"
        second_column = f"mean_log_{component}_{second_cohort}"
        if component.startswith("bf.in."):
            combined[f"mean_log_{component}"] = (
                combined[first_column].clip(lower=0)
                + combined[second_column].clip(lower=0)
            )
        else:
            combined[f"mean_log_{component}"] = combined[first_column] + combined[second_column]

    for variant in VARIANT_TYPES:
        combined[f"log_bf.{variant}"] = (
            combined[f"mean_log_bf.dn.{variant}"]
            + combined[f"mean_log_bf.in.{variant}"]
        )

    combined["log_bf.coding"] = combined[
        [f"log_bf.{variant}" for variant in VARIANT_TYPES]
    ].sum(axis=1)
    combined["bf_coding"] = np.exp(combined["log_bf.coding"])
    combined["fdr_coding"] = bayesian_fdr(combined["bf_coding"])
    combined["heat_score_coding"] = combined["log_bf.coding"].clip(lower=0)
    return combined


def combine_coding_and_noncoding(
    coding: pd.DataFrame,
    noncoding: pd.DataFrame,
    cohorts=DEFAULT_COHORTS,
) -> pd.DataFrame:
    """
    Combine coding and noncoding evidence into the final Fig3 TADA summary.

    Coding and noncoding log Bayes factors are inner-merged and added by gene.
    The combined Bayes factor and Bayesian FDR are calculated from the unclipped
    total; clipping is applied only to the NetFlow3D heat score.

    Returns:
        A gene-level table containing coding, noncoding, and combined evidence,
        coding and combined FDR values, and a nonnegative combined heat score.
    """
    noncoding_columns = [f"bf.em_{cohort}" for cohort in cohorts]
    noncoding = noncoding[["gene_id"] + noncoding_columns].copy()
    noncoding["log_bf.em_separate"] = np.log(noncoding[noncoding_columns]).sum(axis=1)
    merged = coding.merge(noncoding[["gene_id", "log_bf.em_separate"]], on="gene_id")
    merged["log_bf"] = merged["log_bf.coding"] + merged["log_bf.em_separate"]
    merged["bf"] = np.exp(merged["log_bf"])
    merged["fdr_addEM"] = bayesian_fdr(merged["bf"])
    merged["heat_score"] = merged["log_bf"].clip(lower=0)
    return merged


def write_table(df: pd.DataFrame, path: Path) -> None:
    """
    Write a DataFrame as a tab-delimited table, creating parent directories.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_netflow_input(df: pd.DataFrame, score_col: str, path: Path) -> None:
    """
    Write a gene-score table in the two-column format expected by NetFlow3D.

    The input 'gene_id' and selected score columns are renamed to 'Gene'
    and 'Score' in the output table.
    """
    out = df[["gene_id", score_col]].rename(columns={"gene_id": "Gene", score_col: "Score"})
    write_table(out, path)


def main() -> None:
    """
    Compile cohort replicates and write Fig3 TADA and NetFlow3D input tables.

    The workflow summarizes each cohort independently, combines coding evidence
    across cohorts, adds existing noncoding evidence, and writes both detailed
    summaries and two-column NetFlow3D inputs.
    """
    args = parse_args()

    cohort_summaries = []
    for cohort in args.cohorts:
        summary = compile_cohort(result_dir=args.downsample_dir / cohort / "coding_only" / "TADA_results", 
                                 n_rep=args.n_rep, 
                                 cohort=cohort)
        cohort_summary_path = args.downsample_dir / f"{args.output_prefix}_{cohort}_coding_component_mean_log_bf.tsv"
        write_table(summary, cohort_summary_path)
        cohort_summaries.append(summary)
        print(f"Wrote cohort summary: {cohort_summary_path}")

    coding = combine_cohort_summaries(cohort_summaries)
    coding_path = args.downsample_dir / f"{args.output_prefix}_cohort_separated_coding_summary.tsv"
    write_table(coding, coding_path)
    print(f"Wrote coding summary: {coding_path}")

    noncoding = pd.read_csv(args.noncoding_file, sep="\t")
    combined = combine_coding_and_noncoding(coding, noncoding, cohorts=args.cohorts)
    combined_path = args.downsample_dir / f"{args.output_prefix}_cohort_separated_coding_noncoding_summary.tsv"
    write_table(combined, combined_path)
    print(f"Wrote coding+noncoding summary: {combined_path}")

    coding_netflow_path = args.netflow3d_input_dir / f"{args.output_prefix}_coding_NetFlow3D_Input.txt"
    combined_netflow_path = args.netflow3d_input_dir / f"{args.output_prefix}_coding_noncoding_NetFlow3D_Input.txt"
    write_netflow_input(coding, "heat_score_coding", coding_netflow_path)
    write_netflow_input(combined, "heat_score", combined_netflow_path)
    print(f"Wrote NetFlow3D coding input: {coding_netflow_path}")
    print(f"Wrote NetFlow3D coding+noncoding input: {combined_netflow_path}")


if __name__ == "__main__":
    main()
