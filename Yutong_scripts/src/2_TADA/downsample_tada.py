"""
Generate downsampled coding-only TADA inputs for Fig3.

This script does one thing: binomial-thin Fu coding counts to the sample-size
scale of a target noncoding cohort, then write replicate TADA coding inputs.
The shell wrapper runs TADA on those inputs; compile_fig3_tada.py summarizes
the TADA outputs.
"""

import argparse
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd


DEFAULT_DN_CASE_COLS = ("dn.ptv", "dn.misb", "dn.misa")
DEFAULT_DN_CONTROL_COLS = ("dn.ptv.sib", "dn.misb.sib", "dn.misa.sib")
DEFAULT_IN_COUNT_COLS = (
    "in.t.ptv",
    "in.t.misb",
    "in.t.misa",
    "in.u.ptv",
    "in.u.misb",
    "in.u.misa",
)


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate downsampled coding-only TADA input replicates.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--coding-file",
        type=Path,
        required=True,
        help="Clean coding TADA input table to downsample.",
    )
    parser.add_argument(
        "--sample-file",
        type=Path,
        required=True,
        help=("Sample table containing both coding and target noncoding sample "
            "sizes for one cohort."),
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Existing directory where Downsampled_repXXXX_TADA_input.txt files are written.",
    )
    parser.add_argument(
        "--run-sample-file",
        type=Path,
        required=True,
        help="Path for the TADA sample file used by all downsampled replicates.",
    )
    parser.add_argument(
        "--n-rep",
        type=int,
        default=1000,
        help="Total number of deterministic downsampling replicates.",
    )
    parser.add_argument(
        "--rep-start",
        type=int,
        default=0,
        help="First replicate index to write, 0-based and inclusive.",
    )
    parser.add_argument(
        "--rep-end",
        type=int,
        help="One-past-last replicate index to write. Defaults to --n-rep.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2026,
        help="Seed for the NumPy RNG used by binomial thinning.",
    )
    parser.add_argument(
        "--dn-case-count-cols",
        nargs="+",
        default=list(DEFAULT_DN_CASE_COLS),
        help="De novo case coding count columns to downsample.",
    )
    parser.add_argument(
        "--dn-control-count-cols",
        nargs="+",
        default=list(DEFAULT_DN_CONTROL_COLS),
        help="De novo sibling/control coding count columns to downsample.",
    )
    parser.add_argument(
        "--in-case-count-cols",
        nargs="+",
        default=list(DEFAULT_IN_COUNT_COLS),
        help="Inherited transmitted/untransmitted coding count columns to downsample.",
    )
    return parser.parse_args(argv)


def read_sample(sample_file: Path) -> pd.DataFrame:
    """
    Read the cohort sample-size table.
    The sample table must contain the following variables:
        n_dn_case_noncoding, n_dn_control_noncoding, n_in_case_noncoding,
        n_dn_case_coding, n_dn_control_coding, n_in_case_coding
    """
    sample = pd.read_csv(sample_file, sep="\t")
    sample["value"] = pd.to_numeric(sample["value"], errors="raise")
    return sample


def compute_probabilities(sample: pd.DataFrame) -> dict[str, float]:
    """
    Compute binomial thinning probabilities from one sample table.

    Probabilities are computed as the ratio of noncoding to coding sample sizes
    for de novo case, de novo control, and inherited case counts.
    """
    values = sample.set_index("variable")["value"]
    required = (
        "n_dn_case_noncoding",
        "n_dn_control_noncoding",
        "n_in_case_noncoding",
        "n_dn_case_coding",
        "n_dn_control_coding",
        "n_in_case_coding",
    )
    missing = [key for key in required if key not in values.index]
    if missing:
        raise ValueError(f"Sample file is missing required variables: {', '.join(missing)}")

    probabilities = {
        "p_dn_case": float(values["n_dn_case_noncoding"] / values["n_dn_case_coding"]),
        "p_dn_control": float(values["n_dn_control_noncoding"] / values["n_dn_control_coding"]),
        "p_in_case": float(values["n_in_case_noncoding"] / values["n_in_case_coding"]),
    }
    return probabilities


def write_run_sample(sample: pd.DataFrame, run_sample_file: Path) -> Path:
    """
    Write the sample table TADA should use for downsampled coding inputs.
    """
    run_sample = sample.copy()
    values = run_sample.set_index("variable")["value"]
    copy_pairs = (
        ("n_dn_case_coding", "n_dn_case_noncoding"),
        ("n_dn_control_coding", "n_dn_control_noncoding"),
        ("n_in_case_coding", "n_in_case_noncoding"),
        ("n_in_control_coding", "n_in_control_noncoding"),
    )

    for coding_key, noncoding_key in copy_pairs:
        run_sample.loc[run_sample["variable"] == coding_key, "value"] = values[noncoding_key]

    run_sample.to_csv(run_sample_file, sep="\t", index=False)
    return run_sample_file


def downsample_counts(
    coding: pd.DataFrame,
    probabilities: dict[str, float],
    seed: Optional[int] = None,
    rng: Optional[np.random.Generator] = None,
    dn_case_count_cols: Union[tuple[str, ...], list[str]] = DEFAULT_DN_CASE_COLS,
    dn_control_count_cols: Union[tuple[str, ...], list[str]] = DEFAULT_DN_CONTROL_COLS,
    in_case_count_cols: Union[tuple[str, ...], list[str]] = DEFAULT_IN_COUNT_COLS,
) -> pd.DataFrame:
    """
    Return one binomial-thinned coding table.
    """
    count_cols = list(
        dict.fromkeys(
            list(dn_case_count_cols)
            + list(dn_control_count_cols)
            + list(in_case_count_cols)
        )
    )
    missing = [column for column in count_cols if column not in coding.columns]
    if missing:
        raise ValueError(f"Coding table is missing count columns: {', '.join(missing)}")

    if rng is None:
        rng = np.random.default_rng(seed)

    output = coding.copy()
    groups = (
        (dn_case_count_cols, "p_dn_case"),
        (dn_control_count_cols, "p_dn_control"),
        (in_case_count_cols, "p_in_case"),
    )
    for columns, probability_name in groups:
        for column in columns:
            counts = pd.to_numeric(output[column], errors="raise").to_numpy(dtype=int)
            output[column] = rng.binomial(counts, probabilities[probability_name])
    return output


def write_downsampled_inputs(
    coding: pd.DataFrame,
    probabilities: dict[str, float],
    input_dir: Path,
    n_rep: int,
    rep_start: int = 0,
    rep_end: Optional[int] = None,
    seed: int = 2026,
    dn_case_count_cols: Union[tuple[str, ...], list[str]] = DEFAULT_DN_CASE_COLS,
    dn_control_count_cols: Union[tuple[str, ...], list[str]] = DEFAULT_DN_CONTROL_COLS,
    in_case_count_cols: Union[tuple[str, ...], list[str]] = DEFAULT_IN_COUNT_COLS,
) -> list[Path]:
    """
    Write a deterministic replicate range while preserving global RNG order.
    """
    if rep_end is None:
        rep_end = n_rep

    rng = np.random.default_rng(seed)
    written = []
    for replicate in range(rep_end):
        downsampled = downsample_counts(
            coding,
            probabilities,
            rng=rng,
            dn_case_count_cols=dn_case_count_cols,
            dn_control_count_cols=dn_control_count_cols,
            in_case_count_cols=in_case_count_cols,
        )
        if replicate < rep_start:
            continue

        output_file = input_dir / f"Downsampled_rep{replicate:04d}_TADA_input.txt"
        downsampled.to_csv(output_file, sep="\t", index=False)
        written.append(output_file)
    return written


def main() -> None:
    args = parse_args()
    sample = read_sample(args.sample_file)
    probabilities = compute_probabilities(sample)
    write_run_sample(sample, args.run_sample_file)

    coding = pd.read_csv(args.coding_file, sep="\t")
    written = write_downsampled_inputs(
        coding,
        probabilities,
        args.input_dir,
        n_rep=args.n_rep,
        rep_start=args.rep_start,
        rep_end=args.rep_end,
        seed=args.seed,
        dn_case_count_cols=args.dn_case_count_cols,
        dn_control_count_cols=args.dn_control_count_cols,
        in_case_count_cols=args.in_case_count_cols,
    )

    print("Downsampling probabilities:")
    for name, value in probabilities.items():
        print(f"  {name}={value:.6f}")
    print(f"Wrote run sample: {args.run_sample_file}")
    print(f"Wrote {len(written)} downsampled coding inputs to {args.input_dir}")


if __name__ == "__main__":
    main()
