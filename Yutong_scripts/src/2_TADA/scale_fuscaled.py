import argparse
import pandas as pd
from pathlib import Path


BASE_DIR = Path("/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD")
DEFAULT_COHORTS = ["SSC_KOR_STAR", "SPARK"]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scale FuRaw coding counts to cohort-specific sample sizes."
    )
    parser.add_argument(
        "--input",
        default=str(BASE_DIR / "TADA_input_coding_FuRaw.txt"),
        help="FuRaw coding table to scale.",
    )
    parser.add_argument(
        "--sample-tag",
        default="InhAll_DNVAll",
        help="Sample suffix used to resolve cohort sample files.",
    )
    parser.add_argument(
        "--cohort",
        nargs="+",
        help="Cohort to scale to. Can be passed multiple times.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(BASE_DIR),
        help="Directory for scaled output tables.",
    )
    return parser.parse_args()


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"ERROR: file not found: {path}")


def load_params(sample_file: Path) -> pd.DataFrame:
    params = pd.read_csv(sample_file, 
                        sep="\t", 
                        header=0,
                        dtype={"variable": str, "value": int})
    params.set_index("variable", inplace=True)
    return params


def require_param(params: pd.DataFrame, key: str) -> int:
    if key not in params.index:
        raise ValueError(f"ERROR: missing '{key}' in sample file.")
    return params.loc[key, "value"]


def calc_scale(numerator: int, denominator: int) -> float:
    if denominator == 0:
        raise ValueError("ERROR: denominator is 0 while computing scale.")
    return numerator / denominator


def label_from_sample_file(sample_file: Path) -> str:
    label = sample_file.name
    label = label[7:-18]  # cohort
    return label


def scale_table(input_file: Path, output_path: Path, scales: dict[str, float]) -> None:
    df = pd.read_table(input_file)
    scaled_cols = df.columns[16:]

    scale_by_col = {}
    for idx, column_name in enumerate(scaled_cols, start=17):
        if column_name.startswith("dn.") and not column_name.endswith(".sib"):
            scale_by_col[column_name] = scales["dn_case"]
        elif column_name.startswith("dn.") and column_name.endswith(".sib"):
            scale_by_col[column_name] = scales["dn_control"]
        elif column_name.startswith("in.t."):
            scale_by_col[column_name] = scales["in_case"]
        elif column_name.startswith("in.u."):
            scale_by_col[column_name] = scales["in_control"]
        else:
            raise ValueError(
                f"ERROR: unexpected column header after col16: col={idx} header={column_name}"
            )

    for column_name, scale in scale_by_col.items():
        df[column_name] = pd.to_numeric(df[column_name], errors="raise").mul(scale).round().astype(int)

    df.to_csv(output_path, sep="\t", index=False)


def main():
    args = parse_arguments()

    fu_raw = Path(args.input)
    output_dir = Path(args.output_dir)
    cohorts = args.cohort if args.cohort is not None else DEFAULT_COHORTS

    source_sample = BASE_DIR / f"sample_FuRaw_{args.sample_tag}.txt"
    target_samples = [BASE_DIR / f"sample_{cohort}_{args.sample_tag}.txt" for cohort in cohorts]

    require_file(fu_raw)
    require_file(source_sample)
    output_dir.mkdir(parents=True, exist_ok=True)

    source_params = load_params(source_sample)
    source_dn_case = require_param(source_params, "n_dn_case_coding")
    source_dn_control = require_param(source_params, "n_dn_control_coding")
    source_in_case = require_param(source_params, "n_in_case_coding")
    source_in_control = require_param(source_params, "n_in_control_coding")

    for target_sample in target_samples:
        require_file(target_sample)
        target_params = load_params(target_sample)

        scales = {
            "dn_case": calc_scale(require_param(target_params, "n_dn_case_noncoding"), source_dn_case),
            "dn_control": calc_scale(
                require_param(target_params, "n_dn_control_noncoding"), source_dn_control
            ),
            "in_case": calc_scale(require_param(target_params, "n_in_case_noncoding"), source_in_case),
            "in_control": calc_scale(
                require_param(target_params, "n_in_control_noncoding"), source_in_control
            ),
        }

        label = label_from_sample_file(target_sample)
        scaled_output = output_dir / f"TADA_input_coding_FuScaled_{label}.txt"
        scale_table(fu_raw, scaled_output, scales)

        print(f"Generated {scaled_output}")
        print(f"  target sample: {target_sample}")
        print(f"  dn case scale: {scales['dn_case']:.6f}")
        print(f"  dn control scale: {scales['dn_control']:.6f}")
        print(f"  in case scale: {scales['in_case']:.6f}")
        print(f"  in control scale: {scales['in_control']:.6f}")


if __name__ == "__main__":
    main()
