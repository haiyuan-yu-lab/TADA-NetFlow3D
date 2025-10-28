import os
import argparse
from variants import VCFHandler, VCFBEDOverlapper, StatisticalTests
from plotting import plot_logOR_horizontal
from utils import dir_exists, ensure_parent_dir, load_enhancer_dict, run_command

def args_parser():
    parser = argparse.ArgumentParser(description="Enhancer variants enrichment analysis.")
    parser.add_argument("--vcf", type=str, required=True, help="Path to the input VCF file.")
    parser.add_argument("--output-dir", type=str, required=True, help="Directory to the output files.")
    parser.add_argument("--file-prefix", type=str, required=True, help="Prefix for the output files.")
    parser.add_argument("--config", type=str, required=True, help="Path to the config file.")
    parser.add_argument("--enhancer-dir", type=str, default="data/enhancer/", 
                        help="Directory containing enhancer BED files.")
    parser.add_argument("--names", type=str, nargs='+', 
                        help="List of enhancer names to analyze. \
                            Need to match labels in the config file.")
    parser.add_argument("--reference", type=str, 
                        default="/fs/cbsuhy01/storage/yz2676/data/DetailedProCapNet/wdir/reference/hg38.22XY.fasta",
                        help="Path to the reference genome file.")
    parser.add_argument("--cadd-min", type=int, default=None,
                        help="Minimum CADD_PHRED for SNVs. If omitted, CADD filter is skipped.")
    parser.add_argument("--phylop-min", type=int, default=None,
                        help="Minimum phyloP100way for INDELs. If omitted, phyloP filter is skipped.")
    parser.add_argument("--gnomad-null-only", action="store_true",
                        help="If set, require gnomADg_v3_AF == '.'")
    parser.add_argument("--gnomad-af-max", type=float, default=None,
                        help="If set, require gnomADg_v3_AF <= this value. Ignored if --gnomad-null-only is used.")
    return parser.parse_args()

def case_control_separation(vcf_handler: VCFHandler, case_vcf: str, control_vcf: str, 
                            case_outpth: str, control_outpth: str):
    """
    Separate variants into case-only and control-only VCF files.
    """
    vcf_handler.pass_filters(
        output_path=case_vcf,
        vcf_path=vcf_handler.get_vcf(),
        include=[('STATUS', '==', "case")],
        threads=4,
        bgzip_index=True,
        keep_intermediate=False)
    vcf_handler.pass_filters(
        output_path=control_vcf,
        vcf_path=vcf_handler.get_vcf(),
        include=[('STATUS', '==', "control")],
        threads=4,
        bgzip_index=True,
        keep_intermediate=False)

    vcf_handler.isec_unique(
        a_vcf_path=case_vcf,
        b_vcf_path=control_vcf,
        output_path=case_outpth
    )
    vcf_handler.isec_unique(
        a_vcf_path=control_vcf,
        b_vcf_path=case_vcf,
        output_path=control_outpth
    )

def build_filters():
        base_include = []
        include_any_groups = []

        # gnomAD logic
        if args.gnomad_null_only:
            base_include.append(('gnomADg_v3_AF', '==', "."))
        elif args.gnomad_af_max is not None:
            base_include.append(('gnomADg_v3_AF', '<=', args.gnomad_af_max))

        if args.cadd_min is not None:
            include_any_groups.append([
                ("TYPE", "=", "SNV"),
                ("CADD_PHRED", ">=", args.cadd_min)
            ])

        if args.phylop_min is not None:
            include_any_groups.append([
                ("TYPE", "=", "INDEL"),
                ("phyloP100way", ">=", args.phylop_min)
            ])

        return base_include, include_any_groups

def clear_intermediate_files(*file_paths: str):
    """
    Remove intermediates (file + possible index files) independently.
    """
    exts = ("", ".csi", ".tbi")
    for base in file_paths:
        for ext in exts:
            try:
                os.remove(f"{base}{ext}")
            except FileNotFoundError:
                pass

def collect_variant_counts(vcf_handler: VCFHandler, args: argparse.Namespace,
                           total_case_path: str, total_control_path: str) -> dict:
    """
    Collect variant counts for plotting.
    Returns a flat dictionary or a nested dictionary:
          A) Flat dict:
              { "EnhancerType": (logOR, SE, pval), ... }
          B) Nested dict (grouped):
              { "EnhancerType": {
                    "combined": (logOR, SE, pval),
                    "shared":   (logOR, SE, pval),
                    "unique":   (logOR, SE, pval)
                 }, ... }
    Parameters
    ----------
    vcf_handler : VCFHandler
        An instance of the VCFHandler class.
    args : argparse.Namespace
        Parsed command-line arguments.
    total_case_path : str
        Path to the total case VCF file.
    total_control_path : str
        Path to the total control VCF file.
    """
    var_stats = StatisticalTests()
    results = {}

    # Precompute totals once if they’re the same for every enhancer
    total_case = int(vcf_handler.bcftools_count(total_case_path))
    total_control = int(vcf_handler.bcftools_count(total_control_path))

    for enh in args.names:
        # Construct contingency table
        case_in_enh = int(run_command(
            ["wc", "-l", f"{args.output_dir}/{enh}/case_only.txt"], capture=True
        ).split()[0])
        control_in_enh = int(run_command(
            ["wc", "-l", f"{args.output_dir}/{enh}/control_only.txt"], capture=True
        ).split()[0])
        odds_ratio, log_odds_ratio, SE_log_OR, p_value = var_stats.fisher_exact_onesided(
            case_in_enh, total_case, control_in_enh, total_control
        )
        print(case_in_enh, total_case, control_in_enh, total_control)
        results[enh] = (log_odds_ratio, SE_log_OR, p_value)
        print(f"Enhancer: {enh}, logOR: {log_odds_ratio}, SE: {SE_log_OR}, p-value: {p_value}")
    return results


def main(args):
    vcf_handler = VCFHandler()
    vcf_handler.ensure_index(args.vcf)
    vcf_handler.set_vcf(args.vcf)
    dir_exists(args.output_dir)

    # Separate to case/control-only variants
    case_vcf_path = f"{args.output_dir}/{args.file_prefix}.case.vcf.gz"
    control_vcf_path = f"{args.output_dir}/{args.file_prefix}.control.vcf.gz"
    case_only_vcf_path = f"{args.output_dir}/{args.file_prefix}.case_only.vcf.gz"
    control_only_vcf_path = f"{args.output_dir}/{args.file_prefix}.control_only.vcf.gz"
    case_only_filtered_vcf_path = f"{args.output_dir}/{args.file_prefix}.case_only.filtered.vcf.gz"
    control_only_filtered_vcf_path = f"{args.output_dir}/{args.file_prefix}.control_only.filtered.vcf.gz"

    case_control_separation(
        vcf_handler=vcf_handler,
        case_vcf=case_vcf_path,
        control_vcf=control_vcf_path,
        case_outpth=case_only_vcf_path,
        control_outpth=control_only_vcf_path
    )
    print("Finished separating case/control-only variants.")

    # Filter variants based on provided criteria
    base_inc, inc_any = build_filters()
    print(f"Applied filters - Base include: {base_inc}, Include any groups: {inc_any}")
    vcf_handler.pass_filters(
        output_path=case_only_filtered_vcf_path,
        vcf_path=case_only_vcf_path,
        include=tuple(base_inc),
        include_any=tuple(inc_any),
        threads=4,
        bgzip_index=True,
        keep_intermediate=False
    )

    vcf_handler.pass_filters(
        output_path=control_only_filtered_vcf_path,
        vcf_path=control_only_vcf_path,
        include=tuple(base_inc),
        include_any=tuple(inc_any),
        threads=4,
        bgzip_index=True,
        keep_intermediate=False
    )

    # Count final variants
    case_count = vcf_handler.bcftools_count(case_only_filtered_vcf_path)
    control_count = vcf_handler.bcftools_count(control_only_filtered_vcf_path)
    print(f"Final case-only filtered variants: {case_count}")
    print(f"Final control-only filtered variants: {control_count}")

    # Overlap with enhancers
    vcf_overlapper = VCFBEDOverlapper(config=load_enhancer_dict(args.config),
                                      output_dir=args.output_dir,
                                      enhancer_dir=args.enhancer_dir,
                                      ref=args.reference)
    vcf_overlapper.map_variants_to_bed(vcf_path=case_only_filtered_vcf_path,
                                       names=args.names,
                                       file="case_only.txt", 
                                       write_overlap_vcf=True)
    vcf_overlapper.map_variants_to_bed(vcf_path=control_only_filtered_vcf_path,
                                       names=args.names,
                                       file="control_only.txt", 
                                       write_overlap_vcf=True)
    print("Finished overlapping variants with enhancers.")

    # Statistical analysis
    variant_counts = collect_variant_counts(
        vcf_handler=vcf_handler,
        args=args,
        total_case_path=case_only_filtered_vcf_path,
        total_control_path=control_only_filtered_vcf_path
    )

    # Plotting
    plot_logOR_horizontal(
        data=variant_counts,
        se_scale=1,
        annotate="stars",
        add_values=True,
        x_label="Log Odds Ratio",
        title="Enhancer Variant Enrichment",
        fig_size=(8, len(variant_counts)*0.5),
        output=f"{args.output_dir}/{args.file_prefix}_enhancer_variant_enrichment.png"
    )

    clear_intermediate_files(
        case_vcf_path,
        control_vcf_path,
        case_only_vcf_path,
        control_only_vcf_path,
        case_only_filtered_vcf_path,
        control_only_filtered_vcf_path
    )

if __name__ == "__main__":
    args = args_parser()
    main(args)