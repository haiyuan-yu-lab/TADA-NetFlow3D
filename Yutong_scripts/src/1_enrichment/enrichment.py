import argparse
from variants import VCFHandler, VCFBEDOverlapper
from utils import dir_exists, ensure_parent_dir, load_enhancer_dict

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
    vcf_handler.pass_filters(
        output_path=case_only_filtered_vcf_path,
        vcf_path=case_only_vcf_path,
        include=[('gnomADg_v3_AF', '==', ".")],
        include_any=[
            [("TYPE","=","SNV"),  ("CADD_PHRED",">=",12)],
            [("TYPE","=","INDEL"), ("phyloP100way",">=",2)],
        ],
        threads=4,
        bgzip_index=True,
        keep_intermediate=False)
    
    vcf_handler.pass_filters(
        output_path=control_only_filtered_vcf_path,
        vcf_path=control_only_vcf_path,
        include=[('gnomADg_v3_AF', '==', ".")],
        include_any=[
            [("TYPE","=","SNV"),  ("CADD_PHRED",">=",12)],
            [("TYPE","=","INDEL"), ("phyloP100way",">=",2)],
        ],
        threads=4,
        bgzip_index=True,
        keep_intermediate=False)

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
                                       write_overlap_vcf=False)
    vcf_overlapper.map_variants_to_bed(vcf_path=control_only_filtered_vcf_path,
                                       names=args.names,
                                       file="control_only.txt", 
                                       write_overlap_vcf=False)

if __name__ == "__main__":
    args = args_parser()
    main(args)