import argparse
from variants import VCFHandler
from utils import ensure_parent_dir

def args_parser():
    parser = argparse.ArgumentParser(description="Enhancer variants enrichment analysis.")
    parser.add_argument("--vcf", type=str, required=True, help="Path to the input VCF file.")
    parser.add_argument("--output", type=str, required=True, help="Path to the output file.")
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
    ensure_parent_dir(args.output)

    # Separate to case/control-only variants
    case_vcf_path = args.output.replace(".vcf.gz", ".case.vcf.gz")
    control_vcf_path = args.output.replace(".vcf.gz", ".control.vcf.gz")
    case_only_vcf_path = args.output.replace(".vcf.gz", ".case_only.vcf.gz")
    control_only_vcf_path = args.output.replace(".vcf.gz", ".control_only.vcf.gz")
    case_only_filtered_vcf_path = args.output.replace(".vcf.gz", ".case_only.filtered.vcf.gz")
    control_only_filtered_vcf_path = args.output.replace(".vcf.gz", ".control_only.filtered.vcf.gz")

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


if __name__ == "__main__":
    args = args_parser()
    main(args)