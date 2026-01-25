import os
import sys
import argparse
import logging
import pybedtools

import pandas as pd
sys.path.append(os.path.abspath("src/1_enrichment"))
from variants import VCFHandler, VCFBEDOverlapper
from utils import load_enhancer_dict, build_filters

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate TADA input files from variant data.")
    
    # global arguments
    parser.add_argument("--outdir", type=str, required=True, help="Directory for the output TADA input files.")
    parser.add_argument("--config", type=str, required=True, help="Path to the enhancer configuration file.")

    subparsers = parser.add_subparsers(dest="mode", required=True)

    # ----- build mode ------
    p_build = subparsers.add_parser("build",
        help="Build DN-only, INH-only, and UNTRANS-only intermediate tables from raw variant sources."
    )

    p_build.add_argument("--vcf", type=str, required=True, help="Path to the input variant file (VCF format).")
    p_build.add_argument("--type", type=str, choices=["dnv", "inh", "untrans"], required=True, help="Type of variants to process: 'dnv', 'inh', or 'untrans'.")
    p_build.add_argument("--tss-ref", type=str, required=True, help="Path to the gene annotation file.")
    p_build.add_argument("--enhancer-dir", type=str, default="data/enhancer/", help="Directory containing enhancer BED files.")
    p_build.add_argument("--enhancer-prefix", type=str, default="all_peaks_element_200bp_e", help="Prefix for enhancer BED files. " \
                                                "If provided, each enhancer file is named as <prefix>.bed")
    p_build.add_argument("--sex", type=str, choices=["male", "female", "any"], default="any", help="Sex of samples to consider.")
    p_build.add_argument("--gnomad-null-only", action="store_true", help="If set, require gnomADg_AF == '.'")
    p_build.add_argument("--gnomad-af-max", type=float, default=None, help="If set, require gnomADg_AF <= this value. Ignored if --gnomad-null-only is used.")
    p_build.add_argument("--cadd-min", type=int, default=None, help="Minimum CADD_PHRED for SNVs. If omitted, CADD filter is skipped.")
    p_build.add_argument("--phylop-min", type=int, default=None, help="Minimum phyloP100way for INDELs. If omitted, phyloP filter is skipped.")
    p_build.add_argument("--phastcons-min", type=float, default=None, help="Minimum phastCons100way for INDELs. If omitted, phastCons filter is skipped.")
    
    # ----- merge mode ------
    p_merge = subparsers.add_parser("merge",
        help="Merge De novo, Inherited, and Untransmitted TADA input files into a combined file."
    )

    p_merge.add_argument("--mut-rate-file", type=str, default=None, help="Path to the mutation rate file.")
    
    return parser.parse_args()

def enhancer_linked_genes(outdir: str, config: dict, tss_ref: str):
    """
    Link enhancers to genes based on proximity to TSS regions.
    """
    for name in config.keys():
        enhancer_bed_path = os.path.join(outdir, name, f"combined.merged.bed")

        enhancers_df = pd.read_csv(enhancer_bed_path, sep="\t", header=None
                        ).rename({0: "chrom", 1: "start", 2: "end"}, axis=1)
        enhancers_df["name"] = enhancers_df.apply(
            lambda row: f"{row['chrom']}:{row['start']}-{row['end']}", axis=1
        )
        enhancers_bt = pybedtools.BedTool.from_dataframe(enhancers_df)
        centers = enhancers_bt.each(
            lambda x: pybedtools.Interval(
                x.chrom,
                int((x.start + x.end) // 2),
                int((x.start + x.end) // 2) + 1,
                x.name
            )
        )
        slopped_bt = centers.slop(b=50000, genome='hg38'
                                  ).cut([0, 1, 2, 3])

        tss_df = pd.read_csv(tss_ref, sep="\t", header=None
                    )[[0,1,2,6,7]].rename(
                    columns={0: "chrom", 1: "start", 2: "end", 
                             6: "gene_name", 7: "gene_id"}
                    ).drop_duplicates().reset_index(drop=True)
        tss_bt = pybedtools.BedTool.from_dataframe(tss_df)
        
        linked_genes_df = slopped_bt.intersect(
                tss_bt, wa=True, wb=True
            ).cut([0,1,2,3,7,8]).to_dataframe(names=[
                "chrom", "start", "end", "enh_id",
                "gene_name", "gene_id"]).drop_duplicates().reset_index(drop=True)
        
        linked_genes_grouped = linked_genes_df.groupby(
                ["chrom", "start", "end", "enh_id"]
            ).agg({
                "gene_id": lambda x: ",".join(sorted(set(id.split(".")[0] for id in x))),
                "gene_name": lambda x: ",".join(sorted(set(x)))
            }).reset_index()

        renamed_enh_regions = os.path.join(outdir, name, "enh.merged.bed")
        enhancers_bt.saveas(renamed_enh_regions)
        # for variants mapping, keep original enhancer regions
        linked_genes_grouped["start"] = linked_genes_grouped["enh_id"].apply(lambda x: int(x.split(":")[1].split("-")[0]))
        linked_genes_grouped["end"] = linked_genes_grouped["enh_id"].apply(lambda x: int(x.split(":")[1].split("-")[1]))
        linked_genes_grouped.to_csv(enhancer_bed_path, sep="\t", index=False, header=True)

def generate_TADA_input(var_type, config, outdir):
    """
    Generate TADA input files from case and control variant files.

    Args:
        var_type (str): Type of variants to process ('dnv' or 'inh' or 'untrans').
        config (dict): Enhancer configuration dictionary.
        outdir (str): Output directory for TADA input files.
    """
    vart = "dn" if var_type == "dnv" else "in.t" if var_type == "inh" else "in.u"
    for name in config.keys():
        case_variants_path = os.path.join(outdir, name, "case.txt")
        ctrl_variants_path = os.path.join(outdir, name, "ctrl.txt")

        case_variants_df = pd.read_csv(case_variants_path, sep="\t", header=None
                            )[[0,1,2,3,4,5]].rename(
                            columns={0:"chrom", 1: "start", 2: "end", 
                                     3:"enh_id", 4: "gene_id", 5: "gene_name"}
                            ).groupby(["chrom", "start", "end", "enh_id"]).size().reset_index(name=f"{vart}.em")
        case_variants_df["chrom"] = case_variants_df["chrom"].apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else str(x))

        ctrl_variants_df = pd.read_csv(ctrl_variants_path, sep="\t", header=None
                            )[[0,1,2,3,4,5]].rename(
                            columns={0:"chrom", 1: "start", 2: "end", 
                                     3:"enh_id", 4: "gene_id", 5: "gene_name"}
                            ).groupby(["chrom", "start", "end", "enh_id"]).size().reset_index(name=f"{vart}.em.sib")
        ctrl_variants_df["chrom"] = ctrl_variants_df["chrom"].apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else str(x))

        merged_enh = pd.read_csv(os.path.join(outdir, name, "combined.merged.bed"), sep="\t", header=0)

        tada_input_df = merged_enh.merge(
            case_variants_df, on=["chrom", "start", "end", "enh_id"], how="left"
        ).merge(
            ctrl_variants_df, on=["chrom", "start", "end", "enh_id"], how="left"
        ).fillna(0)

        tada_input_df[[f"{vart}.em", f"{vart}.em.sib"]] = tada_input_df[[f"{vart}.em", f"{vart}.em.sib"]].astype(int)
        tada_input_df.to_csv(
            os.path.join(outdir, name, f"{var_type}.txt"), sep="\t", index=False, header=True
        )

def run_build(args):
    logging.info("Running in build mode.")
    outdir = f"{args.outdir}/{args.type}"

    vcf_overlapper = VCFBEDOverlapper(
        config=load_enhancer_dict(args.config),
        output_dir=outdir,
        enhancer_dir=args.enhancer_dir,
        enhancer_prefix=args.enhancer_prefix
    )

    logging.info("Generating merged enhancer BED files.")
    vcf_overlapper.build_merged_bed(names=list(vcf_overlapper.config.keys()))
    
    logging.info("Linking enhancers to genes - used 100kb TSS regions.")
    enhancer_linked_genes(outdir=outdir, config=vcf_overlapper.config, tss_ref=args.tss_ref)
    
    logging.info("Generating TADA input files in directory: %s", outdir)
    vcf_handler = VCFHandler()
    vcf_handler.set_vcf(args.vcf)
    case_vcf_path = os.path.join(outdir, "case.vcf.gz")
    ctrl_vcf_path = os.path.join(outdir, "ctrl.vcf.gz")
    base_inc_case, any_inc_case = build_filters(sex=args.sex, status="case", 
                                     gnomad_null_only=args.gnomad_null_only, 
                                     gnomad_af_max=args.gnomad_af_max,
                                     cadd_min=args.cadd_min, phylop_min=args.phylop_min,
                                     phastcons_min=args.phastcons_min)
    base_inc_ctrl, any_inc_ctrl = build_filters(sex=args.sex, status="control", 
                                     gnomad_null_only=args.gnomad_null_only, 
                                     gnomad_af_max=args.gnomad_af_max,
                                     cadd_min=args.cadd_min, phylop_min=args.phylop_min,
                                     phastcons_min=args.phastcons_min)
    logging.info("Filtering case variants - Base include: %s, " \
                 "Include any groups: %s", base_inc_case, any_inc_case)
    vcf_handler.pass_filters(output_path=case_vcf_path, vcf_path=args.vcf, include=base_inc_case,
        include_any=any_inc_case, threads=4, bgzip_index=True)
    logging.info("Filtering control variants - Base include: %s, " \
                 "Include any groups: %s", base_inc_ctrl, any_inc_ctrl)
    vcf_handler.pass_filters(output_path=ctrl_vcf_path, vcf_path=args.vcf, include=base_inc_ctrl,
        include_any=any_inc_ctrl, threads=4, bgzip_index=True)

    vcf_overlapper.map_variants_to_bed(
        vcf_path=case_vcf_path,
        names=list(vcf_overlapper.config.keys()),
        file="case.txt",
        write_overlap_vcf=True
    )

    vcf_overlapper.map_variants_to_bed(
        vcf_path=ctrl_vcf_path,
        names=list(vcf_overlapper.config.keys()),
        file="ctrl.txt",
        write_overlap_vcf=True
    )
    
    generate_TADA_input(
        var_type=args.type,
        config=vcf_overlapper.config,
        outdir=outdir
    )
    logging.info("Generated counts for %s.", args.type)

def run_merge(args):
    logging.info("Running in merge mode.")
    config = load_enhancer_dict(args.config)
    dnv_root = os.path.join(args.outdir, "dnv")
    inh_root = os.path.join(args.outdir, "inh")
    untrans_root = os.path.join(args.outdir, "untrans")
    merged_outdir = os.path.join(args.outdir, "input")
    os.makedirs(merged_outdir, exist_ok=True)

    if args.mut_rate_file:
        mut_rate_df = pd.read_csv(args.mut_rate_file, sep="\t", header=0)

    for name in config.keys():
        dnv_input_path = os.path.join(dnv_root, name, "dnv.txt")
        inh_input_path = os.path.join(inh_root, name, "inh.txt")
        unt_input_path = os.path.join(untrans_root, name, "untrans.txt")

        out_name = name.strip().replace(" ", "_")
        merged_input_path = os.path.join(merged_outdir, f"TADA_{out_name}.txt")

        dnv_df = pd.read_csv(dnv_input_path, sep="\t", header=0
                             ).drop(columns=["chrom", "start", "end"])
        inh_df = pd.read_csv(inh_input_path, sep="\t", header=0
                             ).drop(columns=["chrom", "start", "end", "in.t.em.sib"])
        unt_df = pd.read_csv(unt_input_path, sep="\t", header=0
                             ).drop(columns=["chrom", "start", "end", "in.u.em.sib"])

        keys = ["enh_id", "gene_id", "gene_name"]
        merged_df = dnv_df.merge(
                        inh_df, on=keys, how="outer"
                    ).merge(
                        unt_df, on=keys, how="outer"
                    ).fillna(0)
        
        merged_df[["chrom", "start", "end"]] = merged_df["enh_id"].str.split("[:-]", expand=True)
        merged_df["start"] = merged_df["start"].astype(int)
        merged_df["end"] = merged_df["end"].astype(int)

        if args.mut_rate_file:
            output_df = merged_df.merge(
                mut_rate_df, on=["chrom", "start", "end"], how="left"
            ).fillna(0).sort_values(by=["chrom", "start", "end"]).reset_index(drop=True)
            
            output_df = output_df[[
                "enh_id", "gene_id", "baselinerate_samocha14",
                "dn.em", "dn.em.sib",
                "in.t.em", "in.u.em",
            ]].rename(columns={"baselinerate_samocha14": "mut_rate"})
        else:
            output_df = merged_df.sort_values(by=["chrom", "start", "end"]).reset_index(drop=True)
            output_df = output_df[[
                "chrom", "start", "end",
                "enh_id", "gene_id",
                "dn.em", "dn.em.sib",
                "in.t.em", "in.u.em",
            ]]
        output_df.to_csv(merged_input_path, sep="\t", index=False, header=True)
        logging.info("Merged TADA input for enhancer set: %s", name)


def main(args):
    logging.basicConfig(level=logging.INFO)
    
    if args.mode == "build":
        run_build(args)
    elif args.mode == "merge":
        run_merge(args)
    else:
        logging.error("Invalid mode selected.")
        sys.exit(1)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)