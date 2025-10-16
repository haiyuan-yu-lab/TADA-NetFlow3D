import os
import argparse
import pandas as pd
from utils import dir_exists
from pints_pipeline import PintsPipeline

def args_parser():
    parser = argparse.ArgumentParser(description="Run PINTS peak calling and boundary extension.")
    parser.add_argument("--bw-pl", required=True, help="Path to the plus strand bigwig file.")
    parser.add_argument("--bw-mn", required=True, help="Path to the minus strand bigwig file.")
    parser.add_argument("--output-path", required=True, help="Directory to save output files.")
    parser.add_argument("--job-name", required=True, help="Prefix for output files.")
    parser.add_argument("--alpha", type=float, default=0.3, help="Stringency for clustering TSSs into a peak (default: 0.3).")
    parser.add_argument("--disable-sticks", action='store_true', help="Remove stick-like peaks (default: False).")
    parser.add_argument("--thread", type=int, default=16, help="Number of threads to use (default: 16).")
    parser.add_argument("--div-ext-left", type=int, default=60, help="Divergent extension left (default: 60).")
    parser.add_argument("--div-ext-right", type=int, default=60, help="Divergent extension right (default: 60).")
    parser.add_argument("--single-ext-left", type=int, default=60, help="Unidirectional extension left (default: 60).")
    parser.add_argument("--single-ext-right", type=int, default=60, help="Unidirectional extension right (default: 60).")
    parser.add_argument("--promoter-file", help="Path to the promoter BED file for boundary extension.")
    return parser.parse_args()

def main(args):
    dir_exists(args.output_path)

    # Initialize and run PINTS caller and boundary extender
    pints = PintsPipeline(bw_pl=args.bw_pl, bw_mn=args.bw_mn, peak_dir=args.output_path,
                          job_name=args.job_name, alpha=args.alpha, disable_sticks=args.disable_sticks,
                          thread=args.thread, promoter_file=args.promoter_file)
    pints.run_pints_caller()
    pints.run_pints_boundary_extender(div_ext_left=args.div_ext_left,
                                     div_ext_right=args.div_ext_right,
                                     single_ext_left=args.single_ext_left,
                                     single_ext_right=args.single_ext_right,
                                     promoter_file=args.promoter_file)

    # Merge elements
    pints.merge_and_save_elements(element_dir=args.output_path,
                                  div_ext_left=args.div_ext_left,
                                  div_ext_right=args.div_ext_right,
                                  single_ext_left=args.single_ext_left,
                                  single_ext_right=args.single_ext_right,
                                  enhancer=True, promoter=True)
    
    # Clear intermediate files
    pints.clear_intermediate_files()

if __name__ == "__main__":
    args = args_parser()
    main(args)