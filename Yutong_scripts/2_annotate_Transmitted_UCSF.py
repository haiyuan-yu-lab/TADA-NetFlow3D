import os
import logging
import argparse

import pandas as pd
from multiprocessing import Pool


def set_parser(parser):
    parser.add_argument("--enhancer_dir",
                        help="Path to the folder containing enhancer subfolders.",
                        type=str,
                        required=True,
                        )
    
    parser.add_argument("--preprocess_dir",
                        help="Path to preprocessed UCSF data.",
                        type=str,
                        required=True,
                        )
    
    parser.add_argument("--cell_line",
                        help="Cell line to be processed.",
                        type=str,
                        required=True,
                        nargs="+",
                        )
    
    parser.add_argument("--out_path", 
                        help="Path to save output.",
                        type=str,
                        required=True,
                        )


def get_all_enhancer_elements(args):
    """
    Get all enhancer elements from different cell lines and save them into one file.
    """
    all_enhancer_elements = []
    for label in args.cell_line:
        print("Processing " + label)
        enhancer_dir = args.enhancer_dir + label + "/promoter.v43.coding"
        enhancer = pd.read_csv(f"{enhancer_dir}/all_peaks_element_200bp_e.bed", sep="\t", header=0)
        all_enhancer_elements.append(enhancer)

    enhancer_df = pd.concat(all_enhancer_elements)

    command = "mkdir -p " + args.enhancer_dir + "/ALL_ENHANCER"
    os.system(command)
    enhancer_df.to_csv(f"{args.enhancer_dir}/ALL_ENHANCER/all_peaks_element_200bp_e.bed", sep="\t", header=True, index=False)

    return enhancer_df


def main(args):
    """
    Save the transmitted variants that overlapped with enhancers in different cell lines.
    Extra step compared to DNV variants: VEP not saved in total file??

    Required files:
        1. PROcap elements in the format of chr, start, end (with header)
        under the folder /fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/ 
        with corresponding cell line names:
            "Astrocytes": {"SC3"}, 
            "Cerebral Cortex": {"CHTN22", "HN21"},
            "Cerebellum": {"CHTN13"},
            "Excitatory Cortical Neurons": {"SC2"}, 
            "H9-derived neural crest cells": {"SC8"},
            "Microglia": {"iMG_WT"},
            "Neuron4R": {"Neuron4R"}, 
            "Neuronal Progenitor Cells": {"LUHMES.Li"}, 
            "ReNcell CX Neural Progenitor cells": {"ReNcellCX"},

            "HEK293": {"HEK293"}, 
            "K562": {"K562"}
        
        2. Variants files in the format of chr, start, end, ref, alt (with header)

    """
    enhancer_df = get_all_enhancer_elements(args)
    enhancer_dir = args.enhancer_dir + "/ALL_ENHANCER"

    # map_variants_to_elements("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Variants/UCSF_Parent.txt.gz", \
    #                             f"{enhancer_dir}/all_peaks_element_200bp_e.bed", enhancer_dir)
    # print("Finished mapping variants to elements")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    set_parser(parser)
    args = parser.parse_args()
    main(args)