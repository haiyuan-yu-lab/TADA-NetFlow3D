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
    
    parser.add_argument("--out_path", 
                        help="Path to save output.",
                        type=str,
                        required=True,
                        )
    

def run_one_chr(df_p, df_var, output):
    df = pd.merge(df_var, df_p[["#Uploaded_variation", "Overlapped_elements_enhancer"]].drop_duplicates(), on = "#Uploaded_variation")
    df.to_csv(output, sep = "\t", header = None, index = None)
    return


def run_one_chr_multirun(args):
    return run_one_chr(*args)


def expand_column_to_rows(df_file):
    df = {"chr": [], "start" : [], "end" : [], "Overlapped_elements_enhancer": []}
    for key, df_key in df_file.groupby(["chr", "start", "end"]):
        element = []
        for ele in df_key["Overlapped_elements"]:
            for x in ele.split(","):
                element.append(x)
        for ele in set(element):
            df["chr"].append(key[0])
            df["start"].append(key[1])
            df["end"].append(key[2])
            df["Overlapped_elements_enhancer"].append(ele)

    df = pd.DataFrame(df)
    return df


def map_variants_to_elements(var_file, element_file, output_path):
    # Run get_overlapped_elements to map variants onto elements
    command = "python /fs/cbsuhy01/storage/NetFlow_ASD/Yutong_output/Enrichment/src/get_overlapped_elements.py -t 20 -o " + output_path + "/Parent_tmp_enhancer.txt -f " + var_file + " -r " + element_file
    os.system(command)
    var = pd.read_csv(output_path + "/Parent_tmp_enhancer.txt", sep = "\t", dtype = str)
    en = expand_column_to_rows(var).astype(str)
    var = pd.merge(var.drop(columns = ["Overlapped_elements"]), en, on = ["chr", "start", "end"])
    var.to_csv(output_path + "/UCSF_Parent_enhancer.txt", sep = "\t", header = True, index = None)
    os.remove(output_path + "/Parent_tmp_enhancer.txt")
    return


def main():
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
            "Microglia": {"Microglia"}
            "Neuron4R": {"Neuron4R"}, # needs update from Xiuqi file #TODO
            "Neuronal Progenitor Cells": {"LUHMES.Li"}, 
            "ReNcell CX Neural Progenitor cells": {"ReNcellCX"},

            "HEK293": {"HEK293"}, 
            "K562": {"K562"}
        
        2. Variants files in the format of chr, start, end, ref, alt (with header)
    - /fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/HEK293/promoter.v43.coding/UCSF_Parent_enhancer.txt
    - /fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/UCSF_Transmitted_ExcludeCoding_damaging.txt

    """
    # label2sample = {"Cerebral Cortex": {"CHTN22", "HN21"}, "Cerebellum": {"CHTN13"}, 
    # "Neuronal Progenitor Cells": {"LUHMES.Li"}, "Astrocytes": {"SC3"}, "Excitatory Cortical Neurons": {"SC2"}, 
    # "Combined": {"CHTN22", "HN21", "CHTN13", "LUHMES.Li", "SC3", "SC2", "Neuron4R", "Microglia"}, 
    # "K562": {"K562"}, "HEK293": {"HEK293"}, "Neuron4R": {"Neuron4R"}, "Microglia": {"Microglia"}}

    for label in ["HEK293"]:#, "ReNcellCX", ]:#["Combined", "K562"]:
        print("Processing " + label)
        y = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/" + label + "/promoter.v43.coding"
        map_variants_to_elements("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Variants/UCSF_Parent.txt.gz", y + "/all_peaks_element_200bp_e.bed", y)
        # print("Finished mapping variants to elements")

        df_parent = pd.read_csv(y + "/UCSF_Parent_enhancer.txt", sep = "\t", dtype = str)
        df_parent["#Uploaded_variation"] = df_parent["chr"].apply(lambda x: x[3:]) + "_" + df_parent["end"] + "_" + df_parent["ref"] + "/" + df_parent["alt"]
        df_parent_header = df_parent[0:1]
        df_parent = df_parent.groupby("chr")

        file = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/UCSF_Transmitted_ExcludeCoding_damaging.txt"
       
        df_var = pd.read_csv(file, sep = "\t", dtype = str)
        df_var_header = df_var[0:1]
        df_var = df_var.groupby("chr")
        header = pd.merge(df_var_header, df_parent_header[["#Uploaded_variation", "Overlapped_elements_enhancer"]].drop_duplicates(), on = "#Uploaded_variation").columns.tolist()
        with open(y + "/UCSF_Transmitted_ExcludeCoding_damaging_enhancer.txt", "w") as file:
            file.write("\t".join(header) + "\n")
        df_input = []
        for x in ["chr" + str(i) for i in range(1, 23)]:
            logging.info(x)
            df_p = df_parent.get_group(x)
            df_c = df_var.get_group(x)
            df_input.append([df_p, df_c, y + "/UCSF_Transmitted_ExcludeCoding_damaging_enhancer_" + x + ".txt"])
        p = Pool(11)
        p.map(run_one_chr_multirun, df_input)
        p.close()
        os.system("cat " + y + "/UCSF_Transmitted_ExcludeCoding_damaging_enhancer_" + "chr{1..22}.txt >> " + y + "/UCSF_Transmitted_ExcludeCoding_damaging_enhancer.txt")
        # os.system("rm " + y + "/" + os.path.basename(file).replace(".txt", "_" + z + "_chr*.txt"))

if __name__ == "__main__":
    main()