import copy
from collections import defaultdict
import os
import pandas as pd
import numpy as np
import sys
sys.path.append("/home/yz2296/bin")
import all_function_py3
# sys.path.append（./Run_NetFlow/all_function_py3.py)

# place this file outside of Run_NetFlow folder
dirname = os.path.dirname(os.path.abspath(__file__))
# S6_file = dirname + "/Run_NetFlow/pDenovoWest.xlsx"
# netflow_path = "/local/storage/yz2296/cancer_hotspot_new/website/website/pipeline/"
netflow_path = "/lustre2/home/js2788/Desktop/Lab/pipeline/"

# Function to read pDenovoWEST data into dataframe S6
# Input filename for pDenovoWEST data


# def read_S6_DeNovo(S6_file):

#     S6 = S6_file
#     S6 = S6[S6["pDenovoWEST"] != "."]
#     S6["pDenovoWEST"] = S6["pDenovoWEST"].astype(float)
#     S6["plog_DenovoWEST"] = S6["pDenovoWEST"].apply(lambda x: -np.log10(x))
#     S6.to_csv(dirname + "/Run_NetFlow/DeNovo_Input/pDNW_log.txt",
#               sep="\t", header=True, index=None)
#     return S6


# change from ensp to uniprots


def get_uniprot(ensp, REVEL_uniprot):
    UniProt = []
    for protein in ensp:
        uniprots = REVEL_uniprot[protein]
        if len(uniprots) == 0:
            UniProt.append(np.nan)
        else:
            uni = ",".join(sorted(uniprots))
            UniProt.append(uni)
    return UniProt

# Create dictionary for uniprot ID and corresponding residues


def res_dict(revel):
    residue_list = revel["Protein_position"].tolist()
    uniprot_list = revel["UniProt"].tolist()
    uniprot2pos = defaultdict(set)
    for i in range(revel.shape[0]):
        for prot in str(uniprot_list[i]).split(","):
            ###################################
            # uniprot2pos[prot] = uniprot2pos[prot].union({residue_list[i]})
            uniprot2pos[prot] = uniprot2pos[prot].union(
                set(all_function_py3.split_consecutive_pos(residue_list[i])))
            ###################################
    return uniprot2pos


# Function to determine whether the protein and its amino acid residues are on the interfaces
def whether_on_interface(row, uniprot2pos):

    if row["interface_residues1"] != "[]":
        prot1_ires = set([int(x)
                         for x in row["interface_residues1"][1:-1].split(",")])
    else:
        prot1_ires = set()
    if row["interface_residues2"] != "[]":
        prot2_ires = set([int(x)
                         for x in row["interface_residues2"][1:-1].split(",")])
    else:
        prot2_ires = set()

    if prot1_ires.intersection(uniprot2pos[row["uniprot1"]]) != set():
        return True
    elif prot2_ires.intersection(uniprot2pos[row["uniprot2"]]) != set():
        return True
    else:
        return False

# Function to produce interface table
# Include protein uniprot ID


def Interface_Table(Interface, uniprot2pos, output):
    Interface = Interface[Interface.apply(
        lambda x: whether_on_interface(x, uniprot2pos), axis=1)]
    Interface[["uniprot1", "uniprot2"]].rename(columns={"uniprot1": "P1", "uniprot2": "P2"}).to_csv(
        output, sep="\t", header=True, index=None)

    return Interface


# Function to create datatables for DeNovo based on pDenovoWEST p-values
# Input S6 = S6 data; p_value = p-value cutoff; score = "p" or "FDR"
# def p_filter_table(S11, p_value = None, score = "p", output_path = dirname):

# 	pscore = copy.deepcopy(S11)
# 	column = str(score)+ "_heat_score"
# 	pscore = S11[S11["FDR_TADA_NDD"] < p_value]

# 	output_path = output_path + "/Run_NetFlow/DeNovo_Input/DN_" + column + "_p_" + str(p_value) + ".txt"
# 	pscore[["gene_id",column]].rename(columns = {"gene_id":"Gene", column :"Score"}).to_csv(output_path, sep = "\t", header = True, index = None)
# 	return output_path

# def denovo_filter_table(S6, p_value, output_path = dirname):
# 	pscore = copy.deepcopy(S6)
# 	pscore["pDenovoWEST"] = pscore["pDenovoWEST"].astype(float)
# 	pscore = pscore[pscore["pDenovoWEST"] < p_value]

# 	output_path = output_path + "/Run_NetFlow/DeNovo_Input/DN_pDenovoWEST_" + str(p_value) + ".txt"
# 	pscore[["GeneID","plog_DenovoWEST"]].rename(columns = {"GeneID":"Gene", "plog_DenovoWEST" :"Score"}).to_csv(output_path, sep = "\t", header = True, index = None)
# 	return output_path

# Function to call NetFlow based on De Novo data and interface
# Input dataframe S6 which includes pDenovoWEST data and p values for ur, ar, and denovo
# Input dirname for path
def run_DN(dirname):
    # Parse DeNovo Variants
    revel = pd.read_csv(
        dirname + "/Run_NetFlow/new_final_combined_RunNF.txt", sep="\t", dtype=str)
    REVEL_uniprot = all_function_py3.gene_ID_conversion(
        revel["ENSP"], convert_from="ensembl.protein", convert_to="uniprot", species="human")
    revel["UniProt"] = get_uniprot(revel["ENSP"].tolist(), REVEL_uniprot)
    ###################################
    print("# of variants: " + str(revel.shape[0]))
    print("# of variants without UniProt annotation: " +
          str(revel["UniProt"].isna().sum()))
    # uniprot2pos = res_dict(revel)
    uniprot2pos = res_dict(revel.dropna(subset=["UniProt"]))
    ###################################

    # Import Pioneer data
    Interface = pd.read_csv(
        dirname + "/Run_NetFlow/PIONEER_human_veryhigh.txt", sep="\t")

    Interface_table = Interface_Table(
        Interface, uniprot2pos, dirname + "/Run_NetFlow/DeNovo_Input/new_Interface_DN.txt")
    return


run_DN(dirname)
