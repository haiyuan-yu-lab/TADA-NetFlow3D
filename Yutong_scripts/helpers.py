""" 
Helper functions.
"""
import statistics
import logging
import os

import numpy as np
import pandas as pd

from collections import defaultdict

def consequence2varclass_coding(consequence):
    """
    consequence: "Consequence" column in VEP output 
    """
    if consequence in {"-", np.nan}:
        return np.nan
    else:
        conseq = set(consequence.replace("&", ",").split(","))
        if "frameshift_variant" in conseq:
            return "Frame_Shift_Del"
        elif "stop_gained" in conseq:
            return "Nonsense_Mutation"
        elif "stop_lost" in conseq:
            return "Nonstop_Mutation"
        elif "start_lost" in conseq:
            return "Translation_Start_Site"
        elif {"splice_acceptor_variant", "splice_donor_variant", "splice_donor_5th_base_variant"}.intersection(conseq) != set():
            return "Splice_Site"
        elif "missense_variant" in conseq:
            return "Missense_Mutation"
        elif "inframe_deletion" in conseq:
            return "In_Frame_Del"
        elif "inframe_insertion" in conseq:
            return "In_Frame_Ins"
        elif {"start_retained_variant", "synonymous_variant", "stop_retained_variant"}.intersection(conseq) != set():
            return "Silent"
        else:
            # logging.warning("Unknown consequence " + consequence + "!")
            return np.nan


def transform_variant(variant):
    """
    Function to transform variant to the desired format.
    """
    chrom, pos, ref, alt = variant.split(':')
    chrom = chrom.replace('chr', '')
    
    if len(alt) < len(ref):
        # Deletion case
        transformed = f"{chrom}_{str(int(pos)+1)}_{ref[1:]}/-"
    elif len(alt) > len(ref):
        # Insertion case
        transformed = f"{chrom}_{str(int(pos)+1)}_-/{alt[1:]}"
    else:
        # SNP or complex variant
        transformed = f"{chrom}_{pos}_{ref}/{alt}"
    
    return transformed


def Exclude_Sib(df, roles):
    """
    Function to exclude siblings with ASD risk
    """
    # get siblings who are too young for a reliable ASD diagnosis
    unknown_fam = set(roles[roles["family_type"] == "Unknown"]["family_sf_id"].tolist()) - {"SF0153305"} # SF0153305 is a special case, see slides on Mar 7, 2022
    pheno = pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/SSC+SPARK+KOR-STAR_parsed.txt", sep = "\t", dtype = str)
    young_sib = set(pheno[(pheno["Father"] != "0") & (pheno["Mother"] != "0") & (pheno["Pheno"] == "control") & pheno["Fam"].isin(unknown_fam)]["SampleID"])

    # get siblings with potential developmental concerns	  
    regis = pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/individuals_registration.csv")
    # regis = pd.read_csv("/fs/cbsuhy02/storage/yz2296/Burden_test_for_Autism/Phenotypes/SPARK_Collection_version5/individuals_registration.csv")
    concern_sib = set(regis[~regis["control_sibling_flag"].isna()]["subject_sp_id"].tolist())

    # exclude siblings with ASD risk
    return df[~df["SampleID"].isin(young_sib.union(concern_sib))]


def Exclude_Mul(df, roles):
    """
    Function to exclude Multiplex individuals
    """
    multiplex_fam = set(roles[roles["family_type"] == "Multiplex"]["family_sf_id"].tolist())
    return df[~df["FAM"].isin(multiplex_fam)]


# Function to parse REVEL scores after annotation
def parse_annotated_Rev(revel):
    s = []

    for x in revel.split(","):
        '''
        x = "ENST1;ENST2;ENST3:0.5"
        '''
        scores = x.split(":")
        if scores[-1] != "NaN":
            s.append(float(scores[-1]))
    if len(s) > 0:
        return max(s)

    else:
        return np.nan
    

def annotate_conservation_to_element(df, outdir, output = None, chr_col = "chr", start_col = "start", end_col = "end", start_col_index = 0, flag = "max"):
    """
    df: At least have 3 columns: chr, start, end. Specify the column names through **chr_col**, **start_col**, **end_col**.
    end_col: 0-indexed, not included by the element (hg38)
    start_col: 0-indexed or 1-indexed (hg38). If 1-indexed, sets **start_col_index** to 1. 
    flag: 1. max; 2. min; 3. median; (the way to aggregate the conservation scores of all positions covered by the same element)
    """
    def get_pos(row):
        for i in range(row[start_col], row[end_col] + 1):
            chr_list.append(row[chr_col])
            pos_list.append(i)
        return True
    def get_conserv(row):
        pos2phylop[(row["chr"], row["pos"])] = row["phyloP100way"]
        pos2phastcons[(row["chr"], row["pos"])] = row["phastCons100way"]
        return True
    def annot_conserv(row, score_dict):
        output = []
        for i in range(row[start_col], row[end_col] + 1):
            output.append(score_dict[(row[chr_col], i)])
        if flag == "max":
            return max(output)
        elif flag == "min":
            return min(output)
        else:
            return statistics.median(output)

    if flag not in {"max", "min", "median"}:
        logging.error("Wrong input flag!")
    else:
        # prepare element file for annotation
        df[start_col] = df[start_col].astype(int)
        df[end_col] = df[end_col].astype(int)
        if start_col_index == 0:
            df[start_col] = df[start_col] + 1

        # running conservation annotation pipeline
        idx = (df[start_col] > 0) & (df[end_col] > 0)
        df1 = df[idx]
        df2 = df[~idx]
        df2["phyloP100way"] = -np.inf
        df2["phastCons100way"] = -np.inf
        chr_list = []
        pos_list = []
        df_tmp = df1[[chr_col, start_col, end_col]].drop_duplicates()
        df_tmp["status"] = df_tmp.apply(get_pos, axis = 1)
        tmp_file = outdir+"/var_tmp"
        pd.DataFrame({"chr": chr_list, "pos": pos_list}).drop_duplicates().to_csv(tmp_file+".txt", sep = "\t", header = True, index = None)
        os.system("python3 annotate_conserv_score.py -v "+tmp_file+".txt -g hg38 -o "+tmp_file+"_output.txt")

        # annotate the input element file
        pos2phylop = defaultdict(lambda : -np.inf)
        pos2phastcons = defaultdict(lambda: -np.inf)
        df_output = pd.read_csv(tmp_file + "_output.txt", sep = "\t", dtype = {"pos": int})
        df_output["status"] = df_output.apply(get_conserv, axis = 1)
        df1["phyloP100way"] = df1.apply(lambda row: annot_conserv(row, pos2phylop), axis = 1)
        df1["phastCons100way"] = df1.apply(lambda row: annot_conserv(row, pos2phastcons), axis = 1)
        df = pd.concat([df1, df2])
        if start_col_index == 0:
            df[start_col] = df[start_col] - 1
        if output != None:
            df.to_csv(output, sep = "\t", header = True, index = None)

        # remove intermediate files
        os.remove(tmp_file + ".txt")
        os.remove(tmp_file + "_output.txt")
    return df