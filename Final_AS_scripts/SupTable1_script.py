description = "Python Script to generate Supplementary Table 1"

import pandas as pd
import numpy as np
import argparse
import sys
import importlib
import os
from os import path
from collections import defaultdict

if __name__ == "__main__":
	# user input
    parser = argparse.ArgumentParser(description = description)

    parser.add_argument('-p','--pheno_file', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8/SSC+SPARK+KOR-STAR_parsed.txt", help = 'Path to phenotype file.')
    parser.add_argument('-R','--raw_file_path', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/UCSF_raw/", help = 'Path to raw mutation files folder.')
    parser.add_argument('-o','--output_path', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/coding_analysis/test4/", help = '[OPTIONAL] Path to the job folder.')
    args = parser.parse_args()

    if os.path.isdir(args.output_path) == False:
        sys.exit("No such directory: " + args.output_path)
    if args.output_path[-1] == "/":
        dirname = args.output_path 
    else:
        dirname = args.output_path + "/" 

    if path.exists(args.pheno_file) == False:
        sys.exit("No such directory: " + args.pheno_file)
    if path.exists(args.raw_file_path) == False:
        sys.exit("No such directory: " + args.raw_file_path)



    def get_Table1(df, cohort="Total", dirname="/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/coding_analysis/test4"):
        if cohort == "Total":
            df_one = df
            
        elif cohort == "KOR":
            df_one = df[(df["SampleID"].str.startswith(str(cohort))) | (df["mother_id"].str.startswith("STAR"))]
            
        else:
            df_one = df[(df["SampleID"].str.startswith(str(cohort)))| (df["mother_id"].str.startswith(str(cohort)))]

        # Calculate unique counts
        proband = len(df_one[df_one["Pheno"] == "case"][["SampleID", "fam_id"]].drop_duplicates())
        sib = len(df_one[df_one["Pheno"] == "control"][["SampleID", "fam_id"]].drop_duplicates())
        mother = len(df_one["mother_id"].unique())
        father = len(df_one["father_id"].unique())

        # Create a dictionary for DataFrame
        df_total = {
            "Cohort" : [cohort],
            "Proband": [proband],
            "Sibling": [sib],
            "Mother": [mother],
            "Father": [father]
        }


        result_df = pd.DataFrame(df_total)
        return result_df

    Inherited = []

    for file in ["KOR-STAR.rare.inherited.tsv.gz", "SSC_WGS9209s.rare.inherited.tsv.gz", "SPARK_WGS_1.rare.inherited.tsv.gz", "SPARK_WGS_2.rare.inherited.tsv.gz", "SPARK_WGS_3.rare.inherited.tsv.gz"]:
        Inherited.append(pd.read_csv(args.raw_file_path + file, sep = "\t",  dtype = str))
    Inherited = pd.concat(Inherited).rename(columns = {"sample_id": "SampleID"}).drop(columns = ["gnomad_AF"])
    pheno = pd.read_csv(args.pheno_file, sep = "\t", dtype = str)
    Inherited = pd.merge(Inherited, pheno, on = ["SampleID"], how = "left")
    table_total = pd.DataFrame()
    for cohort in ["KOR", "SSC", "SP", "Total"]:
        table_total = pd.concat([table_total, get_Table1(Inherited, cohort = cohort)])

    table_total["Total"] = table_total[["Proband", "Sibling", "Mother", "Father"]].sum(axis=1)
    table_total.to_csv(dirname + "SupplmentaryTable1.txt", sep = "\t", header = True, index = None)