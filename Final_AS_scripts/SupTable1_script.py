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



    # Function to generate a summary table (Table 1) based on family information from a given cohort
    def get_Table1(df, cohort="Total", dirname="/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/coding_analysis/test4"):
        # Filter the dataframe based on the cohort type
        if cohort == "Total":
            df_one = df  # If cohort is "Total", use the entire dataframe
        elif cohort == "KOR":
            # If cohort is "KOR", select samples that either start with "KOR" or have mothers starting with "STAR"
            df_one = df[(df["SampleID"].str.startswith(str(cohort))) | (df["mother_id"].str.startswith("STAR"))]
        else:
            # For other cohorts, select samples that either start with the cohort name or have mothers starting with the cohort name
            df_one = df[(df["SampleID"].str.startswith(str(cohort))) | (df["mother_id"].str.startswith(str(cohort)))]

        # Calculate the unique counts for different family members in the cohort
        proband = len(df_one[df_one["Pheno"] == "case"][["SampleID", "fam_id"]].drop_duplicates())  # Count unique probands (cases)
        sib = len(df_one[df_one["Pheno"] == "control"][["SampleID", "fam_id"]].drop_duplicates())  # Count unique siblings (controls)
        mother = len(df_one["mother_id"].unique())  # Count unique mothers
        father = len(df_one["father_id"].unique())  # Count unique fathers

        # Create a dictionary to store the cohort summary information
        df_total = {
            "Cohort": [cohort],
            "Proband": [proband],
            "Sibling": [sib],
            "Mother": [mother],
            "Father": [father]
        }

        # Convert the dictionary into a DataFrame
        result_df = pd.DataFrame(df_total)
        return result_df  # Return the summary DataFrame

    # Initialize an empty list to store the inherited variant data
    Inherited = []

    # Loop through each inherited variant file and read it into a dataframe
    for file in ["KOR-STAR.rare.inherited.tsv.gz", "SSC_WGS9209s.rare.inherited.tsv.gz", "SPARK_WGS_1.rare.inherited.tsv.gz", "SPARK_WGS_2.rare.inherited.tsv.gz", "SPARK_WGS_3.rare.inherited.tsv.gz"]:
        # Read each file into a dataframe and append to the Inherited list
        Inherited.append(pd.read_csv(args.raw_file_path + file, sep="\t", dtype=str))

    # Concatenate all dataframes in the Inherited list into a single dataframe and rename the column "sample_id" to "SampleID"
    Inherited = pd.concat(Inherited).rename(columns={"sample_id": "SampleID"}).drop(columns=["gnomad_AF"])

    # Read the phenotype data from a file
    pheno = pd.read_csv(args.pheno_file, sep="\t", dtype=str)

    # Merge the inherited variant data with the phenotype data based on "SampleID"
    Inherited = pd.merge(Inherited, pheno, on=["SampleID"], how="left")

    # Initialize an empty dataframe to store the summary table for all cohorts
    table_total = pd.DataFrame()

    # Generate a summary table for each cohort and concatenate to table_total
    for cohort in ["KOR", "SSC", "SP", "Total"]:
        table_total = pd.concat([table_total, get_Table1(Inherited, cohort=cohort)])

    # Calculate the total counts for Proband, Sibling, Mother, and Father across all cohorts and add it as a new column "Total"
    table_total["Total"] = table_total[["Proband", "Sibling", "Mother", "Father"]].sum(axis=1)

    # Save the summary table to a file as "SupplmentaryTable1.txt"
    table_total.to_csv(dirname + "SupplmentaryTable1.txt", sep="\t", header=True, index=None)
