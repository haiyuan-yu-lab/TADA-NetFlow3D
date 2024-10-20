description = "Python script to preprocess three ASD datasets"

from io import StringIO
import pysam
import os
import sys
sys.path.append("/home/yz2296/bin")
import all_function_py3
import pandas as pd
import numpy as np
import importlib
import argparse
from os import path
# dirname = os.path.dirname(os.path.abspath(__file__))
dirname = "/fs/cbsuhy01/storage/NetFlow_ASD/Total_RunNetFlow"

if __name__ == "__main__":
    # user input
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-o','--output_path', type = str, default = "./Test_Output/", help = '[OPTIONAL] Path to the job folder. If not specified, the job folder will be stored in the default output folder')
    args = parser.parse_args()

    # Check input files
	# -----------------
    if args.output_path[-1] == "/":
        output_path = args.output_path
    else:
        output_path = args.output_path + "/"
    if path.exists(output_path) == False:
        os.mkdir(output_path)
    run_netflow = output_path + "Run_NetFlow/"
    if path.exists(run_netflow) == False:
        os.mkdir(run_netflow)


    roles = pd.read_excel(
        dirname + "/Run_NetFlow/roles_index.xlsx", index_col=None)

    # Function to parse MPC and CADD scores
    def filter_DMis(row):
        if row["CADD_PHRED"] != "-" and float(row["CADD_PHRED"]) >= 25:
            return True
        if row["MPC"] != "-" and float(row["MPC"]) >= 2:
            return True
        return False

    # Function to parse REVEL scores after annotation


    def parse_annotated_Rev(revel):
        s = []
        for x in revel.split(","):
            scores = x.split(":")
            if scores[-1] != "NaN":
                s.append(float(scores[-1]))
        if len(s) > 0:
            return max(s)
        else:
            return np.nan

    # Function to exclude Multiplex individuals
    def Exclude_Mul(df, roles, col):
        multiplex_asd = roles[roles["family_type"]
                            == "Multiplex"]["subject_sp_id"].tolist()
        return df[~df[col].isin(multiplex_asd)]

    # Function to parse revel scores
    def parse_revel(scores):
        output = []
        for x in scores.split(";"):
            if x != ".":
                output.append(float(x))
        if len(output) > 0:
            return max(output)
        else:
            return np.nan


    def vcf2txt(dirname, file_path, MPC=False):
        if MPC == True:
            cols = ["Consequence", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "HGVSc", "HGVSp", "Protein_position",
                    "Amino_acids", "Codons", "ENSP", "SIFT", "PolyPhen", "AF", "gnomADe_AF", "CADD_PHRED", "CADD_RAW", "MPC"]
        else:
            cols = ["Consequence", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "HGVSc", "HGVSp",
                    "Protein_position", "Amino_acids", "Codons", "ENSP", "SIFT", "PolyPhen", "AF", "gnomADe_AF", "CADD_PHRED", "CADD_RAW"]
        f = open(dirname + file_path + ".txt", "w")
        f.write("\t".join(["#Uploaded_variation"] + cols) + "\n")
        vcf = pysam.VariantFile(dirname + file_path + ".vcf")
        header = vcf.header.info["CSQ"].description.split("|")
        for variant in vcf:
            for annot in variant.info["CSQ"]:
                annot = annot.split("|")
                f.write("\t".join([variant.id] +
                        [annot[header.index(x)] for x in cols]) + "\n")
        f.close()


    def create_VEP(dirname):
        # Read in three supplementary tables
        DMis1 = pd.read_excel(
            dirname + "/Run_NetFlow/PMID35982159_supplementary_data1.xlsx", dtype=str)
        DMis1["DataSet"] = 1
        DMis2 = pd.read_excel(
            dirname + "/Run_NetFlow/PMID35982159_supplementary_data2.xlsx", dtype=str)
        DMis2["DataSet"] = 2
        DMis3 = pd.read_excel(
            dirname + "/Run_NetFlow/PMID35982159_supplementary_data3.xlsx", dtype=str)
        DMis3["DataSet"] = 3

        # Merge three tables and run vep
        DMis = pd.concat([DMis1, DMis2, DMis3])
        DMis["End_Position"] = DMis.apply(lambda x: int(
            x["Position"]) + len(x["Ref"]) - 1, axis=1)
        DMis["Mutation"] = DMis.apply(lambda x: str(
            x["Ref"]) + "/" + str(x["Alt"]), axis=1)
        DMis = DMis.sort_values(by=["Chrom", "End_Position"])
        VEP = DMis[["Chrom", "Position", "End_Position",
                    "Mutation"]].drop_duplicates()
        # Save dataframe as txt file
        VEP.to_csv(output_path + "/Run_NetFlow/Zhou_vep_input.txt",
                sep=" ", header=False, index=None)
        return DMis

    # Function to process hg19 data from Zhou


    def Parse_Zhou(dirname):
        # Read in merged raw supplementary table data
        DMis = create_VEP(dirname)
        DMis["#Uploaded_variation"] = DMis.apply(
            lambda x: x["Chrom"] + "_" + x["Position"] + "_" + x["Ref"] + "/" + x["Alt"], axis=1)
        DMis = DMis.drop(columns=["CADD13", "MPC"])

        # Read in vep results
        Zhou_file = "/Run_NetFlow/Zhou_vep_results"
        vcf2txt(dirname, Zhou_file, MPC=True)
        hg19 = pd.read_csv(
            dirname + "/Run_NetFlow/Zhou_vep_results.txt", sep="\t", dtype=str)
        hg19 = hg19[hg19["Protein_position"] != "-"]

        # Merge raw data and vep results
        merged = pd.merge(DMis, hg19, on="#Uploaded_variation", how="left")
        merged = merged[merged["Pheno"] == "Affected"]
        merged = merged.rename(
            columns={"Chrom": "chr", "Position": "pos", "Ref": "ref", "Alt": "alt"})

        # Filter supplementary table 3 using CADD and MPC scores
        merged_CADD = merged[merged["DataSet"] == 3]
        merged_CADD = merged_CADD[merged_CADD.apply(filter_DMis, axis=1)]

        # Save supplementary table 1 and 2 into separate file to annotate revel scores
    #     merged_REVEL = merged[merged["DataSet"] != 3].rename(columns = {"Chrom": "chr", "Position":"pos", "Ref":"ref", "Alt":"alt"})
        merged_REVEL = merged[merged["DataSet"] != 3]
        path = output_path + "/Run_NetFlow/hg19_preRevel.txt"
        merged_REVEL.to_csv(path, sep="\t", header=True, index=None)

        # Run parse_REVEL
        command = "python /fs/cbsuhyfs1/storage/resources/variant_effects/tools/Parse_REVEL.py -m " + path + \
            " -c hg19 -o " + dirname + \
            "/Run_NetFlow/parsed_REVEL_data/REVEL_annotated_file_Sup_hg19_REVEL.txt"
        os.system(command)

        # Read in revel-annotated supplementary table 1 and 2 and filter using revel scores
        merged_REVEL = pd.read_csv(
            dirname + "/Run_NetFlow/parsed_REVEL_data/REVEL_annotated_file_Sup_hg19_REVEL.txt", sep="\t", dtype=str)
        merged_REVEL = merged_REVEL.astype(str)
        merged_REVEL["parsed_revel"] = merged_REVEL["revel"].apply(
            parse_annotated_Rev)
        merged_REVEL = merged_REVEL[merged_REVEL["parsed_revel"].notna()]
        merged_REVEL["parsed_revel"] = merged_REVEL["parsed_revel"].astype(float)
        merged_REVEL = merged_REVEL[merged_REVEL["parsed_revel"] >= 0.5]

        # Merged all three parsed tables
        df_merge = pd.concat([merged_CADD, merged_REVEL])
        df_merge["GenomeBuild"] = "hg19"

        # Sort out the multiplex family
        df_merge = Exclude_Mul(df_merge, roles, "IID")
    #     df_merge = df_merge.rename(columns = {"Chrom": "chr", "Position":"pos", "Ref":"ref", "Alt":"alt"})
        df_merge.to_csv(output_path + "/Run_NetFlow/hg19_merged.txt",
                        sep="\t", header=True)
        return df_merge


    # Function to process hg38 data from Fu
    def Parse_Fu(dirname):
        # Read in raw data
        Sup20 = pd.read_excel(dirname + "/Run_NetFlow/PMID35982160_supplementary_table20.xlsx",
                            sheet_name="Supplementary Table 20", dtype=str)
        Sup20 = Sup20.dropna(subset=["Role"])

        Sup20[["Chrom", "Position", "Ref", "Alt"]] = pd.DataFrame(
            Sup20["Variant"].apply(lambda x: x.split(":")).tolist(), index=Sup20.index)
        Sup20["End_Position"] = Sup20.apply(lambda x: int(
            x["Position"]) + len(x["Ref"]) - 1, axis=1)
        Sup20["Mutation"] = Sup20.apply(lambda x: str(
            x["Ref"]) + "/" + str(x["Alt"]), axis=1)
        Sup20["#Uploaded_variation"] = Sup20.apply(
            lambda x: x["Chrom"] + "_" + x["Position"] + "_" + x["Ref"] + "/" + x["Alt"], axis=1)

        # Generate vep file to run vep
        VEP_20 = Sup20[["Chrom", "Position",
                        "End_Position", "Mutation"]].drop_duplicates()
        VEP_20 = VEP_20.sort_values(by=["Chrom", "End_Position"])
        VEP_20.to_csv(output_path + "/Run_NetFlow/Fu_vep_input.txt",
                    sep=" ", header=False, index=None)

        # Read in vep results and merge raw data with vep results
        Fu_file = "/Run_NetFlow/Fu_vep_results"
        vcf2txt(dirname, Fu_file)
        Sup20_results = pd.read_csv(
            dirname + "/Run_NetFlow/Fu_vep_results.txt", sep="\t", dtype=str)
        Sup20_merge = pd.merge(Sup20, Sup20_results,
                            on="#Uploaded_variation", how="left")
        Sup20_merge = Sup20_merge[Sup20_merge["Role"] == "Proband"]
        Sup20_merge = Sup20_merge[Sup20_merge["Protein_position"] != "-"]
        Sup20_merge["GenomeBuild"] = "grch38"

        # Exlude multiplex individuals
        Sup20_merge = Exclude_Mul(Sup20_merge, roles, "Sample")
        Sup20_merge = Sup20_merge[Sup20_merge.apply(filter_DMis, axis=1)]
        Sup20_merge = Sup20_merge.rename(
            columns={"Chrom": "chr", "Position": "pos", "Ref": "ref", "Alt": "alt"})

        Sup20_merge.to_csv(output_path + "/Run_NetFlow/Sup20_grch38_merged.txt",
                        sep="\t", header=True, index=None)
        return Sup20_merge


    def Parse_UCSF(dirname):
        # Read raw UCSF data
        UCSF_hg38 = pd.read_csv(
            dirname + "/Run_NetFlow/UCSF_20211001_dnms_WGS_hg38.tsv", sep="\t", dtype=str)
        UCSF_hg38["CHROM"] = UCSF_hg38["CHROM"].apply(lambda x: x[3:])
        UCSF_hg38 = UCSF_hg38.rename(
            columns={"CHROM": "chr", "POS": "pos", "REF": "ref", "ALT": "alt"})
        UCSF_hg38["End_Position"] = UCSF_hg38.apply(
            lambda x: int(x["pos"]) + len(x["ref"]) - 1, axis=1)
        UCSF_hg38["Mutation"] = UCSF_hg38.apply(
            lambda x: str(x["ref"]) + "/" + str(x["alt"]), axis=1)

        # Parse UCSF data for vep
        VEP_UCSF = UCSF_hg38[["chr", "pos",
                            "End_Position", "Mutation"]].drop_duplicates()
        VEP_UCSF = VEP_UCSF.sort_values(by=["chr", "End_Position"])
        VEP_UCSF.to_csv(output_path + "/Run_NetFlow/UCSF_vep_input.txt",
                        sep=" ", header=False, index=None)

        # Read UCSF VEP results and process into LOF and NT
        UCSF_file = "/Run_NetFlow/UCSF_vep_results"
        vcf2txt(dirname, UCSF_file)
        UCSF_vep = pd.read_csv(
            dirname + "/Run_NetFlow/UCSF_vep_results.txt", sep="\t", dtype=str)
        df_LOF = {"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation",
                "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site"}
        df_NT = {"Missense_Mutation", "In_Frame_Del", "In_Frame_Ins"}

        # Loop through Consequence and rename consequences
        UCSF_vep["Variant_Classification"] = UCSF_vep["Consequence"].apply(
            all_function_py3.consequence2varclass_coding)
        UCSF_hg38["#Uploaded_variation"] = UCSF_hg38.apply(
            lambda x: x["chr"] + "_" + x["pos"] + "_" + x["ref"] + "/" + x["alt"], axis=1)

        # Merge UCSF raw data with vep results
        UCSF_merged = pd.merge(UCSF_hg38, UCSF_vep,
                            on="#Uploaded_variation", how="left")
        UCSF_merged = UCSF_merged[UCSF_merged["STATUS"] == "case"]
        # UCSF_merged = UCSF_merged.rename(columns = {"Feature": "ENST"})

        UCSF_merged["GenomeBuild"] = "grch38"
        UCSF_merged = Exclude_Mul(UCSF_merged, roles, "SAMPLE")
        UCSF_merged = UCSF_merged[UCSF_merged["Protein_position"] != "-"]

        # Separate UCSF genes that are lof and nt
        UCSF_gene_lof = UCSF_merged[UCSF_merged["Variant_Classification"].isin(
            df_LOF)]
        UCSF_gene_nt = UCSF_merged[UCSF_merged["Variant_Classification"].isin(
            df_NT)]

        # Save data to txt files
        UCSF_merged.to_csv(output_path + "/Run_NetFlow/UCSF_merged.txt",
                        sep="\t", header=True, index=None)
        UCSF_gene_lof.to_csv(output_path + "/Run_NetFlow/UCSF_lof.txt",
                            sep="\t", header=True, index=None)
        path = output_path + "/Run_NetFlow/UCSF_nt.txt"
        UCSF_gene_nt.to_csv(path, sep="\t", header=True, index=None)

        # Run parse_REVEL
        command = "python /fs/cbsuhyfs1/storage/resources/variant_effects/tools/Parse_REVEL.py -m " + path + \
            " -c grch38 -o " + dirname + \
            "/Run_NetFlow/parsed_REVEL_data/REVEL_annotated_file_UCSF_hg38_REVEL.txt"
        os.system(command)

        # Read in REVEL annotated UCSF nt file
        final_UCSF = pd.read_csv(
            dirname + "/Run_NetFlow/parsed_REVEL_data/REVEL_annotated_file_UCSF_hg38_REVEL.txt", sep="\t", dtype=str)
        final_UCSF = final_UCSF.astype(str)
        final_UCSF["parsed_revel"] = final_UCSF["revel"].apply(parse_annotated_Rev)
        final_UCSF["parsed_revel"] = final_UCSF["parsed_revel"].astype(float)
        final_UCSF = final_UCSF[final_UCSF["parsed_revel"] >= 0.5]
        return final_UCSF


    final_hg19 = Parse_Zhou(dirname)
    final_Sup20 = Parse_Fu(dirname)
    final_UCSF = Parse_UCSF(dirname)
    final_df = pd.concat([final_hg19, final_Sup20, final_UCSF])
    final_df = final_df[final_df["Consequence"].apply(
        lambda x: x.find("missense_variant") != -1)]
    final_df.to_csv(output_path + "Run_NetFlow/new_final_combined_RunNF.txt",
                    sep="\t", header=True, index=None)
