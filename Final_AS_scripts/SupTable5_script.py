description = "Python Script to generate coding count analysis table"

import argparse
import sys
sys.path.append("/home/yz2296/bin")
import all_function_py3
import os
from os import path
import logging
from multiprocessing import Pool
import time
from joblib import Parallel, delayed
import gc
import pandas as pd
import numpy as np
directory = os.path.dirname(os.path.abspath(__file__))

# Parallelize get_untransmitted_var
def get_untransmitted_var_multirun(kwargs):
    return get_untransmitted_var(**kwargs)


if __name__ == "__main__":
    # user input
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-s','--sib_file', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/individuals_registration.csv", help = 'Path to sibling information file.')
    parser.add_argument('-R','--raw_file_path', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/UCSF_raw/", help = 'Path to raw mutation files folder.')
    parser.add_argument('-r','--roles_file', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8/roles_index.csv", help = 'Path to roles file.')
    parser.add_argument('-g','--sample_information', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8/SSC+SPARK+KOR-STAR_parsed.txt", help = 'Path to sex file.')

    parser.add_argument('-m','--Mutability_file', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/coding_analysis/41588_2022_1104_MOESM3_ESM.xlsx", help = 'Path to Mutability file.')
    parser.add_argument('-o','--output_path', type = str, default = "./Test_Output/", help = '[OPTIONAL] Path to the job folder. If not specified, the job folder will be stored in the default output folder')

    parser.add_argument('-I','--job_name', required = True, type = str, help = 'Please specify a name for your job. The output files will be stored with this name')
    args = parser.parse_args()
############################################################################################################################################
	# Check input files
	# -----------------
    if os.path.isdir(args.output_path) == False:
        sys.exit("No such directory: " + args.output_path)
    if args.output_path[-1] == "/":
        dirname = args.output_path + args.job_name 
    else:
        dirname = args.output_path + "/" + args.job_name
    if path.exists(dirname) == False:
        os.mkdir(dirname)

    output_path = dirname + "/Summary/"
    if path.exists(output_path) == False:
        os.mkdir(output_path)

    if path.exists(args.sib_file) == False:
        sys.exit("No such directory: " + args.sib_file)
    if path.exists(args.roles_file) == False:
        sys.exit("No such directory: " + args.roles_file)
    if path.exists(args.sample_information) == False:
        sys.exit("No such directory: " + args.sample_information)
    if path.exists(args.raw_file_path) == False:
        sys.exit("No such directory: " + args.raw_file_path)
    if path.exists(args.Mutability_file) == False:
        sys.exit("No such directory: " + args.Mutability_file)
    

    lof = {"Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Nonstop_Mutation","Translation_Start_Site","Splice_Site"}


    # Function to exclude Multiplex individuals
    def Exclude_Mul(df, roles):
        multiplex_fam = set(roles[roles["family_type"] == "Multiplex"]["family_sf_id"].tolist())
        return df[~df["FAM"].isin(multiplex_fam)]

    def get_var_type(row):
        # Check if either 'ref' or 'alt' is NaN and return np.nan if true
        if pd.isna(row['ref']) or pd.isna(row['alt']):
            return np.nan
        elif len(row["ref"]) == len(row["alt"]) and set(row["ref"] + row["alt"]) - {"A", "C", "G", "T"} == set():
            return "SNV"
        else:
            return "Indel"


    # Function to exclude siblings with ASD risk
    def Exclude_Sib(df, roles):
        # get siblings who are too young for a reliable ASD diagnosis
        unknown_fam = set(roles[roles["family_type"] == "Unknown"]["family_sf_id"].tolist()) - {"SF0153305"} # SF0153305 is a special case, see slides on Mar 7, 2022
        pheno = pd.read_csv(args.sample_information, sep = "\t", dtype = str)
        young_sib = set(pheno[(pheno["Father"] != "0") & (pheno["Mother"] != "0") & (pheno["Pheno"] == "control") & pheno["Fam"].isin(unknown_fam)]["SampleID"])

        # get siblings with potential developmental concerns
        regis = pd.read_csv(args.sib_file)
        # regis = pd.read_csv("/fs/cbsuhy02/storage/yz2296/Burden_test_for_Autism/Phenotypes/SPARK_Collection_version5/individuals_registration.csv")
        concern_sib = set(regis[~regis["control_sibling_flag"].isna()]["subject_sp_id"].tolist())

        # exclude siblings with ASD risk
        return df[~df["SampleID"].isin(young_sib.union(concern_sib))]

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

    #Function to filter gnomad cutoff
    def gnomad_filter(df, gnomad_col = "gnomADe_AF", gnomAD = 0):
        ############################################################################################################################################
        # df_NC = df_NC[df_NC["gnomADg_AF"] != "-"]
        ############################################################################################################################################
        df_in_gnomad = df[df[gnomad_col] != "-"]
        df_notin_gnomad = df[df[gnomad_col] == "-"]
        
        df_in_gnomad[gnomad_col] = df_in_gnomad[gnomad_col].astype(float)
        df_in_gnomad = df_in_gnomad[df_in_gnomad[gnomad_col] <= float(gnomAD)]
        df_filtered = pd.concat([df_in_gnomad, df_notin_gnomad])
        return df_filtered


    # Function to preprocess de novo variants
    def DNV_preprocess(dirname, coding_flag = "damaging"):
        '''
        coding_flag: 1. "damaging"; 2. "strict"; 3. "loose"; 4. None.
            1. "damaging": exclude LOF & damaging missense variants
            2. "strict": exclude variants with ANY coding consequences
            3. "loose": exclude variants with ONLY coding consequences
            4. None: do not exclude any coding variants
        '''
        # Generate VEP input
        UCSF_hg38 = pd.read_csv(args.raw_file_path + "UCSF_20211001_dnms_WGS_hg38.tsv", sep = "\t").sort_values(by = ["CHROM", "POS"])
        UCSF_hg38["vep"] = UCSF_hg38.apply(lambda x: " ".join([x["CHROM"][3:], str(x["POS"]), str(x["POS"] + len(x["REF"]) - 1), x["REF"] + "/" + x["ALT"]]), axis = 1)
        UCSF_hg38[["vep"]].drop_duplicates().to_csv(dirname + "/UCSF_DNV_vep_input.txt", sep = "\t", header = None, index = None)

        # # Load and process VEP output
        all_function_py3.parse_VEP_output("/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/coding_analysis/UCSF_DNV_vep_output.vcf", dirname + "/UCSF_DNV_vep_output.txt", cols = ["Consequence", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "HGVSc", "HGVSp", "Protein_position", "Amino_acids", "Codons", "ENSP", "SIFT", "PolyPhen", "CADD_PHRED", "CADD_RAW", "gnomADg_AF", "gnomADe_AF"])
        UCSF_vep = pd.read_csv(dirname + "/UCSF_DNV_vep_output.txt", sep = "\t", dtype = str)
        UCSF_vep["Variant_Classification"] = UCSF_vep["Consequence"].apply(all_function_py3.consequence2varclass_coding)

        # load and process raw data
        UCSF_hg38 = pd.read_csv(args.raw_file_path + "UCSF_20211001_dnms_WGS_hg38.tsv", sep = "\t", dtype = str)
        UCSF_hg38["CHROM"] = UCSF_hg38["CHROM"].apply(lambda x: x[3:])
        UCSF_hg38 = UCSF_hg38.rename(columns = {"CHROM": "chr", "POS":"pos", "REF":"ref", "ALT":"alt", "STATUS": "Pheno", "SAMPLE": "SampleID"})
        UCSF_hg38["#Uploaded_variation"] = UCSF_hg38.apply(lambda x: x["chr"] + "_" + x["pos"] + "_" + x["ref"] + "/" + x["alt"], axis = 1)

        # Annotate raw data
        UCSF_merged = pd.merge(UCSF_hg38, UCSF_vep, on = "#Uploaded_variation")

        # Exclude siblings with ASD risk
        roles = pd.read_csv(args.roles_file, index_col = None)
        UCSF_merged = Exclude_Sib(UCSF_merged, roles)
        # Exclude multiplex families
        UCSF_merged = Exclude_Mul(UCSF_merged, roles)
        UCSF_merged["chr"] = UCSF_merged["chr"].apply(lambda x: "chr" + str(x))    
        UCSF_merged["start"] = UCSF_merged.apply(lambda x: int(x["pos"]) -1 , axis=1)
        UCSF_merged["end"] = UCSF_merged.apply(lambda row: int(row["pos"]) + len(row["ref"]) - 1, axis = 1)
    
        if coding_flag == "damaging":
            UCSF_lof = UCSF_merged[UCSF_merged["Variant_Classification"].isin(lof)]
            UCSF_lof.to_csv(dirname + "/UCSF_DNV_preprocessed_LOF.txt", sep = "\t", header = True, index = None)
            UCSF_merged[UCSF_merged["Variant_Classification"] == "Missense_Mutation"].to_csv(dirname + "/UCSF_DNV_preprocessed_missense.txt", sep = "\t", header = True, index = None)
            # command = "python /fs/cbsuhyfs01/storage/resources/variant_effects/tools/Parse_REVEL.py -m mis_tmp.txt -c grch38 -o " + dirname + "/mis_tmp_annotated.txt"
            # os.system(command)
            output = dirname + "/UCSF_DNV_preprocessed_missense_hg19.txt"
            input = dirname + "/UCSF_DNV_preprocessed_missense.txt"
            command = "python /fs/cbsuhy01/storage/resources/liftOver/coord_conversion_between_genome_assemblies.py -o " + output + " -f hg38 -t hg19 -e " + input + " -c chr -s start -E end"
            os.system(command)
            
    #         all_function_py3.parse_VEP_output("/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/coding_analysis/UCSF_missense_hg19_vep_output.vcf", dirname + "/UCSF_missense_hg19_vep_output.txt", cols = ["Consequence", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "HGVSc", "HGVSp", "Protein_position", "Amino_acids", "Codons", "ENSP", "SIFT", "PolyPhen", "CADD_PHRED", "CADD_RAW", "gnomADg_AF", "MPC"])
    #         UCSF_mis = pd.read_csv(dirname + "/UCSF_missense_hg19_vep_output.txt", sep = "\t")
    #         # UCSF_mis["parsed_revel"] = UCSF_mis["revel"].apply(parse_annotated_Rev)
    #         # UCSF_Dmis = UCSF_mis[UCSF_mis["parsed_revel"] >= 0.5]

        return output, str(dirname + "/UCSF_DNV_preprocessed_LOF.txt")

    # Function to load and process file in parallel for transmitted
    def load_and_process_file(file):
        df = pd.read_csv(file, sep="\t", dtype=str)
        df = df.rename(columns={"sample_id": "SampleID"}).drop(columns=["gnomad_AF"])
        df[["chr", "pos"]] = df["locus"].str.split(":", expand=True)
        df["#Uploaded_variation"] = df["chr"].str[3:] + "_" + df["pos"] + "_" + df["ref"] + "/" + df["alt"]
        df["start"] = df["pos"].astype(int) - 1
        df["end"] = df.apply(lambda x: int(x["pos"]) + len(x["ref"]) - 1, axis=1)
        return df

    # Function to preprocess transmitted data
    def Transmitted_preprocess(dirname, coding_flag="damaging"):
        '''
        coding_flag: 1. "damaging"; 2. "strict"; 3. "loose"; 4. None.
            1. "damaging": exclude LOF variants
            2. "strict": exclude variants with ANY coding consequences
            3. "loose": exclude variants with ONLY coding consequences
            3. None: do not exclude any coding variants
        '''
        files = [
            "KOR-STAR.rare.inherited.tsv.gz",
            "SSC_WGS9209s.rare.inherited.tsv.gz",
            "SPARK_WGS_1.rare.inherited.tsv.gz",
            "SPARK_WGS_2.rare.inherited.tsv.gz",
            "SPARK_WGS_3.rare.inherited.tsv.gz"
        ]
        # Load and process raw data in parallel
        Inherited = Parallel(n_jobs=5)(delayed(load_and_process_file)(args.raw_file_path + file) for file in files)
        Inherited = pd.concat(Inherited)
        logging.warning("transmitted line original: " + str(Inherited[["#Uploaded_variation", "SampleID"]].drop_duplicates().shape[0]))

        # Load and process variant annotation
        var_annot = pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Variants/UCSF_parent_var_annotation.txt", sep="\t", dtype=str)
        Inherited = Inherited.merge(var_annot, on="#Uploaded_variation")
        logging.warning("transmitted line after mapping: " + str(Inherited[["#Uploaded_variation", "SampleID"]].drop_duplicates().shape[0]))

        # Exclude siblings with ASD risk
        roles = pd.read_csv(args.roles_file, index_col=None)
        Inherited = Exclude_Sib(Inherited, roles)

        # Merge with phenotype data
        pheno = pd.read_csv(args.sample_information, sep="\t", dtype=str)
        Inherited = Inherited.merge(pheno, on=["SampleID"], how="left")

        Inherited["Variant_Classification"] = Inherited["Consequence"].apply(all_function_py3.consequence2varclass_coding)

        if coding_flag == "damaging":
            Inherited_lof = Inherited[Inherited["Variant_Classification"].isin(lof)]
            Inherited_lof.to_csv(dirname + "/UCSF_Transmitted_preprocessed_LOF.txt", sep="\t", header=True, index=None)
            
            # Prepare input for VEP
            lof_vep_input_path = dirname + "/UCSF_Transmitted_lof_vep_input.txt"
            Inherited_lof["vep_input"] = Inherited_lof.apply(lambda x: f"{x['chr'][3:]} {x['pos']} {x['end']} {x['ref']}/{x['alt']}", axis=1)
            Inherited_lof["vep_input"].drop_duplicates().to_csv(lof_vep_input_path, sep="\t", header=None, index=None)

            # Process Missense Mutations
            missense_path = dirname + "/UCSF_Transmitted_preprocessed_missense.txt"
            Inherited[Inherited["Variant_Classification"] == "Missense_Mutation"].to_csv(missense_path, sep="\t", header=True, index=None)
            
            # External command for liftover
            output = dirname + "/UCSF_Transmitted_preprocessed_missense_hg19.txt"
            command = f"python /fs/cbsuhy01/storage/resources/liftOver/coord_conversion_between_genome_assemblies.py -o {output} -f hg38 -t hg19 -e {missense_path} -c chr -s start -E end"
            os.system(command)

        lof_vep_output = dirname + "/UCSF_Transmitted_lof_vep_output.txt"
        all_function_py3.parse_VEP_output("/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/coding_analysis/UCSF_Transmitted_lof_vep_output.vcf", lof_vep_output, cols = ["Consequence", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "HGVSc", "HGVSp", "Protein_position", "Amino_acids", "Codons", "ENSP", "SIFT", "PolyPhen", "CADD_PHRED", "CADD_RAW",  "gnomADe_AF"])
        
        UCSF_Transmitted_lof_preprocessed = pd.read_csv(dirname + "/UCSF_Transmitted_preprocessed_LOF.txt", sep="\t")
        df_lof_vep_output = pd.read_csv(lof_vep_output, sep="\t")
        Transmitted_merged_lof = UCSF_Transmitted_lof_preprocessed.merge(df_lof_vep_output, on=["#Uploaded_variation"], how="left")
        UCSF_Transmitted_lof_path = dirname + "/UCSF_Transmitted_lof_vep_merged.txt"
        Transmitted_merged_lof.to_csv(UCSF_Transmitted_lof_path, sep="\t", header=True, index=None)
        
        return output, UCSF_Transmitted_lof_path

    def get_untransmitted_var(row, df_parent, inherited_var, output):
        father_var = df_parent[df_parent["SampleID"] == row["Father"]]
        father_var["sex"] = "male"
        mother_var = df_parent[df_parent["SampleID"] == row["Mother"]]
        mother_var["sex"] = "female"
        parent_var = pd.concat([father_var, mother_var])
        # parent_var["TYPE"] = parent_var.apply(get_var_type, axis = 1)
        if len(parent_var["SampleID"].unique()) == 2 and inherited_var.shape[0] > 0:

            # get untransmitted variants
            untransmitted_df = parent_var[~parent_var["#Uploaded_variation"].isin(set(inherited_var["#Uploaded_variation"]))]
            untransmitted_df["Pheno"] = row["Pheno"]
            untransmitted_df["Offspring_SampleID"] = row["SampleID"]

            untransmitted_df.to_csv(output, mode = "a", sep = "\t", header = False, index = False)

            base, extension = os.path.splitext(output)
            vep_input = base + "_vep_input" + extension
            untransmitted_df.apply(lambda x: " ".join([x["chr"][3:], str(x["pos"]), str(int(x["pos"]) + len(x["ref"]) - 1), x["ref"] + "/" + x["alt"]]), axis = 1).drop_duplicates().to_csv(vep_input, mode = "a", sep = "\t", header = None, index = None)
        else:
            logging.warning("The couple " + row["Father"] + " & " + row["Mother"] + " have incomplete variant data!")
        return True

    # Function to preprocess untransmitted data
    def Untransmitted_preprocess(chrom, df_parent, df_annot, dirname, coding_flag = "damaging"):

        def get_function_input(x):
            df_parent_tmp = df_parent[df_parent["SampleID"].isin({x["Father"], x["Mother"]})]
            df_inherited_tmp = df_inherited[df_inherited["SampleID"] == x["SampleID"]]
            function_input.append({"row": dict(x), "df_parent": df_parent_tmp, "inherited_var": df_inherited_tmp, "output": preprocessed_output})
            return True

        # get variants
        df_parent = pd.merge(df_parent, df_annot, on = "#Uploaded_variation")
        df_inherited_missense = pd.read_csv(dirname + "/Transmitted/UCSF_Transmitted_preprocessed_missense.txt", sep = "\t", dtype = str)
        df_inherited_lof = pd.read_csv(dirname + "/Transmitted/UCSF_Transmitted_preprocessed_LOF.txt", sep = "\t", dtype = str)
        df_inherited = pd.concat([df_inherited_missense, df_inherited_lof])

        if coding_flag == "damaging":
            parent_lof = df_parent[df_parent["Variant_Classification"].isin(lof)]
            parent_missense = df_parent[df_parent["Variant_Classification"] == "Missense_Mutation"]
            parent_total = pd.concat([parent_lof, parent_missense])
            parent_lof.to_csv(dirname + "/Untransmitted/UCSF_Parent_LOF_" + chrom + ".txt", sep = "\t", header = True, index = None)
            parent_missense.to_csv(dirname + "/Untransmitted/UCSF_Parent_missense_" + chrom + ".txt", sep = "\t", header = True, index = None)
            # df_parent = df_parent[df_parent["#Uploaded_variation"].isin(set(parent_total["#Uploaded_variation"]))]
        
        # get untransmitted variants
        preprocessed_output = dirname + "/Untransmitted/UCSF_Untransmitted_preprocessed_LOF_missense_" + chrom + ".txt"
        pd.DataFrame({x: [] for x in parent_total.columns.tolist() + ["sex", "TYPE", "Pheno", "Offspring_SampleID"]}).to_csv(preprocessed_output, sep = "\t", header = True, index = None)
        df = pd.read_csv(args.sample_information, sep = "\t", dtype = str)
        df = df[df["Father"] != "0"]
        function_input = []
        df["status"] = df.apply(get_function_input, axis = 1)

        # Parallel processing with controlled resource usage
        with Pool(5) as p:
            results = p.map(get_untransmitted_var_multirun, function_input)

        # Clean up memory
        del df_inherited_missense, df_inherited_lof, parent_lof, parent_missense, parent_total, df_parent, df_inherited
        gc.collect()
        return 

    #Function to parse loss of function mutations
    #mutation_type either Transmitted, Untransmitted
    def parseVEP_lof_output(mutation_type, dirname):
        lof_vep_output = dirname + "/UCSF_" + mutation_type + "_lof_vep_output.txt"
        all_function_py3.parse_VEP_output("/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/coding_analysis/UCSF_" + mutation_type  + "_lof_vep_output.vcf", lof_vep_output, cols = ["Consequence", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "HGVSc", "HGVSp", "Protein_position", "Amino_acids", "Codons", "ENSP", "SIFT", "PolyPhen", "CADD_PHRED", "CADD_RAW",  "gnomADe_AF"])
        UCSF_lof_preprocessed = pd.read_csv(dirname + "/UCSF_Transmitted_preprocessed_LOF.txt", sep = "\t")
        df_lof_vep_output = pd.read_csv(lof_vep_output, sep = "\t")
        merged_lof = pd.merge(UCSF_lof_preprocessed, df_lof_vep_output, on = ["#Uploaded_variation"], how = "left")
        UCSF_lof_path = dirname + "/UCSF_" + mutation_type + "_lof_vep_merged.txt"
        merged_lof.to_csv(UCSF_lof_path, sep = "\t", header = True, index = None)

    #Function to parse missense mutations
    #mutation_type could be Transmitted, Untransimitted, or DNV
    def parseVEP_missense_output(UCSF_missense_hg19_path, mutation_type, dirname):
        missense_hg19_raw = pd.read_csv(UCSF_missense_hg19_path, sep = "\t").sort_values(by = ["chr", "start_hg19"])
        missense_hg19 = missense_hg19_raw[missense_hg19_raw["start_hg19"]!= -1]
        missense_hg19["pos_hg19"] = missense_hg19.apply(lambda x: int(x["start_hg19"]) + 1, axis = 1)
        missense_hg19["vep"] = missense_hg19.sort_values(by = ["chr", "pos_hg19"]).apply(lambda x: " ".join([x["chr"][3:], str(x["pos_hg19"]), str(x["pos_hg19"] + len(x["ref"]) - 1), x["ref"] + "/" + x["alt"]]), axis = 1)
        vep_path = dirname + "/UCSF_" + mutation_type + "_missense_hg19_vep_input.txt"
        missense_hg19[["vep"]].drop_duplicates().to_csv(vep_path, sep = "\t", header = None, index = None)

        vep_output = dirname + "/UCSF_" + mutation_type + "_missense_hg19_vep_output.txt"
        all_function_py3.parse_VEP_output("/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/coding_analysis/UCSF_" + mutation_type + "_missense_hg19_vep_output.vcf", vep_output, cols = ["Consequence", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "HGVSc", "HGVSp", "Protein_position", "Amino_acids", "Codons", "ENSP", "SIFT", "PolyPhen", "CADD_PHRED", "CADD_RAW", "gnomADg_AF", "gnomADe_AF","MPC"])
        missense_hg19_output = pd.read_csv(vep_output, sep = "\t")
        missense_hg19_output = missense_hg19_output[~missense_hg19_output["MPC"].isna()]
        print("df before mpc: " + str(missense_hg19_output.shape[0]))

        missense_hg19_output["MPC"] = missense_hg19_output["MPC"].astype(float)
        if mutation_type == "DNV":
            missense_hg19_raw = missense_hg19_raw[["chr", "pos", "ref", "alt", "TYPE", "LEN", "FAM", "SampleID", "Pheno", "BATCH", "#Uploaded_variation", "start", "end", "start_hg19", "end_hg19"]]
        missense_hg19_raw["#Uploaded_variation_hg19"] = missense_hg19_raw.apply(lambda x: "_".join([x["chr"][3:],str(int(x["start_hg19"]) + 1), x["ref"] + "/" + x["alt"]]), axis = 1)
        
        return missense_hg19_output, missense_hg19_raw

    # Function to generate damaging missense table
    def getDmisTable(UCSF_Dmis, missense_hg19_raw, mis_col, gnomad_score, sex, dirname):
        UCSF_dmis_merged = pd.merge(UCSF_Dmis, missense_hg19_raw, on = ["#Uploaded_variation_hg19"])
        print("df Before gnoamd: " + str(UCSF_dmis_merged.shape[0]))

        UCSF_dmis_merged = pd.merge(UCSF_dmis_merged, sex[["SampleID", "sex"]].drop_duplicates(), on = "SampleID", how = "left")
        UCSF_dmis_filtered = gnomad_filter(UCSF_dmis_merged, gnomad_col = "gnomADe_AF", gnomAD = float(gnomad_score))
        print("df after gnoamd: " + str(UCSF_dmis_filtered.shape[0]))
        missense_coding_table = new_var_analysis(dirname, UCSF_dmis_filtered, var_type = mis_col, filter_status = "gnomad_" + str(gnomad_score), id = "Gene")
        missense_coding_table = missense_coding_table.rename(columns = {"Gene": "gene_id"})
        return missense_coding_table

    # Function to generate overall coding information table
    def getTotalTable(missense_hg19_path, UCSF_lof, mut2merge, gnomad_score, dirname, output_path, mutation):
        
        # missense_hg19_raw is UCSF variant file converted to hg19, with no further filtering/parsing
        missense_hg19_output, missense_hg19_raw = parseVEP_missense_output(missense_hg19_path, mutation_type = mutation, dirname = dirname)
        missense_hg19_output = missense_hg19_output.rename(columns = {"#Uploaded_variation" :"#Uploaded_variation_hg19"})

        UCSF_Dmisb = missense_hg19_output[missense_hg19_output["MPC"] >= float(2)]

        UCSF_Dmisa = missense_hg19_output[(missense_hg19_output["MPC"] < 2) & (missense_hg19_output["MPC"] >= 1)]
        misa_col = mutation + ".misa"
        misb_col = mutation + ".misb"
        lof_col = mutation + ".ptv"

        sex = pd.read_csv(args.sample_information, sep = "\t")


        misa_table = getDmisTable(UCSF_Dmisa, missense_hg19_raw, misa_col, gnomad_score, sex, dirname)
        misa_coding_annotated = pd.merge(misa_table, mut2merge, on = ["gene_id"], how = "right").fillna(0)
        misb_table = getDmisTable(UCSF_Dmisb, missense_hg19_raw, misb_col, gnomad_score, sex, dirname)            
        dmis_table = pd.merge(misb_table, misa_coding_annotated, on = ["gene_id"], how = "right").fillna(0)

        dmis_table.to_csv(dirname + "/UCSF_" + mutation +  "_gnomad_" + str(gnomad_score) + "_mutability_annotated_missense_table.txt", sep = "\t", header = True, index = None)

        UCSF_lof_filtered = gnomad_filter(UCSF_lof, gnomad_col = "gnomADe_AF", gnomAD = float(gnomad_score))
        UCSF_lof_merged = pd.merge(UCSF_lof_filtered, sex[["SampleID", "sex"]].drop_duplicates(), on = "SampleID", how = "left")
        lof_coding_table = new_var_analysis(dirname, UCSF_lof_merged, var_type = lof_col, filter_status = "gnomad_" + str(gnomad_score), id = "Gene")
        lof_coding_table = lof_coding_table.rename(columns = {"Gene": "gene_id"})
        lof_coding_annotated = pd.merge(lof_coding_table, mut2merge, on = ["gene_id"], how = "right").fillna(0)
        lof_coding_annotated.to_csv(dirname + "/UCSF_" + mutation + "_gnomad_" + str(gnomad_score) + "_mutability_annotated_lof_table.txt", sep = "\t", header = True, index = None)
        
        df_total = pd.merge(dmis_table, lof_coding_annotated, on = ["gene", "gene_gencodeV33", "gene_id", "chromosome", "pLI", "LOEUF", "mut.ptv", "mut.misb", "mut.misa"])
        new_order = ["gene", "gene_gencodeV33", "gene_id", "chromosome", "pLI", "LOEUF", "mut.ptv", "mut.misb", "mut.misa", lof_col, lof_col + ".sib", misb_col,  misb_col + ".sib", misa_col,  misa_col + ".sib"]
        df_total = df_total[new_order]

        df_total.to_csv(output_path + "UCSF_" + mutation + "_totalMutability_table.txt", sep = "\t", header = True, index = None)
        return df_total
    
 
    # Function to aggregate counts for each variant type
    def new_var_analysis(output_path, df, var_type, filter_status = "filtered", id = "Gene"):
        df_analysis_list = []
        print("created empty list")
        for gene, df_one in df.groupby(id):
            print("Entered for loop, on gene: " + str(gene))

            case = len(df_one[df_one["Pheno"] == "case"][["SampleID", "#Uploaded_variation"]].drop_duplicates())
            control = len(df_one[df_one["Pheno"] == "control"][["SampleID", "#Uploaded_variation"]].drop_duplicates())

            print("Calculated each count")
            # Append each group's result as a DataFrame to the list
            df_analysis_list.append(pd.DataFrame({id: [gene],
                                                str(var_type) : [case], 
                                                str(var_type) + '.sib': [control]}))
            
            print("appending counts to list")
        # Concatenate all group DataFrames
        df_analysis = pd.concat(df_analysis_list)
        print("Concated list")
        # column_sums = df_analysis.drop('Symbol', axis=1).sum()
        print("Calculated sum of each column")
        df_analysis.to_csv(output_path + "/UCSF_" + str(var_type) + "_" + str(filter_status) + "_analysis.txt", sep = "\t", header = True, index = None)
        return df_analysis

    # Function to process untransmitted data by chromosome
    def process_chromosome(i, df_parent_chr, df_annot, dirname):
        # Filter annotations for the specific chromosome
        chr_annot = df_annot[df_annot["#Uploaded_variation"].astype(str).str.startswith(str(i))]

        # Ensure there are annotations to process
        if chr_annot.empty:
            logging.warning(f"No annotations found for chromosome {i}")
            return

        # Call the preprocessing function for untransmitted variants
        Untransmitted_preprocess("chr" + str(i), df_parent_chr, chr_annot, dirname, coding_flag="damaging")

    DNV_output = dirname + "/DNV"
    if path.exists(DNV_output) == False:
        os.mkdir(DNV_output)

    transmitted_output = dirname + "/Transmitted"
    if path.exists(transmitted_output) == False:
        os.mkdir(transmitted_output)
    
    untransmitted_output = dirname + "/Untransmitted"
    if path.exists(untransmitted_output) == False:
        os.mkdir(untransmitted_output)

    UCSF_DNV_missense_hg19_path, UCSF_DNV_lof_path = DNV_preprocess(DNV_output)
    UCSF_Transmitted_missense_hg19_path, UCSF_Transmitted_lof_path = Transmitted_preprocess(transmitted_output, coding_flag = "damaging")
    # Load the parent data and separate by chromosome
    df_parent = pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Variants/UCSF_Parent.txt.gz", sep="\t", dtype=str)
    df_parent_by_chr = {i: df_parent[df_parent["chr"] == f"chr{i}"] for i in range(1, 23)}
    ######
    df_annot = pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Variants/UCSF_parent_var_annotation.txt", sep = "\t", dtype = str)
    df_annot["Variant_Classification"] = df_annot["Consequence"].apply(all_function_py3.consequence2varclass_coding)
    df_annot = df_annot[(df_annot["Variant_Classification"].isin(lof)) | (df_annot["Variant_Classification"] == "Missense_Mutation")]

    # # Parallel processing with joblib
    # Parallel(n_jobs=2)(delayed(process_chromosome)(i, df_parent_by_chr[i], df_annot, dirname) for i in range(1, 23))
   
   # Process each chromosome without parallelization
    for i in range(1, 23):
        process_chromosome(i, df_parent_by_chr[i], df_annot, dirname)

    coding_mut = pd.read_excel(args.Mutability_file,
                          sheet_name="Supplementary Table 5", dtype=str)
    coding_mut = coding_mut.dropna(subset=['gene_id'])
    mut2merge = coding_mut[["gene", "gene_gencodeV33", "gene_id", "chromosome", "pLI", "LOEUF", "mut.ptv", "mut.misb", "mut.misa"]]

    # Initialize dictionaries to store the data tables for each mutation
    tables = {}

    for mutation in ["DNV", "Transmitted"]:
        # Construct file paths for missense and loss of function (lof) data
        missense_hg19_path = locals()[f"UCSF_{mutation}_missense_hg19_path"]
        lof_path = locals()[f"UCSF_{mutation}_lof_path"]
        UCSF_lof = pd.read_csv(lof_path, sep="\t")
        dir_path = f"{dirname}/{mutation}"
        for gnomad_score in [0]:
            # Call function to get the data table for this mutation and score
            table = getTotalTable(missense_hg19_path, UCSF_lof, mut2merge, gnomad_score, dir_path, output_path, mutation)

            # Store the resulting table in the dictionary
            tables[mutation] = table
    
    # Merge the tables for "DNV" and "Transmitted"
    total_table = pd.merge(tables["DNV"], tables["Transmitted"], on=["gene", "gene_gencodeV33", "gene_id", "chromosome", "pLI", "LOEUF", "mut.ptv", "mut.misb", "mut.misa"])
    total_table.to_csv(output_path + "Total_CodingAnalysis_Table.txt", sep = "\t", header = True, index = None)
