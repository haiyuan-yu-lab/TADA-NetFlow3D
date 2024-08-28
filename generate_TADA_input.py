description = "Python Script to link variants to genes"

import pandas as pd
import numpy as np
import argparse
import sys
import importlib
sys.path.append("/home/yz2296/bin")
import all_function_py3
import os
from os import path
import json
from collections import defaultdict
import scipy.stats
import logging
import warnings

# dirname = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
	# user input
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-A','--All_gene_file', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Promoters/promoters_1kb_tss_centered.coding.v43.bed", help = 'Path to all gene file.')
    parser.add_argument('-D','--DNV_Mutation_file', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/CHTN13_CHTN22_HN21_LUHMES.Li_SC2_SC3/promoter.v43.coding/UCSF_DNV_ExcludeCoding_damaging_enhancer.txt", help = 'Path to DNV mutation file.')
    parser.add_argument('-R','--Raw_DNV_file', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/UCSF_DNV_ExcludeCoding_damaging.txt", help = 'Path to DNV mutation file.')

    parser.add_argument('-T','--Transmitted_Mutation_file', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Variants/UCSF_Transmitted_ExcludeCoding_damaging.txt", help = 'Path to Transmitted mutation file.')
    parser.add_argument('-P','--Transmitted_Parent_file', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/CHTN13_CHTN22_HN21_LUHMES.Li_SC2_SC3/promoter.v43.coding/UCSF_Parent_enhancer.txt", help = 'Path to Transmitted Parent enhancer file.')
    parser.add_argument('-s','--Samocha_Mutability_file', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Mutability/20230919_samocha14_calcBaseLineRate_EG_mapping_allsample.bed", help = 'Path to Samocha Mutability file.')
    parser.add_argument('-o','--output_path', type = str, default = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/", help = '[OPTIONAL] Path to the job folder. If not specified, the job folder will be stored in the default output folder')
    parser.add_argument('-I','--job_name', required = True, type = str, help = 'Please specify a name for your job. The output files will be stored with this name')

    # parser.add_argument('-I','--job_name', required = True, type = str, help = 'Please specify a name for your job. The output files will be stored with this name')
    args = parser.parse_args()

############################################################################################################################################
	# Check input files
	# -----------------
    if os.path.isdir(args.output_path) == False:
        sys.exit("No such directory: " + args.output_path)
    if args.output_path[-1] == "/":
        output_path = args.output_path + args.job_name + "/"
    else:
        output_path = args.output_path + "/" + args.job_name + "/"
    if path.exists(output_path) == False:
        os.mkdir(output_path)

    if path.exists(args.DNV_Mutation_file) == False:
        sys.exit("No such directory: " + args.DNV_Mutation_file)
    if path.exists(args.Transmitted_Mutation_file) == False:
        sys.exit("No such directory: " + args.Transmitted_Mutation_file)
    if path.exists(args.Transmitted_Parent_file) == False:
        sys.exit("No such directory: " + args.Transmitted_Parent_file)
    if path.exists(args.Samocha_Mutability_file) == False:
        sys.exit("No such directory: " + args.Samocha_Mutability_file)
    

	############################################################################################################################################


    #Function to filter gnomad cutoff
    def parse_NC(df_NC, gnomAD):
        ############################################################################################################################################
		# df_NC = df_NC[df_NC["gnomADg_AF"] != "-"]
        ############################################################################################################################################
        df_NC_in_gnomad = df_NC[df_NC["gnomADg_AF"] != "-"]
        df_NC_notin_gnomad = df_NC[df_NC["gnomADg_AF"] == "-"]
        
        df_NC_in_gnomad["gnomADg_AF"] = df_NC_in_gnomad["gnomADg_AF"].astype(float)
        df_NC_in_gnomad = df_NC_in_gnomad[df_NC_in_gnomad["gnomADg_AF"] <= float(gnomAD)]
        df_NC_filtered = pd.concat([df_NC_in_gnomad, df_NC_notin_gnomad])
        return df_NC_filtered

    #Function to filter damaging variants
    def damaging_filter(NC, CADD, phyloP, gnomAD):
        df_NC = parse_NC(NC, gnomAD)
        NC_SNV = df_NC[df_NC["TYPE"] == "SNV"]
        NC_SNV["CADD_PHRED"] = NC_SNV["CADD_PHRED"].astype(float)
        NC_SNV = NC_SNV[NC_SNV["CADD_PHRED"] >= float(CADD)]

        NC_Indel = df_NC[df_NC["TYPE"] == "Indel"]
        NC_Indel["phyloP100way"] = NC_Indel["phyloP100way"].astype(float)
        NC_Indel = NC_Indel[NC_Indel["phyloP100way"] >= float(phyloP)]
        return pd.concat([NC_SNV, NC_Indel])

    #Function to combine all EG linkage from every cell type
    #Converts EG linkage from hg19 to hg38 
    def get_allEG(output_path):
        label2sample = {"Cerebral_cortex": {"CHTN22", "HN21"}, "Cerebellum": {"CHTN13"}, "Neuronal Progenitor Cells": {"LUHMES.Li"}, "Astrocyte": {"SC3"}, "Excitatory Cortical Neurons": {"SC2"}} #, "Combined": {"CHTN22", "HN21", "CHTN13", "LUHMES.Li", "SC3", "SC2"}, "K562": {"K562"}, "HEK293": {"HEK293"}}
        # Combining all the enhancer hg38 files
        total_EG = pd.DataFrame()
        for types in label2sample.keys():
            print("In cell type: " + str(types))
            enh_path = "/fs/cbsuhy01/storage/NetFlow_ASD/EG_linkage/Predictions/" + str(types) + "/pints/0.02/EnhancerPredictions.txt"
            if path.exists(enh_path) == True:
                logging.warning("Using new EG for cell type: " + str(types))
                new_EG = pd.read_csv(enh_path, sep = "\t")
                new_EG["Overlapped_elements_enhancer"] = new_EG["name"].apply(lambda x: x.split("|")[-1])
                new_EG = new_EG.rename(columns = {"TargetGene" : "Symbol", "start" : "start_hg38", "end" : "end_hg38", "chr" : "chr_enhancer"})
                df_EG_merge = pd.concat([total_EG, new_EG])
            else:
                logging.warning("Using old EG for cell type: " + str(types))
                all_conv = pd.DataFrame()
                # Combining all the EG linkage files
                all_enh = pd.DataFrame()
                for cell in label2sample[types]:
                    logging.warning("Enter cell: " + str(cell))
                    all_conv = pd.concat([all_conv, pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/" + cell + "/promoter.v43.coding/all_peaks_element_200bp_e_hg19.bed", sep ="\t")], ignore_index=True)                        
                    
                    # ["CHTN22", "HN21", "CHTN13", "SC2", "SC3", "LUHMES.Li"]
                    ##"CHTN13" --> naming
                    # for cell in label2sample[types]:
                    for x in ["CP", "GZ"]:
                        all_enh = pd.concat([all_enh, pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/" + cell + "/promoter.v43.coding/EG_connections/all_peaks_element_200bp_e_hg19_NoHeader__promoters_1kb_tss_centered.coding.v43.hg19.NoHeader__merged_interaction_profile_" + x + ".txt", sep ="\t")], ignore_index=True)
                    all_enh_clean = all_enh.drop_duplicates().reset_index(drop=True)

                    #Generating df expanding all enhancers
                    # df = all_enh_clean[all_enh_clean["chr"].isin(["chr" + str(i) for i in range(1, 23)])]

                    intermediate_dfs = []  # List to hold intermediate DataFrames

                    for gene, df_one in all_enh_clean.groupby("7"):
                        # Check if all values in the 6th column are the same
                        column_6 = df_one.iloc[:, 6]
                        symbol = column_6.iloc[0] if column_6.nunique() == 1 else None

                        enhancers = []
                        for col in ["Distal_enhancers", "Proximal_enhancers"]:
                            for row in df_one[col]:
                                if type(row) == str:
                                    enhancers.extend([":".join(x.split(":")[0:2]) for x in row.split(",")])
                        enhancers = sorted(set(enhancers))

                        # Create an intermediate DataFrame for each element
                        for element in enhancers:
                            parts = element.split(':')
                            chrom = parts[0]
                            start, end = parts[1].split('-')
                            intermediate_dfs.append(pd.DataFrame({'gene_id': [gene], 'Symbol': [symbol], 'enhancer_hg19': [element], 'chr': [chrom], 'start_hg19' : [start], 'end_hg19' : [end]}))

                    # Concatenate all intermediate DataFrames to form the final DataFrame
                    df_EG_final = pd.concat(intermediate_dfs, ignore_index=True)
                    df_EG_final["CellType"] = str(types)
                    df_EG_final["start_hg19"] = df_EG_final["start_hg19"].astype(int)
                    df_EG_final["end_hg19"] = df_EG_final["end_hg19"].astype(int)
                    df_EG_final = df_EG_final.loc[:,~df_EG_final.columns.duplicated()]
                    # or for all_conv
                    all_conv = all_conv.loc[:,~all_conv.columns.duplicated()]
                    all_conv.rename(columns={'start': 'start_hg19', 'end': 'end_hg19'}, inplace=True)

                    #Getting hg38 info
                    df_EG_merge = pd.merge(df_EG_final, all_conv, on=["chr", "start_hg19", "end_hg19"], how='left')
                    #Generating column of unique enhancers with hg38 
                    df_EG_merge["Overlapped_elements_enhancer"] = df_EG_merge.apply(lambda row : str(row['chr']) +  ":" + str(row['start_hg38']) + "-" + str(row['end_hg38']), axis=1)
                    df_EG_merge = df_EG_merge.rename(columns = {'chr': 'chr_enhancer'})
            
            total_EG = pd.concat([total_EG, df_EG_merge])
                
        total_EG.to_csv(output_path + "all_EG_linkage_hg38.bed", sep = "\t", header = True, index = None)
        return total_EG

    def get_oldEG(output_path):
        # label2sample = {"Cerebral_cortex": {"CHTN22", "HN21"}, "Cerebellum": {"CHTN13"}, "Neuronal Progenitor Cells": {"LUHMES.Li"}, "Astrocyte": {"SC3"}, "Excitatory Cortical Neurons": {"SC2"}} #, "Combined": {"CHTN22", "HN21", "CHTN13", "LUHMES.Li", "SC3", "SC2"}, "K562": {"K562"}, "HEK293": {"HEK293"}}
        # Combining all the enhancer hg38 files
        
        all_conv = pd.DataFrame()
        all_enh = pd.DataFrame()
        for cell in ["CHTN22", "HN21", "SC2", "SC3", "LUHMES.Li"]:
            logging.warning("Enter cell: " + str(cell))
            all_conv = pd.concat([all_conv, pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/" + cell + "/promoter.v43.coding/all_peaks_element_200bp_e_hg19.bed", sep ="\t")], ignore_index=True)
                
            # ["CHTN22", "HN21", "CHTN13", "SC2", "SC3", "LUHMES.Li"]
            ##"CHTN13" --> naming
            for x in ["CP", "GZ"]:
                all_enh = pd.concat([all_enh, pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/" + cell + "/promoter.v43.coding/EG_connections/all_peaks_element_200bp_e_hg19_NoHeader__promoters_1kb_tss_centered.coding.v43.hg19.NoHeader__merged_interaction_profile_" + x + ".txt", sep ="\t")], ignore_index=True)
            all_enh_clean = all_enh.drop_duplicates().reset_index(drop=True)

            #Generating df expanding all enhancers
            # df = all_enh_clean[all_enh_clean["chr"].isin(["chr" + str(i) for i in range(1, 23)])]
        all_conv.rename(columns={'start': 'start_hg19', 'end': 'end_hg19'}, inplace=True)

        intermediate_dfs = []  # List to hold intermediate DataFrames

        for gene, df_one in all_enh_clean.groupby("7"):
            # Check if all values in the 6th column are the same
            column_6 = df_one.iloc[:, 6]
            symbol = column_6.iloc[0] if column_6.nunique() == 1 else None

            enhancers = []
            for col in ["Distal_enhancers", "Proximal_enhancers"]:
                for row in df_one[col]:
                    if type(row) == str:
                        enhancers.extend([":".join(x.split(":")[0:2]) for x in row.split(",")])
            enhancers = sorted(set(enhancers))

            # Create an intermediate DataFrame for each element
            for element in enhancers:
                parts = element.split(':')
                chrom = parts[0]
                start, end = parts[1].split('-')
                intermediate_dfs.append(pd.DataFrame({'gene_id': [gene], 'Symbol': [symbol], 'enhancer_hg19': [element], 'chr': [chrom], 'start_hg19' : [start], 'end_hg19' : [end]}))

        # Concatenate all intermediate DataFrames to form the final DataFrame
        df_EG_final = pd.concat(intermediate_dfs, ignore_index=True)
        df_EG_final["CellType"] = str(cell)
        df_EG_final["start_hg19"] = df_EG_final["start_hg19"].astype(int)
        df_EG_final["end_hg19"] = df_EG_final["end_hg19"].astype(int)
        df_EG_final = df_EG_final.loc[:,~df_EG_final.columns.duplicated()]
        # or for all_conv
        all_conv = all_conv.loc[:,~all_conv.columns.duplicated()]
        #Getting hg38 info
        df_EG_merge = pd.merge(df_EG_final, all_conv, on=["chr", "start_hg19", "end_hg19"], how='left')
        #Generating column of unique enhancers with hg38 
        df_EG_merge["Overlapped_elements_enhancer"] = df_EG_merge.apply(lambda row : str(row['chr']) +  ":" + str(row['start_hg38']) + "-" + str(row['end_hg38']), axis=1)
        df_EG_merge = df_EG_merge.rename(columns = {'chr': 'chr_enhancer'})
        
        # total_EG = pd.concat([total_EG, df_EG_merge])
                
        df_EG_merge.to_csv(output_path + "old_all_EG_linkage_hg38.bed", sep = "\t", header = True, index = None)
        return df_EG_merge

    #Function to filter DNV variants
    def parse_DNV(df_EG_merge, CADD, phyloP, gnomad):
        # CADD = 12
        # phyloP = 2
        # gnomad = 0
        #Read in mutants file and merge with EG
        df_mut2enh = pd.read_csv(args.DNV_Mutation_file, sep = "\t")
        
        # df_mut2enh = pd.read_csv("/local/storage/NetFlow_ASD/Week8_429/CHTN13_CHTN22_HN21_LUHMES.Li_SC2_SC3/promoter.v43.coding/UCSF_DNV_ExcludeCoding_damaging_enhancer.txt", sep ="\t")
        df_mut2gn = pd.merge(df_mut2enh, df_EG_merge, on="Overlapped_elements_enhancer", how='inner')
        df_mut2gn_filterd = damaging_filter(df_mut2gn, CADD, phyloP, gnomad)

        return df_mut2gn_filterd
        

    #Function to filter transmitted variants
    def parse_Transmitted(df_EG_merge, CADD, phyloP, gnomAD):
        # CADD = 12
        # phyloP = 2
        # gnomAD = 0
    
        df_parent = pd.read_csv(args.Transmitted_Parent_file, sep = "\t", dtype = str, usecols=["chr", "#Uploaded_variation", "Overlapped_elements_enhancer"])
        # df_parent = df_parent.drop(["pos", "ref", "alt", "SampleID"], axis=1)

        logging.warning("read in parent")
        df_transmitted_unfiltered = pd.read_csv(args.Transmitted_Mutation_file, sep = "\t", dtype = str, usecols = ["chr", "#Uploaded_variation", "SampleID", "CADD_PHRED", "phyloP100way", "gnomADg_AF", "sex", "Pheno", "TYPE"])
        logging.warning("read in transmitted")
        df_transmitted_filtered = damaging_filter(df_transmitted_unfiltered, CADD, phyloP, gnomAD)
        transmitted_scaling_factor = getScalingFactor(df_transmitted_unfiltered, CADD, phyloP, gnomAD)

        merged_transmitted_filtered = pd.DataFrame()
        # merged_transmitted_unfiltered = pd.DataFrame()
        for i in range(1, 23):
            logging.warning("Enter for loop for chr: " + str(i))
            parent = df_parent[df_parent["chr"].isin(["chr" + str(i)])]
            EG_merge = df_EG_merge[df_EG_merge["chr_enhancer"].isin(["chr" + str(i)])]
            # Merge filtered variants with parent data
            transmitted_enh2gn = pd.merge(parent, EG_merge, on=["Overlapped_elements_enhancer"], how='inner')
            logging.warning("Merged parent with EG linkage")
            df_var_filtered = df_transmitted_filtered[df_transmitted_filtered["chr"].isin(["chr" + str(i)])]
            df_transmitted_filtered_gn2var = pd.merge(df_var_filtered, transmitted_enh2gn, on=["#Uploaded_variation"], how='inner')
            merged_transmitted_filtered = pd.concat([merged_transmitted_filtered, df_transmitted_filtered_gn2var])

        # merged_transmitted_filtered = damaging_filter(merged_transmitted_unfiltered, CADD, phyloP, gnomAD)
        merged_transmitted_filtered.to_csv(output_path + "EGmerged_transmitted_filtered.txt", sep = "\t", header = True, index = None)
        return merged_transmitted_filtered, transmitted_scaling_factor #, merged_transmitted_unfiltered 

   

    #Function to analyze variants in both transmitted and DNV based on pheno and sex
    #type is either transmitted or DNV
    def var_analysis(output_path, df_NC, var_type, filter_status = "filtered"):
        df_analysis_list = []
        print("created empty list")
        for gene, df_one in df_NC.groupby("Symbol"):
            print("Entered for loop, on gene: " + str(gene))

            #check if there is no EG linkage
            if (~df_one["enhancer_regions"].isna()).sum() == 0:
                case = control = male_case = female_case = male_control = female_control = np.nan
            
            #check if there is EG linkage but no mapped variants
            elif ((~df_one["enhancer_regions"].isna()).sum() > 0) & ((~df_one["Overlapped_elements_enhancer"].isna()).sum() == 0):
                case = control = male_case = female_case = male_control = female_control = 0
            
            #otherwise, EG exists and linked to variants
            else:
                case = len(df_one[df_one["Pheno"] == "case"]["SampleID"].unique())
                control = len(df_one[df_one["Pheno"] == "control"]["SampleID"].unique())
                male_case = len(df_one[(df_one["sex"] == "male") & (df_one["Pheno"] == "case")]["SampleID"].unique())
                female_case = len(df_one[(df_one["sex"] == "female") & (df_one["Pheno"] == "case")]["SampleID"].unique())
                male_control = len(df_one[(df_one["sex"] == "male") & (df_one["Pheno"] == "control")]["SampleID"].unique())
                female_control = len(df_one[(df_one["sex"] == "female") & (df_one["Pheno"] == "control")]["SampleID"].unique())

            print("Calculated each count")
            # Append each group's result as a DataFrame to the list
            df_analysis_list.append(pd.DataFrame({'Symbol': [gene], 
                                                'Noncoding_' + str(var_type) + '_all_case': [case], 
                                                'Noncoding_' + str(var_type) + '_all_control': [control], 
                                                'Noncoding_' + str(var_type) + '_male_case' : [male_case], 
                                                'Noncoding_' + str(var_type) + '_male_control' : [male_control], 
                                                'Noncoding_' + str(var_type) + '_female_case' : [female_case], 
                                                'Noncoding_' + str(var_type) + '_female_control': [female_control]}))
            # except:
            #     logging.info(gene)
            print("appending counts to list")
        # Concatenate all group DataFrames
        df_analysis = pd.concat(df_analysis_list)
        print("Concated list")
        # column_sums = df_analysis.drop('Symbol', axis=1).sum()
        print("Calculated sum of each column")
        df_analysis.to_csv(output_path + "new_UCSF_" + str(var_type) + "_" + str(filter_status) + "_analysis.txt", sep = "\t", header = True, index = None)
        return df_analysis #, column_sums

   
    #Function to combine all enhancers to corresponding gene
    #Uses bed tool
    def get_enhancer2gene(df_EG_merge):
        df_final = defaultdict(list)

        for gene,df_one in df_EG_merge.groupby("Symbol"):
            df_final["Symbol"].append(gene)
            elements = df_one["Overlapped_elements_enhancer"].unique().tolist()
            elements = sorted(set(elements))
            if len(elements) > 0:
                df_tmp = pd.DataFrame({"elements": elements})
                df_tmp["chr"] = df_tmp["elements"].apply(lambda x: x.split(":")[0])
                df_tmp["start_hg38"] = df_tmp["elements"].apply(lambda x: x.split(":")[1].split("-")[0])
                df_tmp["end_hg38"] = df_tmp["elements"].apply(lambda x: x.split(":")[1].split("-")[1])
        #         df_tmp = pd.merge(df_tmp, df_conversion, on = ["chr", "start", "end"])
                df_tmp["ID"] = df_tmp.apply(lambda x: x["chr"] + ":" + x["start_hg38"] + "-" + x["end_hg38"], axis = 1)
                df_tmp[["chr", "start_hg38", "end_hg38"]].to_csv(gene + ".bed", sep = "\t", header = None, index = None)
                os.system("/programs/bin/bedtools/bin/bedtools sort -i " + gene + ".bed | /programs/bin/bedtools/bin/bedtools merge -i - > " + gene + "_merged.bed")
                df_merged = pd.read_csv(gene + "_merged.bed", sep = "\t", header = None, dtype = str)
                df_merged["region"] = df_merged.apply(lambda x: x[1] + "-" + x[2], axis = 1)
                df_final["enhancer_regions"].append(df_merged[0].tolist()[0] + ":" + ",".join(df_merged["region"].tolist()))
                os.remove(gene + ".bed")
                os.remove(gene + "_merged.bed")
            else:
                df_final["enhancer_regions"].append(np.nan)
                
        df_final = pd.DataFrame(df_final)

        df_final.to_csv(output_path + "gene2enhancer_regions_hg38.bed", sep = "\t", header = True, index = None)
        return df_final

    # def old_getMutability(dirname, df_EG_merge, df_analysis):
    #     baseline_mut = pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Mutability/20230919_samocha14_calcBaseLineRate_EG_mapping_allsample.bed", sep = "\t")
    #     baseline_mut = baseline_mut.rename(columns = {"chrom": "chr_enhancer", "start": "start_hg38", "end" : "end_hg38"})
    #     df_enh_mut_merge = pd.merge(df_EG_merge, baseline_mut, on = ["chr_enhancer", "start_hg38", "end_hg38"], how = "left")

    #     df_enhmut = []

    #     for gene,df_one in df_enh_mut_merge.groupby("gene_id"):
    #     #     df_enhmut["gene_id"].append(gene)
    #         print("current gene is: " + gene)
    #     #     enh_start_min = min(df_one["start_hg38"])
    #     #     enh_start_max = max(df_one["start_hg38"])
    #     #     enh_end_min = min(df_one["end_hg38"])
    #     #     enh_end_max = max(df_one["end_hg38"])
    #     #     enhancers = []
    #         # Sort the DataFrame by the start column
    #         df = df_one.sort_values(by='start_hg38').reset_index(drop=True)
    #     #     for elements in df_one["Overlapped_elements_enhancer"].unique():
    #     #         enhancers.append(elements)
    #         if len(df_one["Overlapped_elements_enhancer"].unique().tolist())>= 1:
        
    #             # Initialize a list to hold the merged ranges
    #             merged_ranges = []
                
    #             if len(df_one["chr_enhancer"].unique()) == 1:
    #                 chrom = df.iloc[0]['chr_enhancer']
    #                 # Initialize the first range
    #                 current_start, current_end = df.iloc[0]['start_hg38'], df.iloc[0]['end_hg38']
    #                 current_size = 0
                    
    #                 max_size = df_one["size"].max()
    #                 max_size_row = df_one[df_one["size"] == max_size]
    #                 mutability_bp = max_size_row.iloc[0]["baselinerate_samocha14"] / max_size
                    
    #                 # Iterate over the enhancer regions
    #                 for i in range(1, len(df)):
    #                     start, end = df.iloc[i]['start_hg38'], df.iloc[i]['end_hg38']
    #                     if start <= current_end:  # Overlap detected
    #                         print('Found overlap')
    #                         print("current size is: " + str(current_size))
    #                         print("old enhancer is: " + str(chrom) + ":" + str(current_start) + "-" + str(current_end))
    #                         print("currnet enhancer is: " + str(chrom) + ":" + str(start) + "-" + str(end))

    #                         current_end = max(current_end, end)  # Extend the current range
    #                         current_size = current_size + (current_end-current_start)

    #                         print("new size is: " + str(current_size))

    #                     else:
    #                         # No overlap, add the current range to the list and start a new range
    #                         current_start, current_end = start, end
    #                         current_size = current_size + (current_end - current_start)
                        
    #                         print("current enhancer is: " + str(chrom) + ":" + str(current_start) + "-" + str(current_end))
    #                         print("current size is: " + str(current_size))
    #     #                     merged_ranges.append((gene, chrom, current_start, current_end))
                
    #                 mutability = mutability_bp * current_size
                

    #             # Add the last range
    #             merged_ranges.append((gene, current_size, mutability_bp, mutability))

    #             # Convert the merged ranges into a DataFrame if needed
    #             df_enhmut.append(pd.DataFrame(merged_ranges, columns=['gene_id', 'total_enhancer_size', 'Mutability_BP_Samocha', 'Mutability_Samocha_calculated']))

    #     df_mutability = pd.concat(df_enhmut)

    #     df_analysis_mut = pd.merge(df_analysis, df_mutability, on = ["gene_id"], how = "left")
    #     return df_analysis_mut



    # def mutability_scaling(filtered_sum, unfiltered_sum):
    #     scaled_values = {}  # Use a dictionary to collect scaled values
    #     list_of_keys = list(filtered_sum.keys())
    #     for key in list_of_keys:
    #         # Calculate scaled value and store it in the dictionary
    #         scaled_values[str(key)] = filtered_sum[key] / unfiltered_sum[key]
        
    #     # Create a DataFrame from the dictionary, which automatically aligns keys as columns
    #     # The values are put into a single row in the DataFrame
    # #     df_scaling = pd.DataFrame([scaled_values])
    #     return scaled_values

    #Function to calculate and annotate Samocha mutability for each gene
    def get_SamochaMutability(baseline_mut, enhancers2gn, output_path, scaling_factor, variant_type = "DNV"):
        # Convert baseline_mut for faster lookup
        # df_enh_mut_merge = pd.merge(df_EG_merge, baseline_mut, on = ["chr_enhancer", "start_hg38", "end_hg38"], how = "left")

        baseline_mut_dict = baseline_mut.set_index("Overlapped_elements_enhancer")["baselinerate_samocha14"].to_dict()

        mutability = []

        mutability_scaled = []

        for enhancers in enhancers2gn["enhancer_regions"]:
            print(enhancers)
            if str(enhancers) == "nan":
                mutability.append(np.nan)
                mutability_scaled.append(np.nan)

            else:  
                print(str(enhancers))
                chrom, regions = str(enhancers).split(":")
                total_size = 0
                max_size = 0
                max_rate = 0

                for region in regions.split(","):
                    start, end = map(int, region.split("-"))
                    current_size = end - start
                    enhancer = f"{chrom}:{region}"

                    # Lookup in pre-processed dictionary
                    if enhancer in baseline_mut_dict:
                        print("enhancer in samocha")
                        total_size += current_size
                        if current_size > max_size:
                            max_size = current_size
                            max_rate = baseline_mut_dict[enhancer]

                # Avoid division by zero
                if max_size > 0:
                    mutability.append((max_rate / max_size) * total_size)
                    mutability_scaled.append((max_rate / max_size) * total_size * scaling_factor)
                else:
                    mutability.append(0)
                    mutability_scaled.append(0)

        enhancers2gn["Noncoding_Mutability_Samocha"] = mutability
        enhancers2gn["Noncoding_Scaled_Mutability_Samocha"] = mutability_scaled

        enhancers2gn.to_csv(output_path + str(variant_type) + "_gene_annotated_samocha_mutability.txt", sep = "\t", header = True, index = None)
        return enhancers2gn

    def get_RouletteMutability(enhancer2gn, scaling_factor, output_path):
        # dirname = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429"
        # output_path = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis"

        element2file = {"enhancer": "all_peaks_element_200bp_e.bed", "promoter": "all_peaks_element_200bp_p.bed"}
        sex2label = {"all": {"male", "female", np.nan}, "male": {"male"}, "female": {"female"}}
        label2sample = {"Cerebral Cortex": {"CHTN22", "HN21"}, "Cerebellum": {"CHTN13"}, "Neuronal Progenitor Cells": {"LUHMES.Li"}, "Astrocytes": {"SC3"}, "Excitatory Cortical Neurons": {"SC2"}, "Combined": {"CHTN22", "HN21", "CHTN13", "LUHMES.Li", "SC3", "SC2"}, "K562": {"K562"}, "HEK293": {"HEK293"}}

        for score in ["MR", "MC"]:
            enhancer2gn["Noncoding_Mutability_" + score] = enhancer2gn["enhancer_regions"].apply(lambda x: get_gene_mutability_Yingying(x, score))
            enhancer2gn["Noncoding_Scaled_Mutability_" + score] = enhancer2gn["Noncoding_Mutability_" + score].apply(lambda x: float(x)*scaling_factor)
        enhancer2gn[["Symbol"] + [x for x in enhancer2gn.columns if x.find("Mutability") != -1] + [x for x in enhancer2gn.columns if x.find("Mutability") == -1 and x != "Symbol"]].to_csv(output_path + "/Roulette_final.txt", sep = "\t", header = True, index = None)
        return enhancer2gn

    def get_gene_mutability_Yingying(enhancers, score, data_path = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/TADA_2022-main/"):
        if type(enhancers) == str:
            chrom, regions = enhancers.split(":")
    #         print(chrom)
            if chrom[3:].isdigit() and 1 <= int(chrom[3:]) <= 22:
                print(chrom)
                with open(data_path + score + "_" + chrom + ".json") as f:
                    mut_rate = json.load(f)
    #                 print(mut_rate)
                output = 0
                for region in regions.split(","):
                    if chrom + ":" + region in mut_rate:
                        output = output + mut_rate[chrom + ":" + region]
    #                     print(output)
                return output
            else:
                return np.nan

    def getScalingFactor(df_var, CADD, phyloP, gnomAD):
        unfiltered_var = df_var[["#Uploaded_variation", "SampleID"]].drop_duplicates().shape[0]
        logging.warning("unfiltered variants is: " + str(unfiltered_var))


        filtered_var = damaging_filter(df_var, CADD, phyloP, gnomAD)
        filtered_var_count = filtered_var[["#Uploaded_variation", "SampleID"]].drop_duplicates().shape[0]
        logging.warning("filtered variants is: " + str(filtered_var_count))


        scaling_factor = filtered_var_count/unfiltered_var
        return scaling_factor

    def getTotalMutabilityTable(total_EG_merge, total_gn2enh, scaling_factor, variant_type = "DNV"):
        #Calculate Samocha Mutability
        baseline_mut = pd.read_csv(args.Samocha_Mutability_file, sep = "\t")
        # baseline_mut = pd.read_csv("/local/storage/NetFlow_ASD/Mutability/20230919_samocha14_calcBaseLineRate_EG_mapping_allsample.bed", sep = "\t")
        baseline_mut = baseline_mut.rename(columns = {"chrom": "chr_enhancer", "start": "start_hg38", "end" : "end_hg38"})
        samocha2gn = pd.merge(total_EG_merge, baseline_mut, on = ["chr_enhancer", "start_hg38", "end_hg38"], how = "left")

        df_samocha = get_SamochaMutability(samocha2gn, total_gn2enh, output_path, scaling_factor, variant_type)
        
        #Calculate Roulette Mutability
        # df_roulette_samocha = get_RouletteMutability(df_samocha, scaling_factor, output_path, variant_type)
        # df_roulette_samocha.to_csv(output_path + str(variant_type) + "gene_annotated_total_mutability.txt", sep = "\t", header = True, index = None)

        return df_samocha #, df_roulette_samocha


    ##########################
    # df_EG_merge = pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/all_EG_linkage_hg38.bed", sep = "\t")
    # df_DNV_analysis = pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/new_UCSF_DNV_filtered_analysis.txt", sep = "\t")
    
    ##########################
    #SETTINGS
    CADD = 12
    phyloP = 2
    gnomAD = 0
    ##########################


    #Get all EG linkage
    # total_EG_merge = get_allEG(output_path)

    total_EG_merge = get_oldEG(output_path)

    #Get dataframe with unique gene to all its corresponding enhancers
    df_enh2gn = get_enhancer2gene(total_EG_merge)

    #Get list of all genes and symbols with corresponding enhancer 
    #NaN if no EG linkage
    all_gene = pd.read_csv(args.All_gene_file, sep = "\t", header = None)
    all_gene = all_gene.rename(columns = {6 : "Symbol", 7 : "gene_id"})
    all_gene = all_gene.drop_duplicates(subset = ["Symbol", "gene_id"])
    total_gn2enh = pd.merge(all_gene, df_enh2gn, on = "Symbol", how = "left")
    total_gn2enh = total_gn2enh.drop([0, 1, 2, 3, 4, 5, 8, 9], axis = 1)
    
    #Filter DNV variants
    df_DNV_filtered= parse_DNV(total_EG_merge, CADD, phyloP, gnomAD)
   
    #Get all EG annotated with DNV variants
    total_gn2enh2DNV = pd.merge(total_gn2enh, df_DNV_filtered, on = ["Symbol", "gene_id"], how = "left")

    #Get analysis for DNV
    df_DNV_analysis = var_analysis(output_path, total_gn2enh2DNV, var_type = "DNV", filter_status = "filtered")

    #Calculate DNV scaling factor
    raw_DNV_var = pd.read_csv(args.Raw_DNV_file, sep = "\t", dtype = str)
    DNV_scaling_factor = getScalingFactor(raw_DNV_var, CADD, phyloP, gnomAD)
    
    df_DNV_mutability = getTotalMutabilityTable(total_EG_merge, total_gn2enh, DNV_scaling_factor, variant_type = "DNV")
    df_DNV_total = pd.merge(df_DNV_analysis, df_DNV_mutability, on = "Symbol", how = "left")
    df_DNV_total.to_csv(output_path + "DNV_Mutability_total.txt",  sep = "\t", header = True, index = None)
   
    ###########################
    #UNCOMMENT
    # Filter transmitted variants
    df_transmitted_filtered, transmitted_scaling_factor = parse_Transmitted(total_EG_merge, CADD, phyloP, gnomAD)
    logging.warning("Returned from parse_transmitted")
   
    #Get all EG annotated with transmitted variants
    total_gn2enh2transmitted = pd.merge(total_gn2enh, df_transmitted_filtered, on = ["Symbol", "gene_id"], how = "left")
    
    # Get analysis for transmitted 
    df_transmitted_analysis = var_analysis(output_path, total_gn2enh2transmitted, var_type = "transmitted", filter_status = "filtered")

    df_transmitted_mutability = getTotalMutabilityTable(total_EG_merge, total_gn2enh, transmitted_scaling_factor, variant_type = "transmitted")
    df_transmitted_total = pd.merge(df_transmitted_analysis, df_transmitted_mutability, on = "Symbol", how = "left")
    df_transmitted_total.to_csv(output_path + "transmitted_Mutability_total.txt",  sep = "\t", header = True, index = None)
  
    # df_transmitted_filtered, df_transmitted_unfiltered = parse_Transmitted(df_EG_merge)
    # transmitted_filtered_count = len(df_transmitted_filtered["SampleID"].unique())
    # logging.warning("transmitted filtered count is: " + str(transmitted_filtered_count))

    # transmitted_unfiltered_count = len(df_transmitted_unfiltered["SampleID"].unique())
    # logging.warning("transmitted unfiltered count is: " + str(transmitted_unfiltered_count))

    
    ###########################
    #UNCOMMENT
    #Merge transmitted and DNV analysis
    # df_analysis = pd.merge(df_transmitted_analysis, df_DNV_analysis, on=["Symbol", "gene_id"], how='left')
    # df_analysis.to_csv(output_path + "new_UCSF_total_filtered_analysis.txt", sep = "\t", header = True, index = None)

   
    ####################
    #Combine mutability into one dataframe
    # df_total_mutability = pd.merge(df_roulette, df_samocha, on = ["gene_id"], how = "left")
    ####################

    # df_roulette_samocha.to_csv(output_path + "gene_annotated_total_mutability.txt", sep = "\t", header = True, index = None)

    #Merge mutability to variant analysis dataframe
    # df_final_analysis = pd.merge(df_analysis, df_total_mutability, on = ["gene_id"], how = "left")
    # df_final_analysis.to_csv(output_path + "new_UCSF_total_mutability_analysis.txt", sep = "\t", header = True, index = None)

    # ###################
    # df_DNV_analysis = pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/new_UCSF_DNV_filtered_analysis.txt", sep = "\t")
    # df_final_analysis = pd.merge(df_analysis, df_roulette_samocha, on = ["Symbol"], how = "left")
    # df_final_analysis.to_csv(output_path + "new_UCSF_mutability_analysis.txt", sep = "\t", header = True, index = None)
    # ###################

    # # Settings
    # element2file = {"enhancer": "all_peaks_element_200bp_e.bed", "promoter": "all_peaks_element_200bp_p.bed"}
    # sex2label = {"all": {"male", "female", np.nan}, "male": {"male"}, "female": {"female"}}
    # label2sample = {"Cerebral Cortex": {"CHTN22", "HN21"}, "Cerebellum": {"CHTN13"}, "Neuronal Progenitor Cells": {"LUHMES.Li"}, "Astrocytes": {"SC3"}, "Excitatory Cortical Neurons": {"SC2"}, "Combined": {"CHTN22", "HN21", "CHTN13", "LUHMES.Li", "SC3", "SC2"}, "K562": {"K562"}, "HEK293": {"HEK293"}}
    # dirname = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429"
    # # Promoters
    # # promoter2file = {}
    # # for version in ["v43", "v21"]:
    # #     for transcript in ["coding", "all"]:
    # #         output = dirname + "/Promoters/promoters_1kb_tss_centered." + transcript + "." + version + ".bed"
    # #         promoter2file["promoter." + version + "." + transcript] = all_function_py3.get_GENCODE_promoters(output, transcripts = transcript, version = version)

    # # Offspring variants
    # var2file = {}
    # for coding_flag in ["damaging"]:#, "strict", "loose", "None"]:
    #     var2file["DNV_ExcludeCoding_" + str(coding_flag)] = dirname + "/UCSF_DNV_ExcludeCoding_" + coding_flag + ".txt"
    #     var2file["Transmitted_ExcludeCoding_" + coding_flag] = dirname + "/Variants/UCSF_Transmitted_ExcludeCoding_" + coding_flag + ".txt"

    # f_fig = open("/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/Enrichment_analysis_final.txt", "w")
    # f_fig.write("\t".join(["Var_type", "gnomADg_AF", "Sample", "Element", "Sex", "Case_var_in_enhancers", "Case_var", "Control_var_in_enhancers", "Control_var", "Odds_ratio", "P_value", "Case_individuals_in_enhancers", "Control_individuals_in_enhancers"]) + "\n")

    # promoter = "promoter.v43.coding"
    # CADD = 12
    # alternative = "greater"
    # var2mapping = {}
    # var2scale = {}

    # for var in var2file:
    #     if var[0:len("DNV")] == "DNV":
    #         gnomAD_list = [0] #, 1e-5, 1e-4, 1e-3, 1e-2]
    #         gnomAD_version = "r2.1.1"
    #     else:
    #         gnomAD_list = [0, np.inf]
    #         gnomAD_version = "r3.1"

    #     df1 = pd.read_csv(var2file[var], sep = "\t")
    #     case_var = set(df1[df1["Pheno"] == "case"]["#Uploaded_variation"])
    #     control_var = set(df1[df1["Pheno"] == "control"]["#Uploaded_variation"])
    #     case_only = case_var - control_var
    #     control_only = control_var - case_var
    #     for sex in sex2label:
    #         df1_sex = df1[df1["sex"].isin(sex2label[sex])]
    #         all_var_count = df1_sex[df1_sex["Pheno"] == "control"][["#Uploaded_variation", "SampleID"]].drop_duplicates().shape[0]
    #         df1_sex_damaging = Yingying_damaging_filter(df1_sex, CADD, phyloP = 2)
    #         for gnomAD in gnomAD_list:
    #             df1_sex_damaging_rare = df1_sex_damaging[df1_sex_damaging["gnomADg_AF"].astype(float) <= gnomAD]
    #             filtered_var_count = df1_sex_damaging_rare[df1_sex_damaging_rare["Pheno"] == "control"][["#Uploaded_variation", "SampleID"]].drop_duplicates().shape[0]
    #             var2scale[(var, gnomAD_version + ":" + str(gnomAD), sex)] = filtered_var_count / all_var_count
    #             b = df1_sex_damaging_rare[df1_sex_damaging_rare["#Uploaded_variation"].isin(case_only)][["chr", "pos", "ref", "alt", "SampleID"]].drop_duplicates().shape[0]
    #             d = df1_sex_damaging_rare[df1_sex_damaging_rare["#Uploaded_variation"].isin(control_only)][["chr", "pos", "ref", "alt", "SampleID"]].drop_duplicates().shape[0]
    #             for element in ["enhancer"]:
    #                 if var[0:len("DNV")] == "DNV":
    #                     file_name = var2file[var].split("/")[-1].replace(".txt", "_" + element + ".txt")
    #                 else:
    #                     file_name = "UCSF_Parent_" + element + ".txt"
    #                 for label in label2sample:
    #                     if label != "Combined" and label2sample[label] - label2sample["Combined"] == set():
    #                         sample_elements = []
    #                         for sample in label2sample[label]:
    #                             sample_elements.append(pd.read_csv(dirname + "/" + sample + "/" + promoter + "/" + element2file[element], sep = "\t", dtype = str))
    #                         sample_elements = pd.concat(sample_elements)
    #                         sample_elements = sample_elements.apply(lambda x: x["chr"] + ":" + x["start"] + "-" + x["end"], axis = 1).unique().tolist()
    #                         combined_var = pd.read_csv(dirname + "/" + "_".join(sorted(label2sample["Combined"])) + "/" + promoter + "/" + file_name, sep = "\t", dtype = str)
    #                         element_var = combined_var[combined_var["Overlapped_elements_" + element].isin(sample_elements)]
    #                     else:
    #                         element_var = pd.read_csv(dirname + "/" + "_".join(sorted(label2sample[label])) + "/" + promoter + "/" + file_name, sep = "\t", dtype = str)
    #                     if "#Uploaded_variation" not in element_var.columns:
    #                         element_var["#Uploaded_variation"] = element_var.apply(lambda x: "_".join([x["chr"][3:], x["pos"], x["ref"] + "/" + x["alt"]]), axis = 1)
    #                     element_var = element_var[["#Uploaded_variation", "Overlapped_elements_" + element]].drop_duplicates()
    #                     df1_sex_damaging_rare_element = pd.merge(df1_sex_damaging_rare, element_var, on = "#Uploaded_variation")
    #                     var2mapping[(var, gnomAD_version + ":" + str(gnomAD), sex, label, element)] = df1_sex_damaging_rare_element
    #                     a = df1_sex_damaging_rare_element[df1_sex_damaging_rare_element["#Uploaded_variation"].isin(case_only)][["chr", "pos", "ref", "alt", "SampleID"]].drop_duplicates().shape[0]
    #                     c = df1_sex_damaging_rare_element[df1_sex_damaging_rare_element["#Uploaded_variation"].isin(control_only)][["chr", "pos", "ref", "alt", "SampleID"]].drop_duplicates().shape[0]
    #                     a_individuals = df1_sex_damaging_rare_element[df1_sex_damaging_rare_element["#Uploaded_variation"].isin(case_only)][["SampleID"]].drop_duplicates().shape[0]
    #                     c_individuals = df1_sex_damaging_rare_element[df1_sex_damaging_rare_element["#Uploaded_variation"].isin(control_only)][["SampleID"]].drop_duplicates().shape[0]

    #                     OR, pval = scipy.stats.fisher_exact([[a, b - a], [c, d - c]], alternative = alternative)
    #                     f_fig.write("\t".join([var, gnomAD_version + ":" + str(gnomAD), label, element, sex, str(a), str(b), str(c), str(d), str(OR), str(pval) + str(a_individuals) + str(c_individuals)]) + "\n")
    # f_fig.close()
    # # 	with open(dirname + "/Scaling_factors.json", "w") as f:
    # # 		json.dump(var2scale, f)

    # # Summary
    # # promoter = "promoter.v43.coding"
    # output_path = dirname + "/" + "_".join(sorted(label2sample["Combined"])) + "/promoter.v43.coding"
    # df_conversion = pd.read_csv(output_path + "/" + element2file["enhancer"].replace(".bed", "_hg19.bed"), sep = "\t", dtype = str)
    # df1 = pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/CHTN13_CHTN22_HN21_LUHMES.Li_SC2_SC3/promoter.v43.coding/EG_connections/all_peaks_element_200bp_e_hg19_NoHeader__promoters_1kb_tss_centered.coding.v43.hg19.NoHeader__merged_interaction_profile_CP.txt", sep = "\t")
    # df2 = pd.read_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/CHTN13_CHTN22_HN21_LUHMES.Li_SC2_SC3/promoter.v43.coding/EG_connections/all_peaks_element_200bp_e_hg19_NoHeader__promoters_1kb_tss_centered.coding.v43.hg19.NoHeader__merged_interaction_profile_GZ.txt", sep = "\t")
    # # df1 = pd.read_csv(output_path + "/EG_connections/" + element2file["enhancer"].replace(".bed", "_hg19_NoHeader") + "__" + promoter2file[promoter].split("/")[-1].replace(".bed", ".hg19.NoHeader") + "__merged_interaction_profile_CP.txt", sep = "\t")
    # # df2 = pd.read_csv(output_path + "/EG_connections/" + element2file["enhancer"].replace(".bed", "_hg19_NoHeader") + "__" + promoter2file[promoter].split("/")[-1].replace(".bed", ".hg19.NoHeader") + "__merged_interaction_profile_GZ.txt", sep = "\t")
    # df = pd.concat([df1, df2])
 
    # df = df[df["chr"].isin(["chr" + str(i) for i in range(1, 23)])]
    # df_final = defaultdict(list)
    # for gene,df_one in df.groupby("7"):
    #     df_final["gene_id"].append(gene)
    #     elements = []
    #     for col in ["Distal_enhancers", "Proximal_enhancers"]:
    #         for row in df_one[col]:
    #             if type(row) == str:
    #                 elements.extend([":".join(x.split(":")[0:2]) for x in row.split(",")])
    #     elements = sorted(set(elements))
    #     if len(elements) > 0:
    #         df_tmp = pd.DataFrame({"elements": elements})
    #         df_tmp["chr"] = df_tmp["elements"].apply(lambda x: x.split(":")[0])
    #         df_tmp["start"] = df_tmp["elements"].apply(lambda x: x.split(":")[1].split("-")[0])
    #         df_tmp["end"] = df_tmp["elements"].apply(lambda x: x.split(":")[1].split("-")[1])
    #         df_tmp = pd.merge(df_tmp, df_conversion, on = ["chr", "start", "end"])
    #         df_tmp["ID"] = df_tmp.apply(lambda x: x["chr"] + ":" + x["start_hg38"] + "-" + x["end_hg38"], axis = 1)
    #         df_tmp[["chr", "start_hg38", "end_hg38"]].to_csv(gene + ".bed", sep = "\t", header = None, index = None)
    #         os.system("/programs/bin/bedtools/bin/bedtools sort -i " + gene + ".bed | /programs/bin/bedtools/bin/bedtools merge -i - > " + gene + "_merged.bed")
    #         df_merged = pd.read_csv(gene + "_merged.bed", sep = "\t", header = None, dtype = str)
    #         df_merged["region"] = df_merged.apply(lambda x: x[1] + "-" + x[2], axis = 1)
    #         df_final["enhancer_regions"].append(df_merged[0].tolist()[0] + ":" + ",".join(df_merged["region"].tolist()))
    #         os.remove(gene + ".bed")
    #         os.remove(gene + "_merged.bed")
    #         for var in var2file:
    #             if var[0:len("DNV")] == "DNV":
    #                 gnomAD_flag = "r2.1.1:0"
    #             else:
    #                 gnomAD_flag = "r3.1:inf"
    #             if var.split("_")[-1] == "damaging":
    #                 for sex in sex2label:
    #                     df_var_tmp = var2mapping[(var, gnomAD_flag, sex, "Combined", "enhancer")][var2mapping[(var, gnomAD_flag, sex, "Combined", "enhancer")]["Overlapped_elements_enhancer"].isin(df_tmp["ID"].tolist())]
    #                     for pheno in ["case", "control"]:
    #                         df_gene_var = df_var_tmp[df_var_tmp["Pheno"] == pheno][["#Uploaded_variation", "SampleID"]].drop_duplicates()
    #                         # df_gene_var["Linked_gene"] = gene
    #                         # for sample,df_sample in df_gene_var.groupby("SampleID"):
    #                         # 	if df_sample.shape[0] > 1:
    #                         # 		df_sample.to_csv(dirname + "/check.txt", mode = "a", header = False, index = None)
    #                         df_final["_".join([var, gnomAD_flag, sex, pheno])].append(len(df_gene_var["SampleID"].unique()))
    #     else:
    #         df_final["enhancer_regions"].append(np.nan)
    #         for var in var2file:
    #             if var[0:len("DNV")] == "DNV":
    #                 gnomAD_flag = "r2.1.1:0"
    #             else:
    #                 gnomAD_flag = "r3.1:inf"
    #             if var.split("_")[-1] == "damaging":
    #                 for sex in sex2label:
    #                     for pheno in ["case", "control"]:
    #                         df_final["_".join([var, gnomAD_flag, sex, pheno])].append(np.nan)
    # df_final = pd.DataFrame(df_final)
    # df_final.to_csv("/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/Summary_final.txt", sep = "\t", header = True, index = None)