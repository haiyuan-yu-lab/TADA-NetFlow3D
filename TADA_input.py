import pandas as pd
import numpy as np
import argparse
import sys
import importlib
sys.path.append("/home/yz2296/bin")
import all_function_py3
import os
from os import path
from glob import glob
import scipy.stats
import logging
from multiprocessing import Pool
from collections import defaultdict
import pysam
import gzip
import logging
import json
from functools import partial
import subprocess
from pybedtools import BedTool
dirname = os.path.dirname(os.path.abspath(__file__))

def damaging_filter(df_NC, CADD, phyloP):
	NC_SNV = df_NC[df_NC["TYPE"] == "SNV"]
	NC_SNV = NC_SNV[NC_SNV["CADD_PHRED"].astype(float) >= float(CADD)]
	NC_Indel = df_NC[df_NC["TYPE"] == "Indel"]
	NC_Indel = NC_Indel[NC_Indel["phyloP100way"].astype(float) >= float(phyloP)]
	return pd.concat([NC_SNV, NC_Indel])

def get_var_count(row, tag, df_group):
	try:
		df_tmp = df_group.get_group(row["chr"])
		return df_tmp[(df_tmp["pos"] > row["start"]) & (df_tmp["pos"] <= row["end"])][["#Uploaded_variation", "SampleID"]].drop_duplicates()
	except:
		logging.warning("For " + tag + ", no variant is found on " + row["chr"] + "!")
		return pd.DataFrame({"#Uploaded_variation": [], "SampleID": []})

def get_gene2enhancers(row):
	for x in row["Linked_genes"].split(","):
		gene2enhancers[x].append(row["enhancer"])
	return True


# Settings
element2file = {"enhancer": "all_peaks_element_200bp_e.bed"}#, "promoter": "all_peaks_element_200bp_p.bed"}
sex2label = {"all": {"male", "female", np.nan}, "male": {"male"}, "female": {"female"}}
label2sample = {"Cerebral Cortex": {"CHTN22", "HN21"}, "Cerebellum": {"CHTN13"}, "Neuronal Progenitor Cells": {"LUHMES.Li"}, "Astrocytes": {"SC3"}, "Excitatory Cortical Neurons": {"SC2"}, "Combined": {"CHTN22", "HN21", "CHTN13", "LUHMES.Li", "SC3", "SC2", "Neuron4R", "Microglia"}, "K562": {"K562"}, "HEK293": {"HEK293"}, "Neuron4R": {"Neuron4R"}, "Microglia": {"Microglia"}}
var2file = {x: dirname + "/UCSF_" + x + "_ExcludeCoding_damaging.txt" for x in ["DNV", "Transmitted", "Untransmitted"]}
af_cutoff = 0
cadd_cutoff = 12
phylop_cutoff = 2
chroms = ["chr" + str(i) for i in range(1, 23)]
file_path = dirname + "/../Week8_429/" + "_".join(sorted(label2sample["Combined"])) + "/promoter.v43.coding/"
df_shan = pd.read_csv("/local/storage/NetFlow_ASD/Mutability/20230919_samocha14_calcBaseLineRate_EG_mapping_allsample.bed", sep = "\t", dtype = str)
df_shan["enhancer"] = df_shan.apply(lambda x: x["chrom"] + ":" + x["start"] + "-" + x["end"], axis = 1)

# #------------------------------------------------------------
# # TADA coding input
# df = pd.read_excel("/local/storage/NetFlow_ASD/TADA_DAWN/fu_Table.xlsx", sheet_name = "Supplementary Table 5").dropna(subset = "gene_id")
# df_prior = pd.read_excel("/local/storage/NetFlow_ASD/TADA_DAWN/fu_Table.xlsx", sheet_name = "Supplementary Table 8").dropna(subset = ["gene_id"])
# df = df.rename(columns = {"gene_gencodeV33": "Hugo_Symbol"})[["Hugo_Symbol", "gene_id", "LOEUF", "pLI", "mut.ptv", "mut.misb", "mut.misa", "dn.ptv", "dn.ptv.sib", "dn.misb", "dn.misb.sib", "dn.misa", "dn.misa.sib", "in.t.ptv", "in.u.ptv", "in.t.misb", "in.u.misb", "in.t.misa",	"in.u.misa"]]
# df = pd.merge(df, df_prior.drop(columns = ["chromosome", "pLI",	"LOEUF", "gene", "gene_gencodeV33"]), on = "gene_id")
# cols = df.columns.tolist()[0:4]
# for x in ["prior", "mut", "dn", "in"]:
#     cols.extend([y for y in df.columns if y.split(".")[0] == x])
# df = df[cols]
# df.to_csv("/local/storage/NetFlow_ASD/TADA_DAWN/TADA_input_coding_FuRaw.txt", sep = "\t", header = True, index = None)
# for col in df.columns:
#     if col.split(".")[0] == "dn":
#         df[col] = df[col].apply(lambda x: round(float(x) * 4467 / 7898))
#     elif col.split(".")[0] == "in":
#         df[col] = df[col].apply(lambda x: round(float(x) * 4794 / 7898))
# df.to_csv("/local/storage/NetFlow_ASD/TADA_DAWN/TADA_input_coding_FuScaled.txt", sep = "\t", header = True, index = None)

# # # preliminary pass (needs further polishment)
# # df_our = pd.read_csv("/local/storage/NetFlow_ASD/TADA_DAWN/TADA_input_scaled_withsib.txt", sep = "\t", dtype = str)
# # df_our["gene_id"] = df_our["gene_id"].apply(lambda x: x.split(".")[0])
# # df_merge = pd.merge(df, df_our, on = "gene_id", how = "outer")
# # col = df_merge.columns.tolist()[0:7] + ["Noncoding_Mutability_samocha14"] + df_merge.columns.tolist()[7:19] + [x for x in df_merge.columns.tolist()[23:] if x.find("all") != -1]
# # df_merge[col].drop(columns = ["Noncoding_Transmitted_all_control", "Noncoding_Untransmitted_all_control"]).to_csv("/local/storage/NetFlow_ASD/TADA_DAWN/TADA_input_preliminary.txt", sep = "\t", header = True, index = None)


# # #------------------------------------------------------------
# # # TADA noncoding input
# # df = pd.read_csv(var2file["DNV"], sep = "\t", dtype = str)
# # total = df.shape[0]
# # df = df[df["gnomADg_AF"].astype(float) <= af_cutoff]
# # df = damaging_filter(df, cadd_cutoff, phylop_cutoff)
# # filtered = df.shape[0]
# # rare_damaging_frac = filtered / total
# rare_damaging_frac = 0.058980066340464904

# for element,file in element2file.items():
# 	df_bygene = {"gene_id": [], "enhancers": []}
# 	var2df = {}
# 	for var in var2file:
# 		if var == "DNV":
# 			df = pd.read_csv(file_path + "UCSF_" + var + "_ExcludeCoding_damaging_" + element + ".txt", sep = "\t", dtype = str)
# 		else:
# 			df = pd.read_csv(file_path + "UCSF_" + var + "_ExcludeCoding_damaging_" + element + "_6.txt", sep = "\t", dtype = str)
# 		df = df[df["gnomADg_AF"].astype(float) <= af_cutoff]
# 		df = damaging_filter(df, cadd_cutoff, phylop_cutoff)
# 		df["varID"] = df.apply(lambda x: "_".join([str(x[col]) for col in ["#Uploaded_variation", "SampleID","CADD_PHRED","phyloP100way","sex"]]), axis = 1)
# 		for key,df_one in df.groupby("Pheno"):
# 		    df_one[["chr", "start", "end", "varID"]].drop_duplicates().to_csv("/home/yz2296/Burden_test_for_Autism/data/" + "_".join([element, var, key]) + ".bed", sep = "\t", header = None, index = None)
# 		df["pos"] = df["pos"].astype(int)
# 		for sex,label in sex2label.items():
# 			for pheno in ["case", "control"]:
# 				tag = "_".join([var, pheno, sex])
# 				df_tmp = df[(df["Pheno"] == pheno) & df["sex"].isin(label)]
# 				var2df[tag] = df_tmp.groupby("chr")
# 				df_bygene[tag] = []
	
# 	df = pd.read_csv(file_path + file.replace(".bed", "_merged.bed"), sep = "\t")
# 	df = df[df["chr"].isin(chroms)]
# 	df[element] = df.apply(lambda x: x["chr"] + ":" + str(x["start"]) + "-" + str(x["end"]), axis = 1)
# 	for var in var2file:
# 		for sex in sex2label:
# 			for pheno in ["case", "control"]:
# 				tag = "_".join([var, pheno, sex])
# 				df[tag] = df.apply(lambda x: get_var_count(x, tag, var2df[tag]).shape[0], axis = 1)
# 	df.to_csv("/fs/cbsuhy01/storage/NetFlow_ASD/TADA_DAWN/" + element + "_mutcount.txt", sep = "\t", header = True, index = None)
# 	command = "python EG_linkage.py -i /fs/cbsuhy01/storage/NetFlow_ASD/TADA_DAWN/" + element + "_mutcount.txt -o /local/storage/NetFlow_ASD/TADA_DAWN/"
# 	os.system(command)

# 	#(needs further polishment)
# 	def get_samocha(x, df_shan, frac):
# 		if type(x) == str:
# 			df = df_shan[df_shan["enhancer"].isin(x.split(","))]
# 			return df["baselinerate_samocha14"].astype(float).max() * frac
# 		else:
# 			return np.nan
# 	df = pd.read_csv("/local/storage/NetFlow_ASD/TADA_DAWN/linked_enhancers.txt", sep = "\t", dtype = str).dropna(subset = ["Linked_genes"])
# 	df["mut.em"] = df["Overlapped_elements"].apply(lambda x: get_samocha(x, df_shan, rare_damaging_frac))
# 	df = df[["enhancer", "Linked_genes", "mut.em", "DNV_case_all", "DNV_control_all", "Transmitted_case_all", "Untransmitted_case_all"]].rename(columns = {"DNV_case_all": "dn.em", "DNV_control_all": "dn.em.sib", "Transmitted_case_all": "in.t.em", "Untransmitted_case_all": "in.u.em"})
# 	df.to_csv("/local/storage/NetFlow_ASD/TADA_DAWN/TADA_input_formal.txt", sep = "\t", header = True, index = None)

# 	df = pd.read_csv("/local/storage/NetFlow_ASD/TADA_DAWN/TADA_input_formal.txt", sep = "\t").dropna() #(needs further polishment)
# 	gene2enhancers = defaultdict(list)
# 	df["status"] = df.apply(get_gene2enhancers, axis = 1)
# 	for gene,enhancers in gene2enhancers.items():
# 	    enhancers = sorted(set(enhancers))
# 	    df_bygene["gene_id"].append(gene)
# 	    df_bygene["enhancers"].append(",".join(enhancers))
# 	    df_tmp = pd.DataFrame({"enhancers": enhancers})
# 	    df_tmp[["chr", "start", "end"]] = pd.DataFrame(df_tmp["enhancers"].apply(lambda x: [x.split(":")[0], int(x.split(":")[1].split("-")[0]), int(x.split(":")[1].split("-")[1])]).tolist(), index = df_tmp.index)
# 	    for var in var2file:
# 	        for sex in sex2label:
# 	            for pheno in ["case", "control"]:
# 	                tag = "_".join([var, pheno, sex])
# 	                df_bygene[tag].append(pd.concat(df_tmp.apply(lambda x: get_var_count(x, tag, var2df[tag]), axis = 1).tolist()).drop_duplicates().shape[0])
# 	df_bygene = pd.DataFrame(df_bygene)
# 	df_bygene.to_csv("/local/storage/NetFlow_ASD/TADA_DAWN/TADA_count.txt", sep = "\t", header = True, index = None)
# 	df_bygene = pd.read_csv("/local/storage/NetFlow_ASD/TADA_DAWN/TADA_count.txt", sep = "\t")[["gene_id",	"enhancers", "DNV_case_all", "DNV_control_all", "Transmitted_case_all", "Untransmitted_case_all"]].rename(columns = {"DNV_case_all": "dn.em", "DNV_control_all": "dn.em.sib", "Transmitted_case_all": "in.t.em", "Untransmitted_case_all": "in.u.em"})
# 	df2 = pd.read_csv("/local/storage/NetFlow_ASD/TADA_DAWN/FuRaw_in.coding_.txt.txt", sep = "\t", dtype = str)
# 	df_merge = pd.merge(df1, df2, on = "gene_id", how = "right")
# 	cols = ["gene_gencodeV33", "gene_id", "LOEUF", "enhancers"] + [x for x in df_merge.columns if x.split(".")[0].split("_")[0] in {"dn", "in", "FDR"}]
# 	# check: df_merge[col]????


