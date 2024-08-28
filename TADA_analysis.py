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
import itertools
from pybedtools import BedTool
dirname = os.path.dirname(os.path.abspath(__file__))

def extract_uniprot_ids(fasta_file):
	uniprot_ids = []
	with open(fasta_file, 'r') as file:
		for line in file:
			if line.startswith('>'):
				# Split the header line and extract the UniProt ID
				parts = line.split('|')
				uniprot_id = parts[1]
				uniprot_ids.append(uniprot_id)
	return uniprot_ids

def get_uniprot2gene(row):
	uniprot2gene[row[0]] = row[2]
	return True

def get_gene2enhancers(row):
    for x in row["Linked_genes"].split(","):
        gene2enhancers[x].append(row["enhancer"])
    return True

# #------------------------------------------------------------
# settings
flag_list = []
status_list = []
for var_source in ["dn", "in"]:  
  for var_type in ["ptv", "misb", "misa", "em"]:
      flag_list.append("_".join([var_source, var_type]))
      status_list.append(False)



# # #------------------------------------------------------------
# # run TADA
# g_dn_em = 5
# g_in_em = 1.2
# df = pd.read_csv("/local/storage/NetFlow_ASD/TADA_DAWN/prior.txt", sep = "\t").set_index("variable")
# df.loc["gamma.dn.em", "value"] = g_dn_em
# df.loc["gamma.in.em", "value"] = g_in_em
# df.reset_index().to_csv("/local/storage/NetFlow_ASD/TADA_DAWN/prior.txt", sep = "\t", header = True, index = None)
# for dataset in ["FuRaw", "FuScaled"]:
# 	half_command = 'Rscript Yingying_TADA.R --coding_file TADA_input_coding_' + dataset + '.txt \
# 	--noncoding_file TADA_input_formal_tmp.txt --prior_file prior.txt --info_file info.txt \
# 	--sample_file sample_' + dataset + '.txt --output_file TADA_output/' + dataset + '_'
# 	for combo in list(itertools.product([True, False], repeat = len(flag_list))):
# 		flag = ""
# 		for i,x in enumerate(combo):
# 			if x == True:
# 				flag = flag + " --" + flag_list[i]
# 		command = half_command + "_".join([str(y) for y in combo]) + ".txt" + flag
# 		os.system(command)


# # #------------------------------------------------------------
# TADA analysis
df_info = pd.read_csv(dirname + "/info.txt", sep = "\t")[["gene", "gene_gencodeV33", "gene_id", "chromosome", "pLI", "LOEUF"]]
threshold = df_info['LOEUF'].quantile(0.1)
df_info["LOEUF_lowest_decile"] = df_info["LOEUF"] <= threshold
type_list = []
gene_list = set()
for var_type in ["ptv", "misb", "misa", "em"]:
    type_list.append(var_type)
    index = [i for i,x in enumerate(flag_list) if x.split("_")[1] == var_type]
    for x in index:
        status_list[x] = True
    df = pd.read_csv(dirname + "/TADA_output/FuRaw_" + "_".join([str(x) for x in status_list]) + ".txt", sep = "\t") 
    gene_list = gene_list.union(set(df[df["fdr"] < 0.05]["gene_id"]))
    df_info = pd.merge(df_info, df[["gene_id", "fdr"]].rename(columns = {"fdr": "FDR_" + "+".join(type_list)}), on = "gene_id")    
df_info = df_info[df_info["gene_id"].isin(gene_list)]

df = pd.read_csv(dirname + "/TADA_input_formal_tmp.txt", sep = "\t")
gene2enhancers = defaultdict(list)
df["status"] = df.apply(get_gene2enhancers, axis = 1)
cols = ["dn.em", "dn.em.sib", "in.t.em", "in.u.em"]
df_output = {x: [] for x in ["gene_id", "enhancers"] + ["sum." + y for y in cols] + ["sep." + y for y in cols]}
for gene,enhancers in gene2enhancers.items():
    df_tmp = df[df["enhancer"].isin(enhancers)]
    df_output["gene_id"].append(gene)
    df_output["enhancers"].append(",".join(df_tmp["enhancer"].tolist()))
    for col in cols:
        df_output["sep." + col].append(",".join([str(x) for x in df_tmp[col]]))
        df_output["sum." + col].append(df_tmp[col].sum())
df_output = pd.DataFrame(df_output)
df_output = pd.merge(df_info, df_output, on = "gene_id", how = "left")
df = pd.read_csv(dirname + "/TADA_input_coding_FuRaw.txt", sep = "\t")[["gene_id", "dn.ptv", "dn.ptv.sib", "dn.misb", "dn.misb.sib", "dn.misa", "dn.misa.sib", "in.t.ptv", "in.u.ptv", "in.t.misb", "in.u.misb", "in.t.misa", "in.u.misa"]]
df_output = pd.merge(df_output, df, on = "gene_id", how = "left")[df_info.columns.tolist() + ["enhancers", "sep.dn.em", "sep.dn.em.sib", "sum.dn.em", "sum.dn.em.sib", "sep.in.t.em", "sep.in.u.em", "sum.in.t.em", "sum.in.u.em", "dn.ptv", "dn.ptv.sib", "dn.misb", "dn.misb.sib", "dn.misa", "dn.misa.sib", "in.t.ptv", "in.u.ptv", "in.t.misb", "in.u.misb", "in.t.misa", "in.u.misa"]]
df_output.to_excel(dirname + '/gene_discovery.xlsx', index = False)




# df_info = pd.read_csv("/local/storage/NetFlow_ASD/TADA_DAWN/info.txt", sep = "\t")
# known_asd_genes = set(pd.read_csv("/fs/cbsuhy02/storage/yz2296/Burden_test_for_Autism/Paper/SFARI_GeneList/Siwei_GeneList.txt", sep = "\t", header = None)[0])
# all_genes = set(df_info["gene_gencodeV33"])
# c2 = len(known_asd_genes.intersection(all_genes))
# n2 = len(all_genes)
# suffix2flag = {"in.coding": "--in_coding"}#, "": []}
# for tag in ["FuRaw"]:#, "FuScaled"]:
#	 for suffix in suffix2flag:
#		 for analysis in [".txt"]:#, "_sanity.txt"]:
#			 df1 = pd.read_csv("/local/storage/NetFlow_ASD/TADA_DAWN/TADA_output_" + "_".join(["NA", "NA", suffix, tag]) + analysis, sep = "\t")
#			 df_bernie = df1[["gene_id", "fdr"]].rename(columns = {"fdr": "FDR_g.dn.em.NA_g.in.em.NA"})
#			 df1 = df1[df1["fdr"] < 0.05]
#			 final_genes = set(df1["gene_id"])
#			 coding_genes = set(df1["Hugo_Symbol"])
#			 c1 = len(coding_genes.intersection(known_asd_genes))
#			 n1 = len(coding_genes)
#			 enr, log_enr, se_log_enr, pval = all_function_py3.calc_enr(c1, n1, c2, n2)
#			 fig_title = "dn.coding"
#			 if suffix == "in.coding":
#				 fig_title = fig_title + "+" + suffix
#			 df1["lbf_" + fig_title] = df1["bf"].apply(lambda x: np.log(x))
#			 df1["-l10q_" + fig_title] = df1["fdr"].apply(lambda x: -np.log10(x))
#			 df_output = {"category": [fig_title], "num_unique_genes": [0], "c1": [c1], "n1": [n1], "c2": [c2], "n2": [n2], "log_enr": [log_enr], "se_log_enr": [se_log_enr], "pval": [pval], "color": ["#b9b190"]}
#			 for g_dn_em in [15, 10, 5]:#, 100, 200, 500]:
#				 for g_in_em in [1.59, 1.2, "NA"]:#["NA", 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 50]:
#					 df2 = pd.read_csv("/local/storage/NetFlow_ASD/TADA_DAWN/TADA_output_" + "_".join([str(g_in_em), str(g_dn_em), suffix, tag]) + analysis, sep = "\t")
#					 df_bernie = pd.merge(df_bernie, df2[["gene_id", "fdr"]].rename(columns = {"fdr": "FDR_g.dn.em." + str(g_dn_em) + "_g.in.em." + str(g_in_em)}), on = "gene_id")
#					 df2 = df2[df2["fdr"] < 0.05]
#					 final_genes = final_genes.union(set(df2["gene_id"]))
#					 identified_genes = set(df2["Hugo_Symbol"])
#					 title = ["dn.coding"]
#					 if suffix == "in.coding":
#						 title.append(suffix)
#					 category = []
#					 if g_dn_em != "NA":
#						 title.append("dn.noncoding")
#						 category.append("g.dn.em=" + str(g_dn_em))
#					 if g_in_em != "NA":
#						 title.append("in.noncoding")
#						 category.append("g.in.em=" + str(g_in_em))
#					 df2["lbf_" + "+".join(title)] = df2["bf"].apply(lambda x: np.log(x))
#					 df2["-l10q_" + "+".join(title)] = df2["fdr"].apply(lambda x: -np.log10(x))
#					 df_merge = pd.merge(df1[["Hugo_Symbol", "gene_id", "LOEUF","lbf_" + fig_title, "-l10q_" + fig_title]], df2[["gene_id", "lbf_" + "+".join(title), "-l10q_" + "+".join(title)]], on = "gene_id", how = "outer")
#					 for metric in ["-l10q_"]:#["lbf_", "-l10q_"]:
#						 all_function_py3.scatter_plot(df_merge[metric + fig_title].tolist(), df_merge[metric + "+".join(title)].tolist(), output = None, \
#													   title = ",".join([tag, "g.dn.em=" + str(g_dn_em), "g.in.em=" + str(g_in_em)]), xlabel=metric + fig_title, ylabel=metric + "\n+".join(title), fontsize = 15, figsize = (3.5,3.5))
#					 # all_function_py3.venn2_diagram(coding_genes, identified_genes, fig_title.replace("+", "\n+"), "\n+".join(title), output = "/local/storage/NetFlow_ASD/Figures/venn2_" + "_".join([tag, suffix, str(g_in_em), str(g_dn_em)]) + ".png", title=tag, figsize=(2,2), show_text = True)
#					 all_function_py3.venn2_diagram(coding_genes, identified_genes, "coding", "coding+noncoding\n" + "\n".join(category), output = "/local/storage/NetFlow_ASD/Figures/venn2_" + "_".join([tag, suffix, str(g_in_em), str(g_dn_em)]) + ".png", title=tag, figsize=(2,2), show_text = True)
#					 unique_genes = identified_genes - coding_genes
#					 c1 = len(unique_genes.intersection(known_asd_genes))
#					 n1 = len(unique_genes)
#					 enr, log_enr, se_log_enr, pval = all_function_py3.calc_enr(c1, n1, c2, n2)
#					 df_output["category"].append(",".join(category))
#					 # df_output["category"].append("+".join(title))
#					 df_output["num_unique_genes"].append(len(unique_genes))
#					 df_output["color"].append("#75c8ae")
#					 for col in df_output:
#						 if col not in {"category", "num_unique_genes", "color"}:
#							 df_output[col].append(eval(col))
#			 df_output = pd.DataFrame(df_output)
#			 df_bernie = pd.merge(df_info[["gene_gencodeV33", "gene_id", "LOEUF"]], df_bernie, on = "gene_id")
#			 df_bernie = df_bernie[df_bernie["gene_id"].isin(final_genes)]
#			 df_bernie.to_csv(output_path + "_".join([tag, suffix, analysis]) + ".txt", sep = "\t", header = True, index = None)
#			 # df_bernie.to_excel(output_path + "_".join([tag, suffix, analysis]) + ".xlsx", index=False)
#			 # df_output.to_csv("/local/storage/NetFlow_ASD/Figures/enr_" + "_".join([tag, analysis]) + ".txt", sep = "\t", header = True, index = None)
#			 all_function_py3.bar_graph([(str(x),i) for i,x in enumerate(df_output["category"])], \
#										df_output["log_enr"].tolist(), df_output["color"].tolist(), \
#										output="/local/storage/NetFlow_ASD/Figures/enr_" + "_".join([tag, suffix]) + ".png", \
#										y_err=df_output["se_log_enr"].tolist(), figsize = (5,5), \
#										ylabel="log(Enrichment of SFARI genes)", \
#										title = ",".join([tag, fig_title]), xtick_rotation = 30, \
#										xtick_align_ha = "right", xtick_align_va = "top")



# # #------------------------------------------------------------
# # run NetFlow_ASD
# expr2flag = {"without.expr": "", "with.expr": "-X /local/storage/NetFlow_ASD/Week8_429/NetFlow_ASD/metadata/Brain_expressed_genes.txt"}
# for noncoding_status,(g_dn_em, g_in_em) in noncoding_status2gamma.items():
#   for coding_status,flag in coding_status2flag.items():
#	   for dataset in ["FuScaled"]:#, "FuRaw"]:
#		   df = pd.read_csv("/local/storage/NetFlow_ASD/TADA_DAWN/TADA_output/TADA_output_" + "_".join([coding_status, noncoding_status, dataset]) + ".txt", sep = "\t").rename(columns = {"gene_id": "Gene"})
#		   df["Score"] = df["bf"].apply(lambda x: np.log(x))
#		   # df["Score"] = df["fdr"].apply(lambda x: -np.log10(x))
#		   # df["Score"] = df["pval"].apply(lambda x: -np.log10(x))
#		   df[["Gene", "Score"]].to_csv("_".join([coding_status, noncoding_status, dataset]) + "_tmp.txt", sep = "\t", header = True, index = None)
#		   for edge_weight in [5, 1]:
#			   for expr,flag in expr2flag.items():
#				   tag = "_".join([coding_status, noncoding_status, dataset, "weight." + str(edge_weight), expr])
#				   command = "python /local/storage/NetFlow_ASD/Week8_429/NetFlow_ASD/NetFlow_ASD.py -I " + tag + \
#				   " " + flag + " -m " + "_".join([coding_status, noncoding_status, dataset]) + "_tmp.txt \
#				   -e /local/storage/NetFlow_ASD/Week8_429/NetFlow_ASD/input/Interface_DN.txt -w " + \
#				   str(edge_weight) + " -o /local/storage/NetFlow_ASD/TADA_DAWN/NetFlow3D_output/ -t 2"
#				   os.system(command)
#		   os.remove("_".join([coding_status, noncoding_status, dataset]) + "_tmp.txt")


# # #------------------------------------------------------------
# # thesis analysis

# df_annot = pd.read_csv("/local/storage/yz2296/cancer_hotspot_new/website/website/revision/NetFlow3D/metadata/HUMAN_9606_idmapping.dat.gz", sep = "\t", header = None)
# df_annot[0] = df_annot[0].apply(lambda x: x.split("-")[0])
# canonical_isoforms = extract_uniprot_ids("/local/storage/yz2296/cancer_hotspot_new/website/website/revision/NetFlow3D/metadata/UP000005640_9606.fasta")
# df_annot = df_annot[df_annot[0].isin(canonical_isoforms) & (df_annot[1] == "Gene_Name")]
# uniprot2gene = {}
# df_annot["status"] = df_annot.apply(get_uniprot2gene, axis = 1)
# file_path = dirname + "/NetFlow3D_output/"
# df_output = {"dataset": [], "expr_status": [], "coding_status": [], "noncoding_status": [], "network_status": [], "structural_status": [], "number_of_identified_genes": [], "identified_genes": []}
# df_full = {"dataset": [], "expr_status": [], "coding_status": [], "noncoding_status": [], "network_status": [], "structural_status": [], "gene": []}
# for x in os.listdir(file_path):
#	 if x[-len(".txt"):] == ".txt" and x.find("BE") == -1:
#		 df_output["network_status"].append(True)
#		 coding_status, noncoding_status, dataset, weight, expr_status, _ = x.split("_")
#		 df = pd.read_csv(file_path + x, sep = "\t")
#		 df = df[df["Subnetwork_size"] > 5]
#		 module_uniprots = set()
#		 for y in df["Subnetwork_genes"]:
#			 module_uniprots = module_uniprots.union(set(y.split(",")))
#		 module_genes = {uniprot2gene[y] for y in module_uniprots}
#		 df_output["number_of_identified_genes"].append(len(module_genes))
#		 df_output["identified_genes"].append(",".join(sorted(module_genes)))
#		 df_full["gene"].extend(sorted(module_genes))
#		 df_full["network_status"].extend([True] * len(module_genes))
#		 for col in ["coding_status", "noncoding_status", "expr_status"]:
#			 if eval(col).split(".")[0] == "with":
#				 df_output[col].append(True)
#				 df_full[col].extend([True] * len(module_genes))
#			 else:
#				 df_output[col].append(False)
#				 df_full[col].extend([False] * len(module_genes))
#		 df_output["dataset"].append(dataset)
#		 df_full["dataset"].extend([dataset] * len(module_genes))
#		 if weight.split(".")[1] == "1":
#			 df_output["structural_status"].append(False)
#			 df_full["structural_status"].extend([False] * len(module_genes))
#		 else:
#			 df_output["structural_status"].append(True)
#			 df_full["structural_status"].extend([True] * len(module_genes))
# file_path = dirname + "/TADA_output/"
# for x in os.listdir(file_path):
#	 if x[-len(".txt"):] == ".txt":
#		 df_output["network_status"].append(False)
#		 df_output["structural_status"].append(False)
#		 df_output["expr_status"].append(False)
#		 _, _, coding_status, noncoding_status, dataset = x.strip(".txt").split("_")
#		 df = pd.read_csv(file_path + x, sep = "\t")
#		 identified_genes = sorted(df[df["fdr"] < 0.05]["Hugo_Symbol"])
#		 df_output["number_of_identified_genes"].append(len(identified_genes))
#		 df_full["gene"].extend(identified_genes)
#		 df_full["network_status"].extend([False] * len(identified_genes))
#		 df_full["structural_status"].extend([False] * len(identified_genes))
#		 df_full["expr_status"].extend([False] * len(identified_genes))
#		 df_output["identified_genes"].append(",".join(identified_genes))
#		 for col in ["coding_status", "noncoding_status"]:
#			 if eval(col).split(".")[0] == "with":
#				 df_output[col].append(True)
#				 df_full[col].extend([True] * len(identified_genes))
#			 else:
#				 df_output[col].append(False)
#				 df_full[col].extend([False] * len(identified_genes))
#		 df_output["dataset"].append(dataset)
#		 df_full["dataset"].extend([dataset] * len(identified_genes))

# df_output = pd.DataFrame(df_output)
# df_full = pd.DataFrame(df_full)
# df_output.to_csv(dirname + "/summary.txt", sep = "\t", header = True, index = None)
# df_full.to_csv(dirname + "/details.txt", sep = "\t", header = True, index = None)