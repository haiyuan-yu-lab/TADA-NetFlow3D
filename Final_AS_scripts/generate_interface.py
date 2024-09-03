import pandas as pd
import numpy as np
import sys
sys.path.append("/home/yz2296/bin")
# sys.path.append（./Run_NetFlow/all_function_py3.py)
import all_function_py3
import os
from collections import defaultdict
import copy

#place this file outside of Run_NetFlow folder
dirname = os.path.dirname(os.path.abspath(__file__))
# S6_file = dirname + "/Run_NetFlow/pDenovoWest.xlsx"
# netflow_path = "/local/storage/yz2296/cancer_hotspot_new/website/website/pipeline/"
netflow_path = "/lustre2/home/js2788/Desktop/Lab/pipeline/"

#Function to read pDenovoWEST data into dataframe S6
#Input filename for pDenovoWEST data
def read_S6_DeNovo(S6_file):

	S6 = S6_file
	S6 = S6[S6["pDenovoWEST"] != "."]
	S6["pDenovoWEST"] = S6["pDenovoWEST"].astype(float)
	S6["plog_DenovoWEST"] = S6["pDenovoWEST"].apply(lambda x: -np.log10(x))
	S6.to_csv(dirname + "/Run_NetFlow/DeNovo_Input/pDNW_log.txt", sep = "\t", header = True, index = None)
	return S6

	# Score = pd.read_excel(dirname + "/Run_NetFlow/41588_2022_1104_MOESM3_ESM.xlsx", sheet_name = "Supplementary Table 11", dtype = str)
	# Score = S6_file
	# Score = Score[Score["FDR_TADA_NDD"] != "."]
	# Score["FDR_TADA_NDD"] = Score["FDR_TADA_NDD"].astype(float)
	# for score in Score["FDR_TADA_NDD"]:
	# 	if score == 0:
	# 		score = 1e-30
	# Score["FDR_heat_score"] = Score["FDR_TADA_NDD"].apply(lambda x: -np.log10(x))
	# Score["p_TADA_NDD"] = Score["p_TADA_NDD"].astype(float)
	# for score in Score["p_TADA_NDD"]:
	# 	if score == 0:
	# 		score = 1e-30
	# Score["p_heat_score"] = Score["p_TADA_NDD"].apply(lambda x: -np.log10(x))
	# Score = Score.dropna(subset = ["FDR_TADA_NDD"])
	# return Score
	# ###################################
	# S6 = pd.read_excel(S6_file, sheet_name = "Table S6", skiprows = 2)
	
	# ###################################
	# # plog = []
	# # for x in S6["pDenovoWEST"]:
	# #	 log = 0 -(np.log(float(x)))
	# #	 plog.append(log)
	# # S6["plog"] = plog
	# ###################################


def read_FDR_scores(FDR_file):
	Score = FDR_file
	Score = Score[Score["FDR_TADA_NDD"] != "."]
	Score["FDR_TADA_NDD"] = Score["FDR_TADA_NDD"].astype(float)
	Score["FDR_TADA_NDD"] = Score["FDR_TADA_NDD"].replace(0, 1e-30)
	Score["FDR_heat_score"] = Score["FDR_TADA_NDD"].apply(lambda x: -np.log10(x))

	Score["p_TADA_NDD"] = Score["p_TADA_NDD"].astype(float)
	Score["p_TADA_NDD"] = Score["p_TADA_NDD"].replace(0, 1e-30)
	Score["p_heat_score"] = Score["p_TADA_NDD"].apply(lambda x: -np.log10(x))

	Score = Score.dropna(subset = ["FDR_TADA_NDD"])
	Score.to_csv(dirname + "/Run_NetFlow/DeNovo_Input/FDR_log.txt", sep = "\t", header = True, index = None)
	return Score

#change from ensp to uniprots
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

#Create dictionary for uniprot ID and corresponding residues
def res_dict(revel):
	residue_list = revel["Protein_position"].tolist()
	uniprot_list = revel["UniProt"].tolist()
	uniprot2pos = defaultdict(set)
	for i in range(revel.shape[0]):
		for prot in str(uniprot_list[i]).split(","):
			###################################
			# uniprot2pos[prot] = uniprot2pos[prot].union({residue_list[i]})
			uniprot2pos[prot] = uniprot2pos[prot].union(set(all_function_py3.split_consecutive_pos(residue_list[i])))
			###################################
	return uniprot2pos


#Function to determine whether the protein and its amino acid residues are on the interfaces
def whether_on_interface(row, uniprot2pos):
	###################################
	# prot1_ires = set(row["interface_residues1"][1:-1].split(","))
	# prot2_ires = set(row["interface_residues2"][1:-1].split(","))
	if row["interface_residues1"] != "[]":
		prot1_ires = set([int(x) for x in row["interface_residues1"][1:-1].split(",")])
	else:
		prot1_ires = set()
	if row["interface_residues2"] != "[]":
		prot2_ires = set([int(x) for x in row["interface_residues2"][1:-1].split(",")])
	else:
		prot2_ires = set()
	###################################
	on_int = []
	if prot1_ires.intersection(uniprot2pos[row["uniprot1"]]) != set():
		return True
	elif prot2_ires.intersection(uniprot2pos[row["uniprot2"]]) != set():
		return True
	else:
		return False

#Function to produce interface table
#Include protein uniprot ID
def Interface_Table(Interface, uniprot2pos, output):
	Interface = Interface[Interface.apply(lambda x: whether_on_interface(x, uniprot2pos), axis = 1)]
	Interface[["uniprot1", "uniprot2"]].rename(columns = {"uniprot1": "P1", "uniprot2": "P2"}).to_csv(output, sep = "\t", header = True, index = None)
	###################################
	# output_path = output_path + "/Run_NetFlow/DeNovo_Input/Interface_DN.txt"
	# P1 = list(Interface["uniprot1"])
	# P2 = list(Interface["uniprot2"])
	# Interaction_Table = pd.DataFrame({"P1":P1, "P2":P2})
	# Interaction_Table.to_csv(output_path, sep = "\t", header = True, index = None)
	###################################
	return Interface


#Function to create datatables for DeNovo based on pDenovoWEST p-values
#Input S6 = S6 data; p_value = p-value cutoff; score = "p" or "FDR"
def p_filter_table(S11, p_value = None, score = "p", output_path = dirname):
	###################################
	# if p_value == None and score != "p":
	###################################
	# if score == "heat" and p_value == np.inf:
	# 	return False
	# if score == "heat":
	# 	pDW = copy.deepcopy(S6)
	# 	Score = []
	# 	for p in pDW["pDenovoWEST"]:
	# 		if p < p_value:
	# 			Score.append(1)
	# 		else:
	# 			Score.append(0)
	# 	pDW["plog"] = Score  
	pscore = copy.deepcopy(S11)
	column = str(score)+ "_heat_score"
	pscore = S11[S11["FDR_TADA_NDD"] < p_value]

	# pscore[column] = pscore[str(score) + "_TADA_NDD"]
	# pscore.to_csv(dirname + "/Run_NetFlow/" + column + "_.txt", sep = "\t")

	# else:
		# pscore = S6[S6["plog_p"] < p_value]
		# value = "plog_p"
		###################################
		# if p_value == None:
		#	 pDW = S6
		# else:
		#	 pDW = S6[S6["pDenovoWEST"] < float(p_value)]
		###################################
	
	output_path = output_path + "/Run_NetFlow/DeNovo_Input/DN_" + column + "_p_" + str(p_value) + ".txt"
	pscore[["gene_id",column]].rename(columns = {"gene_id":"Gene", column :"Score"}).to_csv(output_path, sep = "\t", header = True, index = None)
	return output_path

def denovo_filter_table(S6, p_value, output_path = dirname):
	pscore = copy.deepcopy(S6)
	pscore["pDenovoWEST"] = pscore["pDenovoWEST"].astype(float)
	pscore = pscore[pscore["pDenovoWEST"] < p_value]

	output_path = output_path + "/Run_NetFlow/DeNovo_Input/DN_pDenovoWEST_" + str(p_value) + ".txt"
	pscore[["GeneID","plog_DenovoWEST"]].rename(columns = {"GeneID":"Gene", "plog_DenovoWEST" :"Score"}).to_csv(output_path, sep = "\t", header = True, index = None)
	return output_path

#Function to call NetFlow based on De Novo data and interface
#Input dataframe S6 which includes pDenovoWEST data and p values for ur, ar, and denovo
#Input dirname for path
def run_DN(S11, dirname):
	# Parse DeNovo Variants
	revel = pd.read_csv(dirname + "/Run_NetFlow/new_final_combined_RunNF.txt", sep = "\t", dtype = str)
	REVEL_uniprot = all_function_py3.gene_ID_conversion(revel["ENSP"], convert_from = "ensembl.protein", convert_to = "uniprot", species = "human")
	revel["UniProt"] = get_uniprot(revel["ENSP"].tolist(), REVEL_uniprot)
	###################################
	print("# of variants: " + str(revel.shape[0]))
	print("# of variants without UniProt annotation: " + str(revel["UniProt"].isna().sum()))
	# uniprot2pos = res_dict(revel)
	uniprot2pos = res_dict(revel.dropna(subset = ["UniProt"]))
	###################################
	
	#Import Pioneer data
	Interface = pd.read_csv(dirname + "/Run_NetFlow/PIONEER_human_veryhigh.txt", sep = "\t")
	###################################
	# Interface["On_Interface"] = Interface.apply(lambda x: whether_on_interface(x, uniprot2pos), axis = 1)
	# Interface_table = Interface_Table(Interface, uniprot2pos)
	###################################
	Interface_table = Interface_Table(Interface, uniprot2pos, dirname + "/Run_NetFlow/DeNovo_Input/new_Interface_DN.txt")
	
	# Run NetFlow_ASD
	# for p_value in [0.05, 0.001]:
	# 	for score in ["p", "FDR"]:
	# 		path = p_filter_table(S11, p_value, score)
	# 		if path != False:
	# 			edge_weight = 5
	# 			command = ("python " + netflow_path + "NetFlow_ASD.py -m " + path + " -e " + dirname + \
	# 				"/Run_NetFlow/DeNovo_Input/Interface_DN.txt -R low -I " + "_".join(["p", str("0.05"), str(edge_weight)]) + \
	# 				" -X " + dirname + "/Run_NetFlow/Brain_expressed_genes.txt -w " + str(edge_weight) + " -o " + dirname + \
	# 				"/Run_NetFlow/DeNovo_Output/")
	# 			os.system(command)
	# command = ("python " + netfl
	
def run_DN_WEST(S6, dirname):
	# Parse DeNovo Variants
	revel = pd.read_csv(dirname + "/Run_NetFlow/new_final_combined_RunNF.txt", sep = "\t", dtype = str)
	REVEL_uniprot = all_function_py3.gene_ID_conversion(revel["ENSP"], convert_from = "ensembl.protein", convert_to = "uniprot", species = "human")
	revel["UniProt"] = get_uniprot(revel["ENSP"].tolist(), REVEL_uniprot)
	###################################
	print("# of variants: " + str(revel.shape[0]))
	print("# of variants without UniProt annotation: " + str(revel["UniProt"].isna().sum()))
	# uniprot2pos = res_dict(revel)
	uniprot2pos = res_dict(revel.dropna(subset = ["UniProt"]))
	###################################
	
	#Import Pioneer data
	Interface = pd.read_csv(dirname + "/Run_NetFlow/PIONEER_human_veryhigh.txt", sep = "\t")
	###################################
	# Interface["On_Interface"] = Interface.apply(lambda x: whether_on_interface(x, uniprot2pos), axis = 1)
	# Interface_table = Interface_Table(Interface, uniprot2pos)
	###################################
	Interface_table = Interface_Table(Interface, uniprot2pos, dirname + "/Run_NetFlow/DeNovo_Input/new_Interface_DN.txt")

	# for p_value in [0.01, 0.001]:
	# 	score = "plog_DenovoWEST"
	# 	path = denovo_filter_table(S6, p_value)
	# 	if path != False:
	# 		edge_weight = 5
	# 		command = ("python " + netflow_path + "NetFlow_ASD.py -m " + path + " -e " + dirname + \
	# 			"/Run_NetFlow/DeNovo_Input/Interface_DN.txt -R low -I " + "_".join([score, str(p_value), str(edge_weight)]) + \
	# 			" -X " + dirname + "/Run_NetFlow/Brain_expressed_genes.txt -w " + str(edge_weight) + " -o " + dirname + \
	# 			"/Run_NetFlow/DeNovo_Output/")
	# 		os.system(command)


# Run using pDenovoWEST
# Score = pd.read_excel(dirname + "/Run_NetFlow/pDenovoWest.xlsx", sheet_name = "Table S6", skiprows = 2)
# S6 = read_S6_DeNovo(Score)
# run_DN_WEST(S6, dirname)

# Run using FDR_TADA_NDD
Score = pd.read_excel(dirname + "/Run_NetFlow/41588_2022_1104_MOESM3_ESM.xlsx", sheet_name = "Supplementary Table 11", dtype = str)
S11 = read_FDR_scores(Score)
run_DN(S11, dirname)



#Function to parse p-value data for all rare and ultra rare transmitted mutations
#Create data tables including gene symbol and heat scores for each gene
#Input S6 = data table S6; ar_ur = all rare or ultra rare; confidence = pExt value; cutoff = p-value cutoff
###################################
# def TDT_stats(S6, ar_ur, confidence, cutoff, output_path = dirname + "/Run_NetFlow/TDT_Input/"):
# 	'''
# 	ar_ur: (1) AllLoF; (2) URLoF;
# 	confidence: (1) 0.1; (2) 0.9; (3) None;
# 	cutoff: (1) None; (2) 0; (3) 1;
# 	'''
# 	if confidence != None:
# 		name = ar_ur + ",pExt>="+ str(confidence) + "TDTStat"
# 	else:
# 		name = ar_ur + "TDTStat"
# 	file = ar_ur + "_pExt"+ str(confidence) + "_TDT" + str(cutoff) + ".txt"
# 	df_pExt = S6[["GeneID", name]].rename(columns = {"GeneID": "Gene", name: "Score"})
# 	if cutoff == 1:
# 		df_pExt = df_pExt[df_pExt["Score"] >= cutoff]
# 	elif cutoff == 0:
# 		df_pExt = df_pExt[df_pExt["Score"] > cutoff]

# 	output_path = output_path + file
# 	df_pExt.to_csv(output_path, sep = "\t", header = True, index = None)
# 	return output_path
###################################


#Function to call NetFlow based on All rare and Ultra rare transmitted data with no Interface
#Input dataframe S6 which includes pDenovoWEST data and p values for ur, ar, and denovo
#Input dirname for path
# def run_TDT(S6, dirname):
# 	No_Interface = pd.DataFrame({"P1": [], "P2": []})
# 	No_Interface.to_csv(dirname + "/Run_NetFlow/No_Interface.txt", sep = "\t", header = True, index = None)
	
# 	for cutoff in [None, 1, 0]:
# 		for confidence in [0.1, 0.9, None]:
# 			for af in ["AllLoF", "URLoF"]:
# 				path = TDT_stats(S6, af, confidence, cutoff)
# 				command = "python " + netflow_path + "NetFlow_ASD.py -m " + path + " -e " + dirname + \
# 					"/Run_NetFlow/No_Interface.txt -R low -I " + "_".join([af, "pExt" + str(confidence), "TDT" + str(cutoff)]) + \
# 					" -X " + dirname + "/Run_NetFlow/Brain_expressed_genes.txt -o " + dirname + "/Run_NetFlow/TDT_Output/"
# 				os.system(command)
# 	return
	

# run_TDT(S6, dirname)



# # test = revel[revel["Protein_position"].apply(lambda x: str(x).find("-")!= -1)]
# # #parse protein_position count both ex: 559-560 or 186-188"



# # #Checking to see if there are more than two ENSP in a row
# # for values in revel["ENSP"]:
# #	 if len(values)>18:
# #		 print (values)


# # df = pd.read_csv("/Users/amandasun/Desktop/Lab/pipeline/test_subnetworks.txt", sep = "\t")



# def subnetwork(df):
# 	df = df[df["Subnetwork_size"] >5]
# 	gene_list = []
# 	for genes in df["Subnetwork_genes"]:
# 		for g in genes.split(","):
# 			gene_list.append(g)
# 	genes = all_function_py3.gene_ID_conversion(gene_list, convert_from = "uniprot", convert_to = "symbol", species = "human")
# 	glist = []
# 	for key in genes:
# 		for values in genes[key]:
# 			glist.append(values)
# 	return glist


# def ASD_count(gene_list, ASD):
# 	i = 0
# 	for gene in gene_list:
# 		if gene in ASD:
# 			i +=1
# 	return i

# def enrichment_num(dirname, df):
# 	#Read known ASD genes data into dataframe ASD
# 	ASD=pd.read_csv(dirname + "/Run_NetFlow/known_ASD_genes.txt", sep = "\t")
# 	#Keep gene symbol for each gene in column gene
# 	ASD["gene"] = ASD["gene"].apply(lambda x: str(x[2:]))
# 	#Turn gene column into list
# 	ASD = ASD["gene"].tolist()
# 	gene_list = subnetwork(df)
# 	result = ASD_count(gene_list, ASD)
# 	return result


# def input_gen(S6_file, dirname):
# 	S6 = read_S6(S6_file)
# 	# print(S6)
# 	run_TDT(S6, dirname)
# 	run_DN(S6, dirname)
# 	return True

# x = input_gen(S6_file, dirname)