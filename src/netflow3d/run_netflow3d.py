#!/usr/bin/env python3
"""Provisional NetFlow3D command-line implementation."""

description = "Identify heat-diffusion subnetworks from gene-level scores."

import os
import argparse
from multiprocessing import Pool
import networkx as nx
import pandas as pd
import numpy as np
import network_functions as funcs
from os import path
import sys

def get_one_delta_multirun(args):
	return funcs.get_one_delta(*args)

file_path = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
	# user input
	parser = argparse.ArgumentParser(description = description)
	parser.add_argument('-I','--job_name', required = True, type = str, help = 'Please specify a name for your job. The output files will be stored in a folder with this name')
	parser.add_argument('-m','--gene_score', required = True, type = str, help = 'Gene score file with two columns: 1. Gene; 2. Score;')
	parser.add_argument('-e','--disrupted_interaction', type = str, default = None, help = 'Interaction file with two columns: 1. P1; 2. P2;')
	parser.add_argument('-X','--expressed_genes', type = str, default = None, help = '[OPTIONAL] A text file with the genes expressed in your context. Each row stands for a gene (non-UniProt IDs are accepted). No header line. By default all genes are considered as expressed.')
	parser.add_argument('-t','--threads', type = int, default = 5, help = '[OPTIONAL] The number of threads to use. By default 5')
	parser.add_argument('-o','--output_path', type = str, default = file_path + "/output/", help = '[OPTIONAL] Path to the job folder. If not specified, the job folder will be stored in the default output folder')

	# native files
	parser.add_argument('-n','--binary_interactome', type = str, default = file_path + "/metadata/HomoSapiens_binary_HINThq.txt", help = '[OPTIONAL] A text file with protein-protein interactions. Each row stands for an interaction (IDs should be separated by tab. Only UniProt IDs are accepted). Please include all existing interactions in your context. By default we use the human binary interactome from HINT (http://hint.yulab.org/download/HomoSapiens/binary/hq/)')
	parser.add_argument('-P','--id_mapping', type = str, default = file_path + "/metadata/HUMAN_9606_idmapping.dat.gz", help = '[OPTIONAL] ID conversion file (native file)')
	
	# default settings
	parser.add_argument('-R','--resolution', type = str, default = "low", help = '[OPTIONAL] high or low. If set to high, smaller subnetworks are favored. If set to low, larger subnetworks are favored. By default low')
	parser.add_argument('-r','--restart_prob', type = float, default = 0.2, help = '[OPTIONAL] Restart probability in the diffusion model. By default 0.2')
	parser.add_argument('-w','--weighted_edge', type = float, default = 5.0, help = '[OPTIONAL] Fold increase of heat transfer rate. By default 5.0')
	parser.add_argument('-x','--max_subnetwork_size', type = int, default = 5, help = '[OPTIONAL] Interconnected module size cutoff for determining appropriate delta (default: 5)')
	parser.add_argument('-D','--delta_trial', type = int, default = 20, help = '[OPTIONAL] Number of permutations for determining the value of delta (default: 20)')
	parser.add_argument('-S','--random_seed', type = int, default = 2062, help = '[OPTIONAL] Seed for randomizing heat diffusion process (default: 2062)')
	args = parser.parse_args()


	############################################################################################################################################
	# Check input files
	# -----------------
	os.makedirs(args.output_path, exist_ok = True)
	if path.exists(args.gene_score) == False:
		sys.exit("No such file or directory: " + args.gene_score)
	if path.exists(args.binary_interactome) == False:
		sys.exit("No such file or directory: " + args.binary_interactome)
	if path.exists(args.id_mapping) == False:
		sys.exit("No such file or directory: " + args.id_mapping)
	if args.disrupted_interaction != None and path.exists(args.disrupted_interaction) == False:
		sys.exit("No such file or directory: " + args.disrupted_interaction)
	if args.expressed_genes != None and path.exists(args.expressed_genes) == False:
		sys.exit("No such file or directory: " + args.expressed_genes)
	if args.resolution not in {"high", "low"}:
		sys.exit("Please check --resolution! Can only be set to high or low.")
	############################################################################################################################################

	
	############################################################################################################################################
	# Create the job
	# --------------
	if args.output_path[-1] == "/":
		output_path = args.output_path + args.job_name + "/"
	else:
		output_path = args.output_path + "/" + args.job_name + "/"
	if os.path.exists(output_path) and os.listdir(output_path):
		sys.exit("Output job directory already exists and is not empty: " + output_path)
	os.makedirs(output_path, exist_ok = True)
	############################################################################################################################################


	############################################################################################################################################
	# Print job info
	# --------------
	print("Job name: " + args.job_name)
	print("Output path: " + output_path)
	print("Gene score file: " + args.gene_score)
	if args.disrupted_interaction != None:
		print("Disrupted interaction file: " + args.disrupted_interaction)
	else:
		print("No disrupted interaction file provided.")
	if args.expressed_genes != None:
		print("Expressed genes file: " + args.expressed_genes)
	else:
		print("No expressed genes file provided.")
	print("Number of threads: " + str(args.threads))
	print("Binary interactome file: " + args.binary_interactome)
	print("ID mapping file: " + args.id_mapping)
	print("Resolution: " + args.resolution)
	print("Restart probability: " + str(args.restart_prob))
	print("Weighted edge fold increase: " + str(args.weighted_edge))
	print("Max subnetwork size: " + str(args.max_subnetwork_size))
	print("Delta trial number: " + str(args.delta_trial))
	print("Random seed: " + str(args.random_seed))
	############################################################################################################################################

	
	############################################################################################################################################
	# Load metadata
	# -----------------------------------------
	df_annot = pd.read_csv(args.id_mapping, sep = "\t", dtype = str, header = None)
	df_annot[0] = df_annot[0].apply(lambda x: x.split("-")[0])
	df_annot[2] = df_annot[2].apply(lambda x: x.split(".")[0])
	df_annot = df_annot[[0, 2]].drop_duplicates().rename(columns = {0: "UniProt"})
	############################################################################################################################################


	############################################################################################################################################
	# Build protein-protein interaction network
	# -----------------------------------------
	try:
		df_interactome = pd.read_csv(args.binary_interactome, sep = "\t", header = None)
		pro1_list = df_interactome[0].tolist()
		pro2_list = df_interactome[1].tolist()
	except:
		funcs.remove_whole_dir(output_path)
		sys.exit("Please check the file with binary interactome! Each row must contain two UniProt IDs separated by tab. See example for details.")
	edges = set()
	for i in range(df_interactome.shape[0]):
		pro1 = pro1_list[i].split("-")[0]
		pro2 = pro2_list[i].split("-")[0]
		edges = edges.union({tuple(sorted([pro1, pro2]))})
	G = nx.Graph(sorted(edges))
	for x in G.nodes():
		G.nodes[x]["heat_score"] = 0
	for u,v in G.edges():
		G[u][v]["weight"] = 1
	############################################################################################################################################

	
	############################################################################################################################################
	# Filter out mutations in unexpressed genes and remove unexpressed proteins from the network
	# ------------------------------------------------------------------------------------------
	if args.expressed_genes != None:
		try:
			expr_genes = {x.split(".")[0] for x in pd.read_csv(args.expressed_genes, sep = "\t", header = None)[0]}
		except:
			funcs.remove_whole_dir(output_path)
			sys.exit("Please check the format of the file with expressed genes! Each row should only contain one gene_symbol/ENSG/ENST ID. See example for details.")
		expr_uniprots = set(df_annot[df_annot[2].isin(expr_genes)]["UniProt"])
		G = nx.Graph(G.subgraph(expr_uniprots))
	############################################################################################################################################


	############################################################################################################################################
	# Define heat scores for proteins
	# ---------------
	df = pd.read_csv(args.gene_score, sep = "\t")
	if not {"Gene", "Score"}.issubset(df.columns):
		funcs.remove_whole_dir(output_path)
		sys.exit("Gene score file must contain columns named Gene and Score.")
	df["Score"] = pd.to_numeric(df["Score"], errors = "coerce")
	if df["Score"].isna().any() or not np.isfinite(df["Score"]).all():
		funcs.remove_whole_dir(output_path)
		sys.exit("Every Score value must be finite and numeric.")
	df = pd.merge(df, df_annot, left_on = "Gene", right_on = 2)
	if df.empty:
		funcs.remove_whole_dir(output_path)
		sys.exit("None of the input gene identifiers could be mapped to UniProt IDs.")
	for x, df_one in df.groupby("UniProt"):
		if x in G:
			G.nodes[x]["heat_score"] = df_one["Score"].max()
	############################################################################################################################################


	############################################################################################################################################
	# Assign edge weights to the protein-protein interaction network
	# --------------------------------------------------------------
	if args.disrupted_interaction != None:
		df_disppi = pd.read_csv(args.disrupted_interaction, sep = "\t")
		pro1 = df_disppi["P1"].tolist()
		pro2 = df_disppi["P2"].tolist()
		for i in range(df_disppi.shape[0]):
			# check why doing this???
			if G.has_edge(pro1[i], pro2[i]) and G.nodes[pro1[i]]["heat_score"] > 0 and G.nodes[pro2[i]]["heat_score"] > 0:
				G[pro1[i]][pro2[i]]["weight"] = args.weighted_edge
	############################################################################################################################################


	############################################################################################################################################
	# Network diffusion
	# -----------------
	nx.write_graphml(G, output_path + "initial_state.graphml.gz")
	print("Initial network state saved to " + output_path + "initial_state.graphml.gz")
	
	# find appropriate delta
	if len(G.edges()) < args.max_subnetwork_size:
		delta = -np.inf
	else:
		f = open(output_path + "choose_delta.txt", "w")
		f.write("\t".join(["delta", "top_edge_cutoff", "subnetwork_sizes"]) + "\n")
		f.close()
		function_input = []
		for i in range(args.delta_trial):
			function_input.append([G, output_path + "choose_delta.txt", args.restart_prob, args.max_subnetwork_size, args.random_seed + i])
		p = Pool(args.threads)
		output = p.map(get_one_delta_multirun, function_input)
		p.close()
		if args.resolution == "low":
			delta = pd.read_csv(output_path + "choose_delta.txt", sep = "\t")["delta"].min()
		else:
			delta = pd.read_csv(output_path + "choose_delta.txt", sep = "\t")["delta"].median()
	print("The chosen delta is " + str(delta) + ".")

	# identify subnetworks
	all_subnetworks = funcs.identify_hot_modules(G, output_path + "final_state.graphml.gz", beta = args.restart_prob, delta = delta)
	f = open(output_path + "subnetworks.txt", "w")
	f.write("\t".join(["Subnetwork_genes", "Subnetwork_size"]) + "\n")
	for x in sorted(all_subnetworks, key = lambda z: len(z), reverse = True):
		f.write(",".join(sorted(x)) + "\t" + str(len(x)) + "\n")
	f.close()
	sizes = [len(x) for x in all_subnetworks]
	if sizes:
		print(str(len(all_subnetworks)) + " subnetworks are identified, whose sizes range from " + str(min(sizes)) + " to " + str(max(sizes)) + ".")
	else:
		print("0 subnetworks were identified.")
	############################################################################################################################################
	
	print("Job finished!")
