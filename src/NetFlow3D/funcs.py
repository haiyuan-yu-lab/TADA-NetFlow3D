description='''
- Functions used by NetFlow3D_ASD_version.py
- Dependencies: None
'''

import numpy as np
import networkx as nx
import copy
import shutil


def get_heat_diffusion_matrix(G, beta):
	'''
	Get the diffusion matrix of a network. 
	G: undirected graph (no isolated nodes)
	beta: restart probability / the fraction of retained heat 
	return:
	  --F: heat diffusion matrix
	  --all_nodes: a list of nodes corresponding to the index of heat diffusion matrix
	'''
	all_nodes = list(G.nodes())
	adjacency_matrix = nx.adjacency_matrix(G).todense()
	W = adjacency_matrix / np.sum(adjacency_matrix, axis = 0)
	F = beta * np.linalg.inv(np.identity(len(all_nodes)) - (1 - beta) * W)
	return [F, all_nodes]


def identify_hot_modules(G, graph_output, beta, delta):
	"""
	Parameters
	----------
	G: undirected graph. Each node has an attribute **heat_score**. heat_score = 0 means that there is no heat from that node.
	beta: restart probability
	delta: minimum amount of heat diffused from one node to another so that a directed edge will be added
	
	Return
	------
	a list of sets, the nodes in each set belong to the same hot module
	"""
	def get_hot_subnetworks(F, nodes, Dh, delta):
		"""
		F: heat diffusion matrix
		nodes: a list of nodes corresponding to the index of heat diffusion matrix
		Dh: diagonal matrix with the heat scores of **nodes**
		return: a list of sets, the nodes in each set belong to the same hot subnetwork (min size of a hot subnetwork is 2)
		"""
		idx2uniprot = {}
		for i in range(len(nodes)):
			idx2uniprot[i] = nodes[i]
		E = np.matmul(F, Dh)
		E[E <= delta] = 0
		G = nx.from_numpy_array(E.transpose(), create_using = nx.DiGraph)
		G = nx.relabel_nodes(G, idx2uniprot)
		return G

	all_modules = []
	G_final = nx.DiGraph()
	for element in nx.connected_components(G):
		heat_scores = [G.nodes[item]["heat_score"] for item in element]
		if max(heat_scores) > 0:
			if len(element) > 1:
				G_sub = G.subgraph(element)
				F, nodes = get_heat_diffusion_matrix(G_sub, beta = beta)
				diag_elements = [G_sub.nodes[item]["heat_score"] for item in nodes]
				G_one = get_hot_subnetworks(F, nodes, np.diag(diag_elements), delta = delta)
				G_final.add_nodes_from(G_one.nodes(data = True))
				G_final.add_edges_from(G_one.edges(data = True))
				for x in nx.strongly_connected_components(G_one):
					if len(x) > 1 or G.nodes[list(x)[0]]["heat_score"] > 0:
						all_modules.append(x)
			else:
				all_modules.append(element)
	nx.write_graphml(G_final, graph_output)
	return all_modules


def format_graph(G):
	"""
	Sort the nodes and edges in a graph.
	"""
	H = nx.Graph()
	H.add_nodes_from(sorted(G.nodes(data = True)))
	H.add_edges_from(sorted([(sorted([u,v])[0], sorted([u,v])[1], info) for u,v,info in G.edges(data = True)]))
	return H


def get_suitable_delta(G, beta, size_cutoff, output, upper_bound = 1, lower_bound = 0):
	"""
	Purpose
	--------- 
	- To find the appropriate delta for a network with labeled heat sources
	- delta: minimum amount of heat diffused from one node to another so that a directed edge will be added
	- "Appropriate" means that using this delta value there won't be giant modules. The maximum module size would be **size_cutoff**
	
	Parameters
	----------
	G: undirected graph. Each node has an attribute **heat_score**. heat_score = 0 means that there is no heat placed on that node.
	beta: restart probability / the fraction of retained heat
	size_cutoff: the maximum size of a hot module
	

	Returns
	-------
	final_delta: the appropriate delta value
	final_size: the list of hot module sizes using **final_delta** as cutoff
	"""
	# get edge values 
	heat_edges = []
	E_list = []
	for element in nx.connected_components(G):
		heat_scores = [G.nodes[item]["heat_score"] for item in element]
		if len(element) > 1 and max(heat_scores) > 0:
			G_sub = G.subgraph(element)
			F, nodes = get_heat_diffusion_matrix(G_sub, beta = beta)
			diag_elements = [G_sub.nodes[item]["heat_score"] for item in nodes]
			E = np.matmul(F, np.diag(diag_elements))
			E_list.append(E)
			# heat_edges.extend(E[E.nonzero()].tolist())
			# heat_edges.extend(E[E > 0].flatten().tolist())
			heat_edges.extend(E[E.nonzero()].ravel().tolist())
	heat_edges = sorted(heat_edges, reverse = True)

	# get the maximum difference between upper_bound and lower_bound
	diff = max(2 / len(heat_edges), 0.00001)

	# find appropriate delta
	module_size = []
	for E_copy in copy.deepcopy(E_list):
		G = nx.from_numpy_array(E_copy.transpose(), create_using = nx.DiGraph)
		for element in list(nx.strongly_connected_components(G)):
			module_size.append(len(element))
	if max(module_size) <= size_cutoff:
		final_delta = -1
		final_midpoint = 1
	else:
		while upper_bound - lower_bound > diff:
			midpoint = (upper_bound + lower_bound) / 2
			delta = heat_edges[int(midpoint * len(heat_edges))]
			module_size = []
			for E_copy in copy.deepcopy(E_list):
				E_copy[E_copy <= delta] = 0
				G = nx.from_numpy_array(E_copy.transpose(), create_using = nx.DiGraph)
				for element in list(nx.strongly_connected_components(G)):
					module_size.append(len(element))
			if max(module_size) <= size_cutoff:
				lower_bound = midpoint
			else:
				upper_bound = midpoint
		final_midpoint = lower_bound
		final_delta = heat_edges[int(final_midpoint * len(heat_edges)) + 1]

	# write output
	f = open(output, "a+")
	f.write("\t".join([str(final_delta), str(final_midpoint), ",".join([str(thing) for thing in module_size])]) + "\n")
	f.close()
	
	return


def get_one_delta(G_original, delta_output, restart_prob, max_subnetwork_size, seed):
	G = copy.deepcopy(G_original)
	nx.double_edge_swap(G, nswap = len(G.edges()), max_tries = min(5 * len(G.edges()), 5000000), seed = seed)
	for u,v in G.edges():
		G[u][v]["weight"] = 1.0
	get_suitable_delta(G, beta = restart_prob, size_cutoff = max_subnetwork_size, output = delta_output)
	return


def remove_whole_dir(dir_name):
	shutil.rmtree(dir_name)
	return


