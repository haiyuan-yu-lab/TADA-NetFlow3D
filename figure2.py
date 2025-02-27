import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.path import Path
import matplotlib.patches as patches
import pandas as pd
import numpy as np
import networkx as nx
from matplotlib.patches import Rectangle
from collections import defaultdict
from itertools import combinations
import math
import matplotlib
import sys
from scipy.integrate import quad
sys.path.append("/home/yz2296/bin")
import all_function_py3

def quadratic_bezier_curves(t, p0, p1, p2):
	'''
	p0: e.g. (0, 0)
	p1: e.g. (1, 1)
	p2: e.g. (2, 0)
	'''
	x = (1 - t) * ((1 - t) * p0[0] + t * p1[0]) + t * ((1 - t) * p1[0] + t * p2[0])
	y = (1 - t) * ((1 - t) * p0[1] + t * p1[1]) + t * ((1 - t) * p1[1] + t * p2[1])
	return (x, y)


def quadratic_bezier_length(P0, P1, P2):
    def integrand(t):
        dx_dt = 2 * (1 - t) * (P1[0] - P0[0]) + 2 * t * (P2[0] - P1[0])
        dy_dt = 2 * (1 - t) * (P1[1] - P0[1]) + 2 * t * (P2[1] - P1[1])
        return np.sqrt(dx_dt**2 + dy_dt**2)
        
    length, _ = quad(integrand, 0, 1)
    return length


def circos_plot(cluster_file, G_PPI, output = None, ordered_genes = None, donut_width = 0.6, donut_radius = 5, edge_width = 0.05, figsize = (30,30), \
	fontsize = 4, text_status = True, text_padding = 1.004, homodimer_padding = 1.0084, figure_padding = 1.02, offcenter = 0, mis_frac = 0.5, \
	highlight_genes = set(), layer_count = 1, layer_zoom = 2):
	'''
	cluster_file: /fs/cbsuhy01/storage/yz2296/cancer_hotspot_new/website/website/download/pancan_clusters.txt
	output: /fs/cbsuhy01/storage/yz2296/cancer_hotspot_new/website/figures/figure_2/circos_plot.pdf
	offcenter: to what extent are the edges away from center
	'''
	# load cluster information
	df = pd.read_csv(cluster_file, sep = "\t", dtype = str)
	if ordered_genes != None:
		df = df[df["Affected_genes"].apply(lambda x: set(x.split(",")) - set(ordered_genes) == set())]
	
	# group the genes
	df_lof = df[df["Type"] == "LoF_IntraProtein_Enriched"]
	df_intra = df[df["Type"] == "Missense_IntraProtein_Cluster"]
	df_inter = df[df["Type"] == "Missense_InterProtein_Cluster"]
	lof_genes = set(df_lof["Affected_genes"].tolist())
	intra_genes = set(df_intra["Affected_genes"].tolist())
	inter_genes = set()
	homodimer_genes = set()
	for x in df_inter["Affected_genes"]:
		inter_genes = inter_genes.union(set(x.split(",")))
		if len(set(x.split(","))) == 1:
			homodimer_genes = homodimer_genes.union(set(x.split(",")))
	
	# sort the genes
	if ordered_genes == None:
		part1 = sorted(lof_genes - intra_genes - inter_genes) + sorted((lof_genes - intra_genes).intersection(inter_genes))
		part21 = sorted(lof_genes.intersection(intra_genes) - inter_genes)
		part22 = sorted(lof_genes.intersection(intra_genes).intersection(inter_genes))
		part3 = sorted((intra_genes - lof_genes).intersection(inter_genes)) + sorted(intra_genes - lof_genes - inter_genes)
		ordered_genes = part1 + part22[0:int(len(part22) / 2)] + part21 + part22[int(len(part22) / 2):] + part3 + sorted(inter_genes - intra_genes - lof_genes)
	
	# get the annotation for each gene
	gene2label = {}
	for i, x in enumerate(ordered_genes):
		label = []
		if x in lof_genes:
			label.append(True)
		else:
			label.append(False)
		if x in intra_genes:
			label.append(True)
		else:
			label.append(False)
		label.extend([0.25 + i / len(ordered_genes), 0.25 + (i + 1) / len(ordered_genes)])
		gene2label[x] = label

	# plot sectors
	fig,ax = plt.subplots(figsize = figsize)
	ax.set_aspect(1)
	gap = 1 / len(ordered_genes) * 0.1
	max_string_len = max([len(x) for x in ordered_genes])
	for gene in gene2label:
		lof, mis, start_angle, end_angle = gene2label[gene]
		if lof == True:
			fov = Wedge((0,0), donut_radius, (start_angle + gap) * 360, (end_angle - gap) * 360, facecolor = "#fed9a6", width = donut_width, \
				edgecolor = "white", linewidth = 0)
			ax.add_artist(fov)
		else:
			fov = Wedge((0,0), donut_radius, (start_angle + gap) * 360, (end_angle - gap) * 360, facecolor = "#252525", width = donut_width, \
				edgecolor = "white", linewidth = 0)
			ax.add_artist(fov)
		if mis == True:
			fov = Wedge((0,0), donut_radius - 0.5 * (1 - mis_frac) * donut_width, \
				(start_angle + gap + 0.5 * (1 - mis_frac) * (end_angle - start_angle - 2 * gap)) * 360, \
				(end_angle - gap - 0.5 * (1 - mis_frac) * (end_angle - start_angle - 2 * gap)) * 360, facecolor = "#fb8072", \
				width = mis_frac * donut_width, edgecolor = "white", linewidth = 0)
			# fov = Wedge((0,0), donut_radius - 1 / 3 * donut_width, (start_angle + gap) * 360, (end_angle - gap) * 360, facecolor = "#fb8072", \
			# 	width = 1 / 3 * donut_width, edgecolor = "white", linewidth = 0)
			ax.add_artist(fov)
		text_angle = (start_angle + end_angle) / 2
		if gene in highlight_genes:
			color = "#80b1d3"
			fontweight = "heavy"
		else:
			color = "#252525"
			fontweight = "light"
		if text_angle < 0.75:
			rotation = text_angle * 360 + 180
			gene = gene.rjust(max_string_len)
		else:
			rotation = text_angle * 360
			gene = gene.ljust(max_string_len)
		if text_status == True:
			plt.text(donut_radius * text_padding * np.cos(text_angle * 2 * np.pi), donut_radius * text_padding * np.sin(text_angle * 2 * np.pi), gene, \
				horizontalalignment = "center", verticalalignment = "center", rotation = rotation, fontsize = fontsize, color = color, \
				fontweight = fontweight) # r = 5.02 for figure 2
	
	# plot lines
	G_sub = G_PPI.subgraph(ordered_genes)
	# codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
	codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
	for x,y,info in G_sub.edges(data = True):
		_, _, start1, end1 = gene2label[x]
		_, _, start2, end2 = gene2label[y]
		angle1 = (start1 + end1) / 2
		angle2 = (start2 + end2) / 2
		angle = (angle1 + angle2) / 2
		if x != y:
			diff = abs(ordered_genes.index(x) - ordered_genes.index(y))
			if diff > int(0.5 * len(ordered_genes)):
				smaller, larger = sorted([ordered_genes.index(x), ordered_genes.index(y)])
				diff = smaller + len(ordered_genes) - larger
				angle = angle + 0.5
			r_marker = 0.5 * 1 / diff * donut_radius
			r_edge = donut_radius - donut_width * 0.96
		else:
			r_marker = donut_radius * 1.2
			r_edge = donut_radius - donut_width * 0.1
		x_marker = r_marker * np.cos(angle * 2 * np.pi)
		y_marker = r_marker * np.sin(angle * 2 * np.pi)
		if x != y:
			verts = [(r_edge * np.cos(angle1 * 2 * np.pi), r_edge * np.sin(angle1 * 2 * np.pi)), # P0
					 (x_marker, y_marker),  # P1
					 # (x_marker, y_marker),  # P2
					 (r_edge * np.cos(angle2 * 2 * np.pi), r_edge * np.sin(angle2 * 2 * np.pi))]  # P3
		else:
			verts = [(r_edge * np.cos((start1 + 2 * gap) * 2 * np.pi), r_edge * np.sin((start1 + 2 * gap) * 2 * np.pi)), # P0
					 (x_marker, y_marker),  # P1
					 # (x_marker, y_marker),  # P2
					 (r_edge * np.cos((end1 - 2 * gap) * 2 * np.pi), r_edge * np.sin((end1 - 2 * gap) * 2 * np.pi))]  # P3
		if info["weight"] > 1:
			ew = edge_width
			color = "#fb8072"
			zorder = 1
		else:
			ew = edge_width / 5
			color = "#737373" #"#bdbdbd"
			zorder = 0
		path = Path(verts, codes)
		patch = patches.PathPatch(path, facecolor = 'none', lw = ew, edgecolor = color, zorder = zorder)
		ax.add_patch(patch)

	# step = 0.05
	# for i in range(layer_count):
	# 	circle = plt.Circle((0, 0), donut_radius * layer_zoom * (1 + i * step), fill = False)
	# 	ax.add_artist(circle)
			
	# for x in homodimer_genes:
	# 	_, _, start, end = gene2label[x]
	# 	angle = (start + end) / 2
	# 	# plt.scatter(donut_radius * homodimer_padding * np.cos(angle * 2 * np.pi), donut_radius * homodimer_padding * np.sin(angle * 2 * np.pi), \
	# 	# 	c = "#fb8072", s = fontsize, marker = "o", linewidths = 0)
	# 	circle = plt.Circle((donut_radius * homodimer_padding * np.cos(angle * 2 * np.pi), donut_radius * homodimer_padding * np.sin(angle * 2 * np.pi)), \
	# 		0.75 * donut_radius * (end - start), color = '#fb8072')
	# 	ax.add_patch(circle)
	
	# ax.set_xlim(-donut_radius * layer_zoom * figure_padding, donut_radius * layer_zoom * figure_padding)
	# ax.set_ylim(-donut_radius * layer_zoom * figure_padding, donut_radius * layer_zoom * figure_padding)
	ax.set_xlim(-donut_radius * figure_padding, donut_radius * figure_padding)
	ax.set_ylim(-donut_radius * figure_padding, donut_radius * figure_padding)
	plt.axis('off')
	if output != None:
		fig.savefig(output, bbox_inches='tight')
	plt.show()
	return




def shell_layout_ASD(G, cluster_marker, ax, x_range, y_range, ordered_genes = None, interval_ratio = 1, \
	edge_magnifier = 1, fontsize = 16, figure_padding = 0.2, tada_genes = set(), tadawgs_genes = set(), \
	center = (0,0), marker_zoom = 0.05, title = None, title_pos = "top", title_fontsize = 50):
	'''
	G: a module whose nodes have attribute **secondary**, and whose edges have attribute **weight**
	cluster_marker: the image representing a cluster
	tada_genes: TADA genes
	tadawgs_genes: genes identified by TADA-WGS
	interval_ratio: 1 means one node could be placed between two neighboring nodes.
	title_pos: 1. top; 2. bottom
	'''
	# get maxiumum edge weight
	max_string_len = max([len(x) for x in ordered_genes])
	max_edge = max([info["weight"] for u,v,info in G.edges(data = True)])
	gene2label = {}
	for i, x in enumerate(ordered_genes):
		gene2label[x] = (0.25 + (i + 0.5) / len(ordered_genes)) % 1
	# cmap = mcolors.LinearSegmentedColormap.from_list("", ["white", "red"])
	# norm = mcolors.Normalize(vmin=0, vmax=10)
	# colors = [cmap(norm(value)) for value in values]

	# plot nodes
	perimeter = (1 + interval_ratio) * len(ordered_genes)
	shell_radius = perimeter / (2 * np.pi)
	r_text = shell_radius + 0.5 + fontsize / 70
	self_loop_weight = []
	for u,v,info in G.edges(data = True):
		if u == v:
			self_loop_weight.append(info["weight"])
	if len(self_loop_weight) > 0:
		r_text = r_text + 0.5 + max(self_loop_weight) * edge_magnifier / 50
	ax.set_aspect(1)
	for gene in gene2label:
		node_color = "black"
		# if gene in tada_genes: # original TADA 能identify到的
		# 	node_color = "#de2d26"
		# elif gene in tadawgs_genes: # 只有TADA-WGS 能identify到的
		# 	node_color = "#fc9272"
		# else:
		# 	node_color = "#fee0d2" # TADA-WGS也identify不到只能用network identify到的
		angle = gene2label[gene]
		x = center[0] + shell_radius * np.cos(angle * 2 * np.pi)
		y = center[1] + shell_radius * np.sin(angle * 2 * np.pi)
		circle = plt.Circle((x, y), 0.5, color = node_color)
		ax.add_artist(circle)
		if angle < 0.75 and angle > 0.25:
			rotation = angle * 360 + 90 + 180
		else:
			rotation = angle * 360 - 90
		plt.text(center[0] + r_text * np.cos(angle * 2 * np.pi), center[1] + r_text * np.sin(angle * 2 * np.pi), gene, horizontalalignment = "center", \
			verticalalignment = "center", rotation = rotation, fontsize = fontsize, color = "#252525", fontweight = "light")
	
	# plot lines
	codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
	for x,y,info in G.edges(data = True):
		angle1 = gene2label[x]
		angle2 = gene2label[y]
		angle = (angle1 + angle2) / 2
		if x != y:
			smaller, larger = sorted([ordered_genes.index(x), ordered_genes.index(y)])
			if larger - smaller > 0.5 * len(ordered_genes):
				diff = smaller + len(ordered_genes) - larger
			else:
				diff = larger - smaller
			r_marker = 0.5 / diff * shell_radius
			if abs(angle1 - angle2) > 0.5:
				angle = angle - 0.5
			x_marker = center[0] + r_marker * np.cos(angle * 2 * np.pi)
			y_marker = center[1] + r_marker * np.sin(angle * 2 * np.pi)
			verts = [(center[0] + shell_radius * np.cos(angle1 * 2 * np.pi), center[1] + shell_radius * np.sin(angle1 * 2 * np.pi)), # P0
					 (x_marker, y_marker),  # P1
					 (center[0] + shell_radius * np.cos(angle2 * 2 * np.pi), center[1] + shell_radius * np.sin(angle2 * 2 * np.pi))]  # P2
			path = Path(verts, codes)
			patch = patches.PathPatch(path, facecolor = 'none', lw = info["weight"] * edge_magnifier, edgecolor = "#252525", zorder = 0)
			ax.add_patch(patch)
			if info["weight"] > 1:
				# length = quadratic_bezier_length(verts[0], verts[1], verts[2])
				# x_vertex, y_vertex = quadratic_bezier_curves(0.5 / length, verts[0], verts[1], verts[2])
				# all_function_py3.display_image(cluster_marker, ax, x_vertex, y_vertex, zoom = marker_zoom)
				# x_vertex, y_vertex = quadratic_bezier_curves(1 - 0.5 / length, verts[0], verts[1], verts[2])
				# all_function_py3.display_image(cluster_marker, ax, x_vertex, y_vertex, zoom = marker_zoom)
				x_vertex, y_vertex = quadratic_bezier_curves(0.5, verts[0], verts[1], verts[2])
				all_function_py3.display_image(cluster_marker, ax, x_vertex, y_vertex, zoom = marker_zoom)
		else:
			circle = plt.Circle((center[0] + (shell_radius + 0.65) * np.cos(angle1 * 2 * np.pi), \
				center[1] + (shell_radius + 0.65) * np.sin(angle1 * 2 * np.pi)), 0.35, \
				edgecolor = "#252525", linewidth = info["weight"] * edge_magnifier, zorder = 0, fill = False)
			ax.add_artist(circle)
			if info["weight"] > 1:
				all_function_py3.display_image(cluster_marker, ax, center[0] + (shell_radius + 1) * np.cos(angle1 * 2 * np.pi), \
					center[1] + (shell_radius + 1) * np.sin(angle1 * 2 * np.pi), zoom = marker_zoom)

	if title != None:
		if title_pos == "top":
			y_title = center[1] + r_text + fontsize / 70 + 0.3
			va = "bottom"
		else:
			y_title = center[1] - r_text - 2 * fontsize / 70 - 0.3
			va = "top"
		plt.text(center[0], y_title, title, horizontalalignment = "center", verticalalignment = va, fontsize = title_fontsize, color = "#252525", \
			fontweight = "normal")

	ax.set_xlim(x_range[0], x_range[1])
	ax.set_ylim(y_range[0], y_range[1])
	plt.axis('off')
	return ax


def example_PPI(cluster_file, G, cluster_marker, lof_marker, ax, x_range, y_range, ordered_genes, gene2color, interval_ratio = 1, \
	edge_magnifier = 1, fontsize = 16, figure_padding = 0.2, center = (0,0), gene2secondary = {}, marker_zoom = 0.05, title = None, \
	title_pos = "top", title_fontsize = 50, with_labels = True):
	'''
	cluster_file: /fs/cbsuhy01/storage/yz2296/cancer_hotspot_new/website/website/download/pancan_clusters.txt
	G: a module whose nodes have attribute **secondary**, and whose edges have attribute **weight**
	cluster_marker: the image representing a cluster
	interval_ratio: 1 means one node could be placed between two neighboring nodes.
	title_pos: 1. top; 2. bottom
	'''
	# load cluster information
	df = pd.read_csv(cluster_file, sep = "\t", dtype = str)

	# get the ordered nodes
	max_string_len = max([len(x) for x in ordered_genes])

	# get maxiumum edge weight
	max_edge = max([info["weight"] for u,v,info in G.edges(data = True)])

	# get the annotations for each gene
	df = df[df["Affected_genes"].apply(lambda x: x in ordered_genes)]
	lof_genes = set(df[df["Type"] == "LoF_IntraProtein_Enriched"]["Affected_genes"])
	intra_genes = set(df[df["Type"] == "NonTrunc_IntraProtein_Cluster"]["Affected_genes"])
	gene2label = {}
	for i, x in enumerate(ordered_genes):
		label = []
		if x in lof_genes:
			label.append(True)
		else:
			label.append(False)
		if x in intra_genes:
			label.append(True)
		else:
			label.append(False)
		label.append((0.25 + (i + 0.5) / len(ordered_genes)) % 1)
		gene2label[x] = label

	# plot nodes
	perimeter = (1 + interval_ratio) * len(ordered_genes)
	shell_radius = perimeter / (2 * np.pi)
	r_text = shell_radius + 0.5 + fontsize / 70
	self_loop_weight = []
	for u,v,info in G.edges(data = True):
		if u == v:
			self_loop_weight.append(info["weight"])
	if len(self_loop_weight) > 0:
		r_text = r_text + 0.5 + max(self_loop_weight) * edge_magnifier / 50
	ax.set_aspect(1)
	for gene in gene2label:
		node_color = gene2color[gene]
		lof, mis, angle = gene2label[gene]
		x = center[0] + shell_radius * np.cos(angle * 2 * np.pi)
		y = center[1] + shell_radius * np.sin(angle * 2 * np.pi)
		circle = plt.Circle((x, y), 0.5, color = node_color)
		ax.add_artist(circle)
		if mis == True and lof == False:
			all_function_py3.display_image(cluster_marker, ax, x, y, zoom = marker_zoom)
		elif lof == True and mis == False:
			all_function_py3.display_image(lof_marker, ax, x, y, zoom = marker_zoom)
		elif mis == True and lof == True:
			x0 = center[0] + (shell_radius + 0.18) * np.cos(angle * 2 * np.pi)
			y0 = center[1] + (shell_radius + 0.18) * np.sin(angle * 2 * np.pi)
			all_function_py3.display_image(cluster_marker, ax, x0, y0, zoom = marker_zoom)
			x1 = center[0] + (shell_radius - 0.29) * np.cos(angle * 2 * np.pi)
			y1 = center[1] + (shell_radius - 0.29) * np.sin(angle * 2 * np.pi)
			all_function_py3.display_image(lof_marker, ax, x1, y1, zoom = marker_zoom)
			
		if gene in gene2secondary:
			x_add = (1.5 - 0.22) * np.cos((angle - 0.16) * 2 * np.pi)
			y_add = (1.5 - 0.22) * np.sin((angle - 0.16) * 2 * np.pi)
			plt.plot([x, x + x_add], [y, y + y_add], linewidth = edge_magnifier, color = "#bdbdbd", zorder = 0, linestyle = "dashed")
			circle = plt.Circle((x + 1.5 * np.cos((angle - 0.16) * 2 * np.pi), y + 1.5 * np.sin((angle - 0.16) * 2 * np.pi)), 0.22, \
				edgecolor = "#bdbdbd", linewidth = edge_magnifier, fill = False)
			ax.add_artist(circle)
			plt.text(x + 1.5 * np.cos((angle - 0.16) * 2 * np.pi), y + 1.5 * np.sin((angle - 0.16) * 2 * np.pi), str(len(gene2secondary[gene])), \
				horizontalalignment = "center", verticalalignment = "center", \
				fontsize = int(fontsize * 0.7), color = "#bdbdbd")
		if angle < 0.75 and angle > 0.25:
			rotation = angle * 360 + 90 + 180
		else:
			rotation = angle * 360 - 90
		if with_labels == True:
			plt.text(center[0] + r_text * np.cos(angle * 2 * np.pi), center[1] + r_text * np.sin(angle * 2 * np.pi), gene, horizontalalignment = "center", \
				verticalalignment = "center", rotation = rotation, fontsize = fontsize, color = "#252525", fontweight = "light")
	
	# plot lines
	codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
	for x,y,info in G.edges(data = True):
		_, _, angle1 = gene2label[x]
		_, _, angle2 = gene2label[y]
		angle = (angle1 + angle2) / 2
		if x != y:
			smaller, larger = sorted([ordered_genes.index(x), ordered_genes.index(y)])
			if larger - smaller > 0.5 * len(ordered_genes):
				diff = smaller + len(ordered_genes) - larger
			else:
				diff = larger - smaller
			r_marker = 0.5 / diff * shell_radius
			if abs(angle1 - angle2) > 0.5:
				angle = angle - 0.5
			x_marker = center[0] + r_marker * np.cos(angle * 2 * np.pi)
			y_marker = center[1] + r_marker * np.sin(angle * 2 * np.pi)
			verts = [(center[0] + shell_radius * np.cos(angle1 * 2 * np.pi), center[1] + shell_radius * np.sin(angle1 * 2 * np.pi)), # P0
					 (x_marker, y_marker),  # P1
					 (center[0] + shell_radius * np.cos(angle2 * 2 * np.pi), center[1] + shell_radius * np.sin(angle2 * 2 * np.pi))]  # P2
			path = Path(verts, codes)
			patch = patches.PathPatch(path, facecolor = 'none', lw = info["weight"] * edge_magnifier, edgecolor = "#252525", zorder = 0)
			ax.add_patch(patch)
			if info["weight"] > 1:
				x_vertex, y_vertex = quadratic_bezier_curves(0.5, verts[0], verts[1], verts[2])
				all_function_py3.display_image(cluster_marker, ax, x_vertex, y_vertex, zoom = marker_zoom)
		else:
			circle = plt.Circle((center[0] + (shell_radius + 0.65) * np.cos(angle1 * 2 * np.pi), \
				center[1] + (shell_radius + 0.65) * np.sin(angle1 * 2 * np.pi)), 0.35, \
				edgecolor = "#252525", linewidth = info["weight"] * edge_magnifier, zorder = 0, fill = False)
			ax.add_artist(circle)
			if info["weight"] > 1:
				all_function_py3.display_image(cluster_marker, ax, center[0] + (shell_radius + 1) * np.cos(angle1 * 2 * np.pi), \
					center[1] + (shell_radius + 1) * np.sin(angle1 * 2 * np.pi), zoom = marker_zoom)

	if title != None:
		if title_pos == "top":
			y_title = center[1] + r_text + fontsize / 70 + 0.3
			va = "bottom"
		else:
			y_title = center[1] - r_text - 2 * fontsize / 70 - 0.3
			va = "top"
		plt.text(center[0], y_title, title, horizontalalignment = "center", verticalalignment = va, fontsize = title_fontsize, color = "#252525", \
			fontweight = "normal")

	ax.set_xlim(x_range[0], x_range[1])
	ax.set_ylim(y_range[0], y_range[1])
	plt.axis('off')
	return ax


def example_HeatDiffusion(G, ax, ordered_genes, gene2color, x_range, y_range, interval_ratio = 1, fontsize = 16, center = (0,0), \
	title = None, title_pos = "top", title_fontsize = 50, specific_edge = {}, edge_radius = 0.25, arrow_size = 6, with_labels = True):
	'''
	cluster_file: /fs/cbsuhy01/storage/yz2296/cancer_hotspot_new/website/website/download/pancan_clusters.txt
	G: a module whose nodes have attribute **secondary**, and whose edges have attribute **weight**
	interval_ratio: 1 means one node could be placed between two neighboring nodes.
	title_pos: 1. top; 2. bottom
	'''
	node_radius = 0.5

	# plot nodes
	perimeter = (1 + interval_ratio) * len(ordered_genes)
	shell_radius = perimeter / (2 * np.pi)
	r_text = shell_radius + 0.5 + fontsize / 70
	ax.set_aspect(1)
	node2pos = {}
	for i, gene in enumerate(ordered_genes):
		angle = (0.25 + (i + 0.5) / len(ordered_genes)) % 1
		x = center[0] + shell_radius * np.cos(angle * 2 * np.pi)
		y = center[1] + shell_radius * np.sin(angle * 2 * np.pi)
		circle = plt.Circle((x, y), node_radius, color = gene2color[gene])
		node2pos[gene] = (x, y)
		ax.add_artist(circle)
		if angle < 0.75 and angle > 0.25:
			rotation = angle * 360 + 90 + 180
		else:
			rotation = angle * 360 - 90
		if with_labels == True:
			plt.text(center[0] + r_text * np.cos(angle * 2 * np.pi), center[1] + r_text * np.sin(angle * 2 * np.pi), gene, horizontalalignment = "center", \
				verticalalignment = "center", rotation = rotation, fontsize = fontsize, color = "#252525", fontweight = "light")

	# plot edges
	for u,v,info in G.edges(data = True):
		if u != v:
			if (u,v) not in specific_edge:
				x1, y1 = node2pos[u]
				x2, y2 = node2pos[v]
				if x1 == x2:
					if y2 > y1:
						theta = 0.25 * 2 * np.pi
					else:
						theta = -0.25 * 2 * np.pi
				else:
					theta = math.atan((y2 - y1) / (x2 - x1))
					if x2 - x1 < 0:
						theta = theta + np.pi
				x1 = x1 + node_radius * 1.1 * math.cos(theta - 0.2)
				y1 = y1 + node_radius * 1.1 * math.sin(theta - 0.2)
				x2 = x2 + node_radius * 1.1 * math.cos(theta + np.pi + 0.2)
				y2 = y2 + node_radius * 1.1 * math.sin(theta + np.pi + 0.2)
				rad = edge_radius
			else:
				pos1, pos2, rad = specific_edge[(u,v)]
				x1, y1 = pos1
				x2, y2 = pos2
			kw = dict(arrowstyle = "Simple, tail_width=0, head_width=" + str(arrow_size) + ", head_length=" + str(arrow_size), color = "red")
			lw = info["weight"]
			a = patches.FancyArrowPatch((x1, y1), (x2, y2), connectionstyle = "arc3,rad=" + str(rad), **kw, lw = lw, zorder = 50)
			ax.add_patch(a)

	if title != None:
		if title_pos == "top":
			y_title = center[1] + r_text + fontsize / 70 + 0.3
			va = "bottom"
		else:
			y_title = center[1] - r_text - 2 * fontsize / 70 - 0.3
			va = "top"
		plt.text(center[0], y_title, title, horizontalalignment = "center", verticalalignment = va, fontsize = title_fontsize, color = "#252525", \
			fontweight = "normal")

	ax.set_xlim(x_range[0], x_range[1])
	ax.set_ylim(y_range[0], y_range[1])
	plt.axis('off')
	return ax


def mutation_distribution_lof(gene, protein_length, maf_file, cluster_file, gene_color = "#a1d99b", output = None, unit_length = 0.05):
	'''
	gene: e.g. "MXD1"
	maf_file: /local/storage/yz2296/cancer_hotspot_new/data/version_5/Mutation_data/Pan_cancer/Zhang_Mutations.txt
	cluster_file: /local/storage/yz2296/cancer_hotspot_new/website/website/download/pancan_clusters.txt
	'''
	def plot_mut(row):
		for i in range(row["Hugo_Symbol"]):
			plt.scatter((int(row["Protein_position"]) - 0.5) * unit_length, 1.55 + i, s = 2000, marker = "X", \
				c = lof_color, linewidth = 2, edgecolors = "k")
		return True

	lof = {"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site"}

	# get clusters where the gene is involved
	df_cluster = pd.read_csv(cluster_file, sep = "\t", dtype = str)
	df_cluster_lof = df_cluster[(df_cluster["Affected_genes"] == gene) & df_cluster["Type"].apply(lambda x: x[0:len("LoF")] == "LoF")]
	print(df_cluster_lof["LoF_enrichment"].tolist()[0], df_cluster_lof["Adjusted_pvalue"].tolist()[0])
	if df_cluster_lof.shape[0] > 0:
		lof_color = "red"
	else:
		lof_color = "#d9d9d9"
		
	# get mutations in the gene
	df_mut = pd.read_csv(maf_file, sep = "\t")
	df_lof = df_mut[df_mut["Variant_Classification"].isin(lof) & (df_mut["Hugo_Symbol"] == gene)].groupby("Protein_position").count().reset_index().sort_values(by = "Protein_position", ascending = False)
	
	# plot mutation distribution
	fig_width = unit_length * protein_length
	fig_height = 16
	fig, ax = plt.subplots(figsize = (fig_width * 1.1, fig_height))
	plt.xlim(-0.05 * fig_width, 1.05 * fig_width)
	plt.ylim(0, fig_height)
	plt.axis('off')
	ax.fill([0, fig_width, fig_width, 0], [0, 0, 1, 1], gene_color)
	df_lof["status"] = df_lof.apply(plot_mut, axis = 1)
	if output != None:
		fig.savefig(output, bbox_inches = "tight", transparent = True)
	return ax


def mutation_distribution_mis_circle(gene, cluster_id, protein_seq, maf_file, cluster_file, gene_color = "#a1d99b", output = None, unitaa_width = 0.1, block_height = 1, \
	unitmut_height = 1, res_size = 100, fontsize = 24):
	'''
	gene: e.g. "MXD1"
	maf_file: /local/storage/yz2296/cancer_hotspot_new/data/version_5/Mutation_data/Pan_cancer/Zhang_Mutations.txt
	cluster_file: /local/storage/yz2296/cancer_hotspot_new/website/website/download/pancan_clusters.txt
	block_width: the width of 10 AAs
	gene_color: the value could be str or dict. If dict, keys are protein positions and values are the corresponding colors (e.g. {1: "blue", 2: "green", ...})
	'''
	mis = {"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"}

	# get clusters where the gene is involved
	df_cluster = pd.read_csv(cluster_file, sep = "\t")
	driver_res = set()
	for x in df_cluster[df_cluster["Signature_ID"] == cluster_id]["Mutation_frequency"]:
		for y in x.split(","):
			res, _ = y.split(":")
			gene_name, pos = res.split("_")
			if gene_name == gene:
				driver_res = driver_res.union({int(pos)})
	
	# get mutations in the gene
	df_mut = pd.read_csv(maf_file, sep = "\t")
	df_mut = df_mut[df_mut["Variant_Classification"].isin(mis)]
	df_tmp = df_mut[df_mut["Hugo_Symbol"] == gene]
	df_mis = defaultdict(int)
	for x in df_tmp["Protein_position"]:
		for y in all_function_py3.split_consecutive_pos(x):
			df_mis[y] = df_mis[y] + 1
	mis_height = max(df_mis.values()) * unitmut_height
	
	# plot mutation distribution
	total_height = mis_height + block_height + 0.5
	total_len = len(protein_seq) * unitaa_width
	fig, ax = plt.subplots(figsize = (total_len * 1.1, total_height * 1.2))
	if type(gene_color) == str:
		ax.add_patch(Rectangle((0, 0), total_len, block_height, color = gene_color))
	else:
		for key in gene_color:
			ax.add_patch(Rectangle(((key - 1) * unitaa_width, 0), unitaa_width, block_height, color = gene_color[key]))
	plt.text(0.5 * unitaa_width, - 0.4 * block_height, "1", fontsize = fontsize, ha = "center", va = "top", c = "#252525")
	plt.text(total_len - 0.5 * unitaa_width, - 0.4 * block_height, str(len(protein_seq)), fontsize = fontsize, ha = "center", va = "top", c = "#252525")
	for x in sorted(df_mis)[::-1]:
		if x in driver_res:
			color_ball = "red"
			zorder = 50
		else:
			color_ball = "#d9d9d9"
			zorder = 10
		x_pos = (x - 0.5) * unitaa_width
		plt.plot([x_pos, x_pos], [block_height, block_height + df_mis[x] * unitmut_height], color = "black", linewidth = 1, zorder = zorder)
		plt.scatter(x_pos, block_height + df_mis[x] * unitmut_height, s = res_size, marker = "o", c = color_ball, zorder = zorder, edgecolors = 'black')
		
	plt.plot([-5 * unitaa_width, -5 * unitaa_width], [block_height, block_height + mis_height], color = "#252525", linewidth = 2)
	for x in range(max(df_mis.values()) + 1):
		plt.plot([-15 * unitaa_width, -5 * unitaa_width], [block_height + x * unitmut_height, block_height + x * unitmut_height], color = "#252525", \
			linewidth = 2)
		plt.text(-20 * unitaa_width, block_height + x * unitmut_height, str(x), fontsize = fontsize, c = "#252525", ha = "right", va = "center")
	plt.xlim(-total_len * 0.05, total_len * 1.05)
	plt.ylim(-total_height * 0.05, total_height * 1.05)
	plt.axis('off')
	
	if output != None:
		fig = plt.gcf()
		fig.savefig(output, bbox_inches = "tight")
	plt.show()
	plt.close()
	return


def mutation_distribution_mis_triangle(gene, clusters, protein_seq, maf_file, cluster_file, gene_color = "#a1d99b", output = None, \
	unitaa_width = 0.1, non_ires = None, mutation_side = "above"):
	'''
	gene: e.g. "RYBP"
	clusters: e.g. ["RING1_RYBP_Missense1", "RNF2_RYBP_Missense1"]
	maf_file: /local/storage/yz2296/cancer_hotspot_new/data/version_5/Mutation_data/Pan_cancer/Zhang_Mutations.txt
	cluster_file: /local/storage/yz2296/cancer_hotspot_new/website/website/download/pancan_clusters.txt
	mutation_side: 1. above; 2. below
	non_ires: non-interface residues, e.g. "12-16,19-33,35,39-41"
	'''
	mis = {"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"}

	# get clusters where the gene is involved
	df_cluster = pd.read_csv(cluster_file, sep = "\t")
	driver_res = set()
	for x in df_cluster[df_cluster["Signature_ID"].isin(clusters)]["Mutation_frequency"]:
		for y in x.split(","):
			res, _ = y.split(":")
			gene_name, pos = res.split("_")
			if gene_name == gene:
				driver_res = driver_res.union({int(pos)})
	
	# get mutations in the gene
	df_mut = pd.read_csv(maf_file, sep = "\t")
	df_mut = df_mut[df_mut["Variant_Classification"].isin(mis) & (df_mut["Hugo_Symbol"] == gene)]
	df_mis = defaultdict(int)
	for x in df_mut["Protein_position"]:
		for y in all_function_py3.split_consecutive_pos(x):
			df_mis[y] = df_mis[y] + 1
	
	# plot mutation distribution
	fig, ax = plt.subplots(figsize = ((len(protein_seq) + 20) * unitaa_width, 4))
	ax.add_patch(Rectangle((0, 0), len(protein_seq) * unitaa_width, 1, color = gene_color))
	for key in driver_res:
		ax.add_patch(Rectangle(((key - 1) * unitaa_width, 0), unitaa_width, 1, color = "red"))
	if non_ires != None:
		for key in non_ires.split(","):
			if len(key.split("-")) == 1:
				ax.add_patch(Rectangle(((int(key) - 1) * unitaa_width, 0), unitaa_width, 1, color = "white", alpha = 0.5))
			else:
				x, y = key.split("-")
				ax.add_patch(Rectangle(((int(x) - 1) * unitaa_width, 0), (int(y) - int(x) + 1) * unitaa_width, 1, color = "white", alpha = 0.5))
	for x in sorted(df_mis, reverse = True):
		if x in driver_res:
			color = "red"
		else:
			color = "#d9d9d9"
		x_pos = (x - 0.5) * unitaa_width
		# if x >= 1150:
		if mutation_side == "above":
			for i in range(df_mis[x]):
				plt.scatter(x_pos, 1.55 + 0.87 * i, s = 2000, marker = "v", c = color, edgecolors = "k", linewidth = 2)
		else:
			for i in range(df_mis[x]):
				plt.scatter(x_pos, -0.55 - 0.87 * i, s = 2000, marker = "^", c = color, edgecolors = "k", linewidth = 2)
	plt.xlim(-10 * unitaa_width, (len(protein_seq) + 10) * unitaa_width)
	# plt.plot([1149 * unitaa_width, 1149 * unitaa_width], [0,0.5], color = "r", linewidth = 0.5)
	if mutation_side == "above":
		plt.ylim(0, 4)
	else:
		plt.ylim(-3, 1)
	plt.axis('off')
	
	if output != None:
		fig = plt.gcf()
		fig.savefig(output, bbox_inches = "tight")
	plt.show()
	plt.close()
	return

# def mutation_distribution(gene, protein_length, protein_seq, maf_file, cluster_file, gene_color = "#a1d99b", output = None, block_size = 10, block_width = 1, \
# 	block_height = 1, connector_width = 0.6, unitmut_height = 1, text_height = 0.7, lof_flag = False, res_size = 100):
# 	'''
# 	gene: e.g. "MXD1"
# 	maf_file: /local/storage/yz2296/cancer_hotspot_new/data/version_5/Mutation_data/Pan_cancer/Zhang_Mutations.txt
# 	cluster_file: /local/storage/yz2296/cancer_hotspot_new/website/website/download/pancan_clusters.txt
# 	block_width: the width of 10 AAs
# 	gene_color: the value could be str or dict. If dict, keys are protein positions and values are the corresponding colors (e.g. {1: "blue", 2: "green", ...})
# 	'''
# 	lof = {"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site"}
# 	mis = {"In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"}

# 	# get clusters where the gene is involved
# 	df_cluster = pd.read_csv(cluster_file, sep = "\t")
# 	df = df_cluster[df_cluster["Affected_genes"].apply(lambda x: gene in x.split(","))]
# 	df_cluster_lof = df[df["Type"].apply(lambda x: x[0:len("LoF")] == "LoF")]
# 	if df_cluster_lof.shape[0] > 0:
# 		lof_color = "red"
# 	else:
# 		lof_color = "#969696"
# 	df_cluster_mis = df[df["Type"].apply(lambda x: x[0:len("Missense")] == "Missense")]
# 	driver_res = set()
# 	if df_cluster_mis.shape[0] > 0:
# 		for x in df_cluster_mis["Mutation_frequency"]:
# 			for y in x.split(","):
# 				res, _ = y.split(":")
# 				gene_name, pos = res.split("_")
# 				if gene_name == gene:
# 					driver_res = driver_res.union({int(pos)})
	
# 	# get mutations in the gene
# 	df_mut = pd.read_csv(maf_file, sep = "\t")
# 	df_mut = df_mut[df_mut["Variant_Classification"].isin(mis.union(lof))]
# 	df_tmp = df_mut[df_mut["Hugo_Symbol"] == gene]
# 	all_res = set()
# 	for x in df_tmp["Protein_position"]:
# 		all_res = all_res.union(set(all_function_py3.split_consecutive_pos(x)))
# 	df_mis = defaultdict(int)
# 	for x in df_tmp[df_tmp["Variant_Classification"].isin(mis)]["Protein_position"]:
# 		for y in all_function_py3.split_consecutive_pos(x):
# 			df_mis[y] = df_mis[y] + 1
# 	df_lof = defaultdict(int)
# 	for x in df_tmp[df_tmp["Variant_Classification"].isin(lof)]["Protein_position"]:
# 		for y in all_function_py3.split_consecutive_pos(x):
# 			df_lof[y] = df_lof[y] + 1

# 	# get sequence blocks with mutation(s)
# 	range_list = []
# 	if protein_length % block_size == 0:
# 		block_num = int(protein_length / block_size)
# 	else:
# 		block_num = int(protein_length / block_size) + 1
# 	for i in range(block_num):
# 		start = i * block_size + 1
# 		end = min((i + 1) * block_size + 1, protein_length + 1)
# 		in_res = all_res.intersection(set(range(start, end)))
# 		if len(in_res) > 0:
# 			range_list.append((start, end - 1))
# 	print(range_list)
# 	concat_range = []
# 	for start, end in range_list:
# 		concat_range.append(",".join([str(x) for x in range(start, end + 1)]))
# 	concat_range = ",".join(concat_range)
# 	zipped_range = all_function_py3.zip_res_range(concat_range)
# 	range_list = zipped_range[1:-1].split(",")
# 	print(zipped_range)
	
# 	# plot mutation distribution
# 	total_len = 0
# 	start_flag = 0
# 	if range_list[0].split("-")[0] != "1":
# 		total_len = total_len + connector_width
# 		start_flag = 1
# 	end_flag = 0
# 	if int(range_list[-1].split("-")[1]) < protein_length:
# 		total_len = total_len + connector_width
# 		end_flag = 1
# 	total_len = total_len + connector_width * (len(range_list) - 1)
# 	for x in range_list:
# 		start, end = x.split("-")
# 		total_len = total_len + (int(end) + 1 - int(start)) / block_size * block_width
# 	mis_height = max(df_mis.values()) * unitmut_height
# 	if lof_flag == True:
# 		lof_height = max(df_lof.values()) * unitmut_height
# 	else:
# 		lof_height = 0
# 	total_height = mis_height + block_height + lof_height
# 	fig, ax = plt.subplots(figsize = (total_len * 1.2, total_height * 1.2))
# 	cursor = 0
# 	if start_flag == 1:
# 		plt.plot([0, connector_width * 0.35], [0.5 * block_height + lof_height, 0.5 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 		plt.plot([connector_width * 0.65, connector_width], [0.5 * block_height + lof_height, 0.5 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 		plt.plot([connector_width * 0.25, connector_width * 0.45], [0.2 * block_height + lof_height, 0.8 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 		plt.plot([connector_width * 0.55, connector_width * 0.75], [0.2 * block_height + lof_height, 0.8 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 		cursor = cursor + connector_width
# 	if end_flag == 1:
# 		plt.plot([total_len - connector_width, total_len - connector_width * 0.65], [0.5 * block_height + lof_height, 0.5 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 		plt.plot([total_len - connector_width * 0.35, total_len], [0.5 * block_height + lof_height, 0.5 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 		plt.plot([total_len - connector_width * 0.75, total_len - connector_width * 0.55], [0.2 * block_height + lof_height, 0.8 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 		plt.plot([total_len - connector_width * 0.45, total_len - connector_width * 0.25], [0.2 * block_height + lof_height, 0.8 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 	start, end = range_list[0].split("-")
# 	plt.text(cursor, - text_height * 0.95, start, fontsize = 10, ha = "center", c = "#252525")
# 	if type(gene_color) == str:
# 		ax.add_patch(Rectangle((cursor, lof_height), block_width * (int(end) - int(start) + 1) / block_size, block_height, color = gene_color))
# 	# else:
# 	for x in df_mis:
# 		if x >= int(start) and x <= int(end):
# 			mutcount = df_mis[x]
# 			if x in driver_res:
# 				color_ball = "red"
# 			else:
# 				color_ball = "#969696"
# 			x_pos = cursor + block_width / block_size * (x - int(start) + 0.5)
# 			plt.plot([x_pos, x_pos], [block_height + lof_height, block_height + lof_height + mutcount * unitmut_height], color = "#252525", linewidth = 0.5)
# 			plt.scatter(x_pos, block_height + lof_height + mutcount * unitmut_height, s = res_size, marker = "o", c = color_ball, zorder = 10)
# 	if lof_flag == True:
# 		for x in df_lof:
# 			if x >= int(start) and x <= int(end):
# 				mutcount = df_lof[x]
# 				x_pos = cursor + block_width / block_size * (x - int(start) + 0.5)
# 				plt.plot([x_pos, x_pos], [lof_height - mutcount * unitmut_height, lof_height], color = "#252525", linewidth = 0.5)
# 				plt.scatter(x_pos, lof_height - mutcount * unitmut_height, s = res_size, marker = "x", c = lof_color, zorder = 10, linewidth = 5)
# 	plt.plot([cursor, cursor + block_width * (int(end) - int(start) + 1) / block_size], [- 0.4 * text_height, - 0.4 * text_height], color = "#252525", linewidth = 0.5)
# 	plt.plot([cursor, cursor], [- 0.6 * text_height, - 0.4 * text_height], color = "#252525", linewidth = 0.5)
# 	plt.plot([cursor + block_width * (int(end) - int(start) + 1) / block_size, cursor + block_width * (int(end) - int(start) + 1) / block_size], [- 0.6 * text_height, - 0.4 * text_height], color = "#252525", linewidth = 0.5)
# 	cursor = cursor + block_width * (int(end) - int(start) + 1) / block_size
# 	plt.text(cursor, - text_height * 0.95, end, fontsize = 10, ha = "center", c = "#252525")
# 	for x in range_list[1:]:
# 		plt.plot([cursor, cursor + connector_width * 0.35], [0.5 * block_height + lof_height, 0.5 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 		plt.plot([cursor + connector_width * 0.65, cursor + connector_width], [0.5 * block_height + lof_height, 0.5 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 		plt.plot([cursor + connector_width * 0.25, cursor + connector_width * 0.45], [0.2 * block_height + lof_height, 0.8 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 		plt.plot([cursor + connector_width * 0.55, cursor + connector_width * 0.75], [0.2 * block_height + lof_height, 0.8 * block_height + lof_height], color = "#252525", linewidth = 0.5)
# 		cursor = cursor + connector_width
# 		start, end = x.split("-")
# 		plt.text(cursor, - text_height * 0.95, start, fontsize = 10, ha = "center", c = "#252525")
# 		ax.add_patch(Rectangle((cursor, lof_height), block_width * (int(end) - int(start) + 1) / block_size, block_height, color = gene_color))
# 		for x in df_mis:
# 			if x >= int(start) and x <= int(end):
# 				mutcount = df_mis[x]
# 				x_pos = cursor + block_width / block_size * (x - int(start) + 0.5)
# 				if x in driver_res:
# 					color_ball = "red"
# 					plt.text(x_pos, block_height + lof_height + mutcount * unitmut_height + 0.3 * unitmut_height, protein_seq[x - 1]  + str(x), fontsize = 24, color = "#252525", rotation = 90, va = "bottom", ha = "center", zorder = 20)
# 				else:
# 					color_ball = "#969696"
# 				plt.plot([x_pos, x_pos], [block_height + lof_height, block_height + lof_height + mutcount * unitmut_height], color = "#252525", linewidth = 0.5)
# 				plt.scatter(x_pos, block_height + lof_height + mutcount * unitmut_height, s = res_size, marker = "o", c = color_ball, zorder = 10)
# 		if lof_flag == True:
# 			for x in df_lof:
# 				if x >= int(start) and x <= int(end):
# 					mutcount = df_lof[x]
# 					x_pos = cursor + block_width / block_size * (x - int(start) + 0.5)
# 					plt.plot([x_pos, x_pos], [lof_height - mutcount * unitmut_height, lof_height], color = "#252525", linewidth = 0.5)
# 					plt.scatter(x_pos, lof_height - mutcount * unitmut_height, s = res_size, marker = "x", c = lof_color, zorder = 10, linewidth = 5)
# 		plt.plot([cursor, cursor + block_width * (int(end) - int(start) + 1) / block_size], [- 0.4 * text_height, - 0.4 * text_height], color = "#252525", linewidth = 0.5)
# 		plt.plot([cursor, cursor], [- 0.6 * text_height, - 0.4 * text_height], color = "#252525", linewidth = 0.5)
# 		plt.plot([cursor + block_width * (int(end) - int(start) + 1) / block_size, cursor + block_width * (int(end) - int(start) + 1) / block_size], [- 0.6 * text_height, - 0.4 * text_height], color = "#252525", linewidth = 0.5)
# 		cursor = cursor + block_width * (int(end) - int(start) + 1) / block_size
# 		plt.text(cursor, - text_height * 0.95, end, fontsize = 10, ha = "center", c = "#252525")
		
# 	plt.plot([- 0.4 * text_height, - 0.4 * text_height], [lof_height + block_height, lof_height + block_height + mis_height], color = "#252525", linewidth = 2)
# 	plt.plot([- 0.4 * text_height, - 0.4 * text_height], [0, lof_height], color = "#252525", linewidth = 2)
# 	for x in range(mis_height + 1):
# 		plt.plot([- 0.6 * text_height, - 0.4 * text_height], [lof_height + block_height + x, lof_height + block_height + x], color = "#252525", linewidth = 2)
# 		plt.text(- 0.8 * text_height, lof_height + block_height + x, str(x), fontsize = 24, c = "#252525", ha = "right", va = "center")
# 	if lof_flag == True:
# 		for x in range(lof_height + 1):
# 			plt.plot([- 0.6 * text_height, - 0.4 * text_height], [x, x], color = "#252525", linewidth = 2)
# 			plt.text(- 0.8 * text_height, x, str(lof_height - x), fontsize = 24, c = "#252525", ha = "right", va = "center")
# 	plt.xlim(- total_len * 0.15, total_len * 1.05)
# 	plt.ylim(- total_height * 0.15, total_height * 1.05)
# 	plt.axis('off')
	
# 	if output != None:
# 		fig = plt.gcf()
# 		fig.savefig(output, bbox_inches = "tight")
# 	plt.show()
# 	plt.close()
# 	return


def DiGraph_layout_FineGrained(G, node2pos, node2color, output = None, specific_edge = {}, node_radius = 0.3, node_rim = 5, edge_radius = 0.25, edge_color = "r", \
	edge_magnifier = 50, arrow_size = 6):
	'''
	[NOTE] This function does not display self loops
	G: directed graph
	node2pos: e.g. {"MAX": (0, 0), "MLX": (2, 0), "MXD1": (1, 1), "MXD4": (1, 0), "MNT": (1,-1)}
	node2color: e.g. {"MAX": "blue", "MLX": "blue", "MXD1": "orange", "MXD4": "orange", "MNT": "orange"}
	specific_edge: in the format of {(node1, node2): [start_pos, end_pos, arc_radius]}, e.g. {("MLX", "MAX"): [(2, 0.3), (0, 0.3), 1.2]}
	edge_magnifier: the actual edge weight might be too small for displaying purpose. This value specifies to what extent you want to amplify the edge weight
	'''
	min_x = min([x - node_radius for x, y in node2pos.values()])
	max_x = max([x + node_radius for x, y in node2pos.values()])
	min_y = min([y - node_radius for x, y in node2pos.values()])
	max_y = max([y + node_radius for x, y in node2pos.values()])
	plt.figure(figsize = ((max_x - min_x) * 1.1 * 1.5 ,(max_y - min_y) * 1.1 * 1.5))
	plt.axes().set_aspect(1)
	plt.xlim(min_x - 0.05 * (max_x - min_x), max_x + 0.05 * (max_x - min_x))
	plt.ylim(min_y - 0.05 * (max_y - min_y), max_y + 0.05 * (max_y - min_y))
	plt.axis('off')
	ax = plt.gca()
	for key in node2pos:
		x, y = node2pos[key]
		circle = plt.Circle((x, y), node_radius, edgecolor = node2color[key], linewidth = node_rim, fill = False)
		ax.add_artist(circle)
		ax.text(x, y, key, va = "center", ha = "center")

	kw = dict(arrowstyle = "Simple, tail_width=0, head_width=" + str(arrow_size) + ", head_length=" + str(arrow_size), color = edge_color)
	for u,v,info in G.edges(data = True):
		if u != v:
			if (u,v) not in specific_edge:
				x1, y1 = node2pos[u]
				x2, y2 = node2pos[v]
				if x1 == x2:
					if y2 > y1:
						theta = 0.25 * 2 * np.pi
					else:
						theta = -0.25 * 2 * np.pi
				else:
					theta = math.atan((y2 - y1) / (x2 - x1))
					if x2 - x1 < 0:
						theta = theta + np.pi
				x1 = x1 + node_radius * 1.05 * math.cos(theta - 0.2)
				y1 = y1 + node_radius * 1.05 * math.sin(theta - 0.2)
				x2 = x2 + node_radius * 1.05 * math.cos(theta + np.pi + 0.2)
				y2 = y2 + node_radius * 1.05 * math.sin(theta + np.pi + 0.2)
				rad = edge_radius
			else:
				pos1, pos2, rad = specific_edge[(u,v)]
				x1, y1 = pos1
				x2, y2 = pos2
			a = patches.FancyArrowPatch((x1, y1), (x2, y2), connectionstyle = "arc3,rad=" + str(rad), **kw, lw = info["weight"] * edge_magnifier, zorder = 50)
			ax.add_patch(a)

	if output != None:
		plt.gcf().savefig(output, bbox_inches = "tight")
	plt.show()
	plt.close()
	return


def Graph_layout_CoarseGrained(G, node2pos, cluster_marker = None, marker_zoom = 0.05, node2radius = None, edge_magnifier = 1, \
	output = None, with_labels = False, node2color = None, edge2color = None, noderim_lw = 5, noderim_color = None):
	'''
	G: undirected graph. All edges have attribute **weight**. Some nodes and edges have attribute **type** whose value is a subset of {"missense", "LoF"}. 
	'''
	if node2radius == None:
		node2radius = {x: 0.35 for x in node2pos}
	if edge2color == None:
		edge2color = {tuple(sorted([u,v])): "#d9d9d9" for u,v in G.edges()}
	if node2color == None:
		node2color = {u: "#d9d9d9" for u in G.nodes()}
	x_list = []
	y_list = []
	for key in node2pos:
		x, y = node2pos[key]
		x_list.extend([x - node2radius[key], x + node2radius[key]])
		y_list.extend([y - node2radius[key], y + node2radius[key]])
	plt.figure(figsize = ((max(x_list) - min(x_list)) * 1.05, (max(y_list) - min(y_list)) * 1.05))
	plt.axes().set_aspect(1)
	plt.xlim(min(x_list) - 0.025 * (max(x_list) - min(x_list)), max(x_list) + 0.025 * (max(x_list) - min(x_list)))
	plt.ylim(min(y_list) - 0.025 * (max(y_list) - min(y_list)), max(y_list) + 0.025 * (max(y_list) - min(y_list)))
	plt.axis('off')
	ax = plt.gca()
	for key in node2pos:
		x, y = node2pos[key]
		# if "type" in G.nodes[key]:
		# 	if "LoF" in G.nodes[key]["type"]:
		# 		edgecolor = "red"
		# 	else:
		# 		edgecolor = boundary_color
		# 		ax = all_function_py3.display_image(cluster_marker, ax, x, y, zoom = marker_zoom)
		# else:
		# 	edgecolor = boundary_color
		# circle = plt.Circle((x, y), node2radius[key], edgecolor = edgecolor, linewidth = noderim_lw, facecolor = node2color[key])
		circle = plt.Circle((x, y), node2radius[key], facecolor = node2color[key], edgecolor = noderim_color, linewidth = noderim_lw, zorder = 50)
		ax.add_artist(circle)
		if "alpha" in G.nodes[key]:
			circle = plt.Circle((x, y), node2radius[key], facecolor = "white", zorder = 50, alpha = 1 - G.nodes[key]["alpha"])
			ax.add_artist(circle)
		if with_labels == True:
			ax.text(x, y, key, fontsize = 20, color = "k", zorder = 50, ha = "center", va = "center")
	for u,v,info in G.edges(data = True):
		x1, y1 = node2pos[u]
		x2, y2 = node2pos[v]
		lw = info["weight"] * edge_magnifier
		ax.plot([x1, x2], [y1, y2], lw = lw, color = edge2color[tuple(sorted([u,v]))], zorder = 1)#"#525252"
		if "type" in info and "missense" in info["type"]:
			ax = all_function_py3.display_image(cluster_marker, ax, (x1 + x2) / 2, (y1 + y2) / 2, zoom = marker_zoom)
	if output != None:
		fig = plt.gcf()
		fig.savefig(output, bbox_inches = "tight", transparent = True)
	return ax


def Graph_layout_FineGrained(G, node2info, node2format, edge2info = {}, node_rim = 5, edge_radius = 0.25, edge_color = "#252525", edge_magnifier = 1, \
	arrow_size = 6, fontsize = 50, zoom = 1):
	'''
	[NOTE] This function does not display self loops
	G: undirected graph
	node2info: e.g. {"MAX": ((0, 0), "/local/storage/yz2296/cancer_hotspot_new/website/figures/figure_4/MAX_cartoon.png"), ...} The position and the image of each node.
	node2format: e.g. {"MAX": ("blue", 2, ("orange", "bold")), "MLX": ("black", 1, ("black", "normal")), ...}. The first element in each tuple stands for node color. The second element stands for node radius. The third element stands for the text color and font weight.
	edge2info: e.g. {("MAX", "MNT"):  (1, None), ("MAX", "MAX"): (2, "left"), ...}. The first element stands for the radius of an edge, and the second element is the position of self-loop (left/right/top/bottom).
	edge_magnifier: the actual edge weight might be too small for displaying purpose. This value specifies to what extent you want to amplify the edge weight
	'''
	x_list = []
	y_list = []
	for key in node2info:
		(x, y), _ = node2info[key]
		_, node_radius, _ = node2format[key]
		if G.has_edge(key, key):
			rad, pos = edge2info[(key, key)]
			if pos == "left":
				x_list.extend([x - 1.5 * node_radius - rad, x + node_radius])
				x_list.extend([y - node_radius, y + node_radius])
			elif pos == "right":
				x_list.extend([x + 1.5 * node_radius + rad, x - node_radius])
				y_list.extend([y - node_radius, y + node_radius])
			elif pos == "top":
				x_list.extend([x - node_radius, x + node_radius])
				y_list.extend([y + 1.5 * node_radius + rad, y - node_radius])
			else:
				x_list.extend([x - node_radius, x + node_radius])
				y_list.extend([y - 1.5 * node_radius - rad, y + node_radius])
		else:
			x_list.extend([x - node_radius, x + node_radius])
			y_list.extend([y - node_radius, y + node_radius])
	min_x = min(x_list)
	max_x = max(x_list)
	min_y = min(y_list)
	max_y = max(y_list)
	plt.figure(figsize = ((max_x - min_x) * 1.05 * zoom, (max_y - min_y) * 1.1 * zoom))
	plt.axes().set_aspect(1)
	plt.xlim(min_x - 0.025 * (max_x - min_x), max_x + 0.025 * (max_x - min_x))
	plt.ylim(min_y - 0.025 * (max_y - min_y), max_y + 0.025 * (max_y - min_y))
	plt.axis('off')
	ax = plt.gca()
	for key in node2info:
		(x, y), image = node2info[key]
		node_color, node_radius, (text_color, fontweight) = node2format[key]
		circle = plt.Circle((x, y), node_radius, edgecolor = node_color, linewidth = node_rim * zoom, facecolor = "white")
		ax.add_artist(circle)
		ax.text(x, y + node_radius * 1.06, key, va = "bottom", ha = "center", fontsize = fontsize, color = text_color, fontweight = fontweight)
		ax = all_function_py3.display_image(image, ax, x, y, zoom = 1)

	# kw = dict(arrowstyle = "<->, head_width=" + str(arrow_size) + ", head_length=" + str(arrow_size), color = edge_color)
	kw = dict(arrowstyle = "Simple, tail_width=0, head_width=0, head_length=0.1", color = edge_color)
	for u,v,info in G.edges(data = True):
		if (v, u) in edge2info:
			u, v = v, u
		(x1, y1), _ = node2info[u]
		(x2, y2), _ = node2info[v]
		# _, n1_radius, _ = node2format[u]
		# _, n2_radius, _ = node2format[v]
		if u != v:
			if (u,v) not in edge2info:
				if x1 == x2 or y1 == y2:
					rad = 0
				else:
					rad = edge_radius
			else:
				rad, _ = edge2info[(u,v)]
			# if x1 == x2:
			# 	if y2 > y1:
			# 		theta = 0.25 * 2 * np.pi
			# 	else:
			# 		theta = -0.25 * 2 * np.pi
			# else:
			# 	theta = math.atan((y2 - y1) / (x2 - x1))
			# 	if x2 - x1 < 0:
			# 		theta = theta + np.pi
			# theta = theta + np.pi / 6
			# x1 = x1 + n1_radius * 1.05 * math.cos(theta - 0.2)
			# y1 = y1 + n1_radius * 1.05 * math.sin(theta - 0.2)
			# x2 = x2 + n2_radius * 1.05 * math.cos(theta + np.pi + 0.2)
			# y2 = y2 + n2_radius * 1.05 * math.sin(theta + np.pi + 0.2)
			a = patches.FancyArrowPatch((x1, y1), (x2, y2), connectionstyle = "arc3,rad=" + str(rad), **kw, \
				lw = info["weight"] * edge_magnifier * zoom, zorder = 0)
			ax.add_patch(a)
		else:
			rad, pos = edge2info[(u, v)]
			_, node_radius, _ = node2format[u]
			if pos == "left":
				x = x1 - node_radius * 1.5
				y = y1
			elif pos == "right":
				x = x1 + node_radius * 1.5
				y = y1
			elif pos == "top":
				x = x1
				y = y1 + node_radius * 1.5
			else:
				x = x1
				y = y1 - node_radius * 1.5
			circle = plt.Circle((x, y), rad, edgecolor = edge_color, linewidth = info["weight"] * edge_magnifier * zoom, zorder = 0, \
				fill = False)
			ax.add_artist(circle)
	return ax


def module_layout(gene2color, gene2uniprot, cluster_file, mutation_file, output = None, fontsize = 32, cancer2color = {"coad": "#fb8072", "ucec": "#b3de69", "luad": "#bc80bd", "brca": "#bf812d", "hnsc": "#1f78b4", "read": "#ff7f00", "blca": "#33a02c", "kirp": "#8c510a", "stad": "#f1b6da", "lusc": "#c51b7d", "paad": "#35978f"}):
	'''
	gene2color: e.g. {"MXD1": "#70d3ff", "MAX": "#fbba1d", "MXD4": "#70d3ff", "MLX": "#fbba1d", "MNT": "#70d3ff"}
	cluster_file: e.g. /local/storage/yz2296/cancer_hotspot_new/website/website/download/pancan_clusters.txt
	mutation_file: e.g. /local/storage/yz2296/cancer_hotspot_new/data/version_5/Mutation_data/Pan_cancer/Zhang_Mutations.txt
	'''
	def get_cluster_info(row):
		for x in row["Mutation_frequency"].split(","):
			res, mutcount = x.split(":")
			cluster2res["-".join(row["Signature_ID"].split("_")[0:-1]) + "(" + row["Structure_source"] + ")"].append(res)
			res2mutcount[res] = int(mutcount)
			gene, pos = res.split("_")
			gene2res[gene] = gene2res[gene].union({int(pos)})
		return True
	def plot_mut(ax, row, tumor2x, res2y, radius = 1):
		for pos in all_function_py3.split_consecutive_pos(row["Protein_position"]):
			x = tumor2x[row["Tumor_Sample_Barcode"]]
			y = res2y[row["Hugo_Symbol"] + "_" + str(pos)]
			circle = plt.Circle((x, y), radius, color = "r", zorder = 30)
			ax.add_artist(circle)
			ax.text(x, y, row["Amino_acids"].split("/")[1], color = "white", zorder = 50, fontsize = fontsize, va = "center", \
					ha = "center")
		return True

	# prepare supporting information
	margin_width = 0
	margin_height = 0.5
	legend_width = 4
	uniprot2seq = all_function_py3.get_proseq([gene2uniprot[x] for x in gene2color])
	
	# load cluster information
	df_cluster = pd.read_csv(cluster_file, sep = "\t")
	module_clusters = df_cluster[df_cluster["Affected_genes"].apply(lambda x: set(x.split(",")) - set(gene2color) == set())]
	module_clusters = module_clusters[module_clusters["Type"] != "LoF_IntraProtein_Enriched"]
	
	# load mutation information
	df_driver = pd.read_csv(mutation_file, sep = "\t", dtype = str)
	df_driver = df_driver[df_driver["Variant_Classification"].isin({"Missense_Mutation", "In_Frame_Del", "In_Frame_Ins"})]
	module_drivers = df_driver[df_driver["Hugo_Symbol"].isin(gene2color)]

	# figure preparation
	res2mutcount = {}
	gene2res = defaultdict(set)
	cluster2res = defaultdict(list)
	module_clusters["status"] = module_clusters.apply(get_cluster_info, axis = 1)
	# gene2res = dict(gene2res)
	module_drivers = module_drivers[module_drivers.apply(lambda x: set(all_function_py3.split_consecutive_pos(x["Protein_position"])).intersection(gene2res[x["Hugo_Symbol"]]) != set(), axis = 1)]
	cancer_list = []
	for cancer,df_cancer in module_drivers.groupby("TUMORTYPE"):
		cancer_list.append((cancer, sorted(df_cancer["Tumor_Sample_Barcode"].unique().tolist())))
	cancer_list = sorted(cancer_list, key = lambda x: len(x[1]), reverse = True)

	# draw figure
	box_width = len(module_drivers["Tumor_Sample_Barcode"].unique()) + 2 * margin_width * len(cancer_list)
	box_height = len(res2mutcount) + 2 * margin_height * len(gene2color)
	x_min = -7 - len(cluster2res)
	x_max = box_width + 1 + max(res2mutcount.values()) + 0.5
	y_min = -box_height - 0.5
	y_max = 7
	fig, ax = plt.subplots(figsize = (x_max - x_min, y_max - y_min))
	plt.xlim(x_min, x_max)
	plt.ylim(y_min, y_max)
	cursor = 0
	res2y = {}
	for x in gene2color:
		ax.fill([0, box_width, box_width, 0], [-cursor, -cursor, -cursor - len(gene2res[x]) - 2 * margin_height, \
			-cursor - len(gene2res[x]) - 2 * margin_height], gene2color[x], alpha = 0.5)
		ax.text(-5 - len(cluster2res), -cursor - 0.5 * len(gene2res[x]) - margin_height, x, fontsize = int(1.5 * fontsize), color = "#252525", \
			ha = "center", va = "center")
		for i,y in enumerate(sorted(gene2res[x])):
			ax.text(-0.5, -cursor - margin_height - i - 0.5, uniprot2seq[gene2uniprot[x]][y - 1] + str(y), fontsize = fontsize, ha = "right", \
				va = "center", color = "#252525")
			res2y[x + "_" + str(y)] = -cursor - margin_height - i - 0.5
			ax.plot([-0.3, 0], [-cursor - margin_height - i - 0.5, -cursor - margin_height - i - 0.5], linewidth = 2, color = "#252525")
		cursor = cursor + len(gene2res[x]) + 2 * margin_height
	cursor = 0
	tumor2x = {}
	for cancer,tumors in cancer_list:
		for i,x in enumerate(tumors):
			ax.fill([cursor + i + 0.1 + margin_width, cursor + i + 0.9 + margin_width, cursor + i + 0.9 + margin_width, \
				cursor + i + 0.1 + margin_width], [1, 1, 2, 2], cancer2color[cancer])
			tumor2x[x] = cursor + i + 0.5 + margin_width
		cursor = cursor + len(tumors) + 2 * margin_width
	module_drivers["status"] = module_drivers.apply(lambda x: plot_mut(ax, x, tumor2x, res2y, radius = 0.4), axis = 1)
	for i,x in enumerate(cluster2res):
		y_list = []
		for y in cluster2res[x]:
			circle = plt.Circle((-3 - i, res2y[y]), 0.15, color = "#252525")
			ax.add_artist(circle)
			y_list.append(res2y[y])
		ax.text(-3 - i, 0, x, rotation = 90, ha = "center", va = "bottom", fontsize = fontsize)
		ax.plot([-3 - i, -3 - i], [min(y_list), max(y_list)], color = "#252525", linewidth = 2)
	for x in res2mutcount:
		ax.fill([box_width + 1, box_width + 1 + res2mutcount[x], box_width + 1 + res2mutcount[x], box_width + 1], [res2y[x] + 0.4, res2y[x] + 0.4, \
			res2y[x] - 0.4, res2y[x] - 0.4], "#252525")
	ax.text(box_width + 1 + 0.5 * max(res2mutcount.values()), 1.5, "Mutation frequency", ha = "center", va = "center", fontsize = fontsize, \
		color = "#252525")
	ax.plot([box_width + 1, box_width + 1 + max(res2mutcount.values())], [-0.2, -0.2], linewidth = 5, color = "#252525")
	for i in range(max(res2mutcount.values()) + 1):
		ax.plot([box_width + 1 + i, box_width + 1 + i], [-0.2, 0.1], linewidth = 5, color = "#252525")
		ax.text(box_width + 1 + i, 0.2, str(i), fontsize = fontsize, va = "bottom", ha = "center", color = "#252525")
	one_row = 6 #int(len(tumor2x) / legend_width)
	for i,x in enumerate(cancer2color):
		i = i + 1
		if i % one_row == 0:
			v_pos = (i / one_row - 1) * 2
			h_pos = (one_row - 1) * legend_width
		else:
			v_pos = int(i / one_row) * 2
			h_pos = (i % one_row - 1) * legend_width
		ax.fill([h_pos + 0.1, h_pos + 0.9, h_pos + 0.9, h_pos + 0.1], [5 - v_pos, 5 - v_pos, 6 - v_pos, 6 - v_pos], \
				cancer2color[x])
		ax.text(h_pos + 1.5, 5.5 - v_pos, x, va = "center", ha = "left", fontsize = fontsize, color = "#252525")
	plt.axis("off")
	fig = plt.gcf()
	if output != None:
		fig.savefig(output, bbox_inches = "tight")
	return [fig, ax]


############################################################################################################################################
# figure 3
############################################################################################################################################
# module2OrderedGenes = {\
# #                        "ALB,APC,B2M,BTRC,CD1A,CD1B,CD1C,CD1D,CD8A,CDH5,COX6A1,COX6B1,COX6C,CTNNB1,CTNNBIP1,FBXW11,FCGRT,HFE,HLA-A,HLA-B,HLA-C,HLA-DRA,HLA-DRB1,HLA-DRB5,HLA-E,HLA-F,HLA-G,HMHA1,IDH2,IRS2,JUP,KDR,LILRB2,LSP1,MR1,MT-CO1,MT-CO2,MT-CO3,NR5A2,PDIA3,PPP1CB,PPP1R16B,PPP1R3B,PYGB,PYGL,SYNM,TAPBP,TAPBPL,VEGFC": [("TAPBPL","MR1","HLA-F","HLA-G","LILRB2","HLA-E","HLA-C","B2M","HFE","CD1A","CD1B","CD1C","CD1D","CD8A","HLA-A","HMHA1","HLA-DRA","HLA-DRB5","HLA-DRB1","TAPBP","HLA-B"), ("COX6B1","MT-CO3","COX6A1","MT-CO1","COX6C","MT-CO2"), ("CTNNB1","APC","JUP","CDH5","CTNNBIP1","FBXW11","BTRC")], \
# #                        "ASF1A,CRCP,ERCC2,ERCC3,GTF2A1,GTF2A2,GTF2B,GTF2E1,GTF2E2,GTF2F1,GTF2F2,GTF2H1,GTF2H3,GTF2H4,HIRA,POLR1D,POLR2A,POLR2B,POLR2F,POLR2H,POLR2I,POLR2K,POLR2L,POLR3A,POLR3B,POLR3C,POLR3D,POLR3E,POLR3F,POLR3H,POLR3K,TAF1,TAF10,TAF11,TAF2,TAF3,TAF4,TAF5,TAF6,TAF7,TAF8,TAF9,TBP": [("GTF2F2","GTF2F1","POLR2I","GTF2B","POLR2B","POLR2K","POLR2L","POLR3B","POLR3E","POLR3D","POLR1D","POLR3F","POLR3C","POLR3A","POLR3H","CRCP","POLR3K","POLR2F","POLR2H","POLR2A","TAF10","TAF3","TAF5","TAF6","TAF8","TAF2","TAF1","TAF7","ASF1A","HIRA","TBP","TAF9","TAF11","TAF4","GTF2A1","GTF2A2","GTF2E2","GTF2E1","ERCC2","ERCC3","GTF2H1","GTF2H3","GTF2H4")], \
# #                        "ICT1,MRP63,MRPL10,MRPL11,MRPL13,MRPL14,MRPL15,MRPL16,MRPL18,MRPL19,MRPL20,MRPL23,MRPL24,MRPL27,MRPL3,MRPL30,MRPL32,MRPL35,MRPL38,MRPL39,MRPL4,MRPL40,MRPL41,MRPL42,MRPL46,MRPL47,MRPL48,MRPL51,MRPL52,MRPL53,MRPL54,MRPS18A,MRPS30": [("MRP63", "MRPS18A", "MRPL13", "MRPL53", "MRPL10", "MRPL16", "MRPL11", "MRPL54", "ICT1", "MRPL15", "MRPL35", "MRPL14", "MRPL19", "MRPL3", "MRPL39", "MRPL32", "MRPL46", "MRPL40", "MRPL48", "MRPL27", "MRPL18", "MRPL38", "MRPL42", "MRPL20", "MRPL52", "MRPL30", "MRPL4", "MRPL51", "MRPL41", "MRPL24", "MRPL47", "MRPL23", "MRPS30")]
# #                        "ARID1A,ARNT,ARNTL,BRD4,CBFB,COPS2,COPS3,COPS4,COPS7B,DPF2,EGLN1,EPAS1,GPS1,HIF1A,NPAS2,RUNX1,SMARCA4,SMARCC2,SMARCE1,SOCS2,TCEB1,TCEB2,VHL,WHSC1L1": [("WHSC1L1", "BRD4", "EPAS1", "ARNT", "ARNTL", "NPAS2", "EGLN1", "HIF1A", "VHL", "SMARCA4", "SMARCE1", "SMARCC2", "DPF2", "ARID1A", "COPS3", "COPS7B", "COPS4", "COPS2", "GPS1", "CBFB", "RUNX1", "TCEB1", "TCEB2", "SOCS2")], \
# #                        "CREBBP,FOSB,FOSL2,IRF3,JUN,JUNB,JUND,KANSL1,KMT2A,MEN1,MTDH,MYB,NCOA1,NCOA3,NCOR1,NR3C1,PABPC1,PAIP1,PPARG,PSIP1,RXRB,SND1,WDR5": [("NCOA1", "CREBBP", "IRF3", "JUN", "FOSL2", "JUNB", "FOSB", "JUND", "MEN1", "PSIP1", "KMT2A", "WDR5", "KANSL1", "MTDH", "SND1", "MYB", "PAIP1", "PABPC1", "NCOA3", "NR3C1", "NCOR1", "PPARG", "RXRB")], \
# #                        "BCL6,BCOR,BCORL1,DOT1L,KDM2B,MLLT3,PCGF1,PCGF3,PCGF5,PCGF6,RING1,RNF2,RYBP": [("PCGF5", "BCL6", "BCOR", "MLLT3", "DOT1L", "KDM2B", "BCORL1", "PCGF1", "RING1", "RYBP", "RNF2", "PCGF6", "PCGF3")], \
#                        "CASP1,CASP4,ELF4,GSDMD,IL18,IL18R1,PYCARD": [("CASP1", "GSDMD", "CASP4", "ELF4", "IL18R1", "IL18", "PYCARD")], \
# #                        "PSMA2,PSMA4,PSMA5,PSMA6,PSMA7,PSMB10,PSMB7,PSMC1,PSMC2,PSMC3,PSMC4,PSMC6,PSMD1,PSMD11,PSMD12,PSMD13,PSMD14,PSMD3,PSMD6,PSMD7,PSMD8,PSME4,SHFM1": [("PSMA2", "PSMA6", "PSMB10", "PSMB7", "PSMD12", "PSMD13", "PSMD7", "PSMA5", "PSME4", "PSMA7", "PSMA4", "PSMC1", "PSMC2", "PSMC3", "PSMC6", "PSMC4", "PSMD8", "PSMD14", "PSMD1", "PSMD3", "SHFM1", "PSMD6", "PSMD11")], \
# #                        "CAND1,CUL1,CUL4A,CUL5,DCUN1D1,GLMN,NAE1,NEDD8,RBX1,RNF7,UBA3,UBE2M": [("RNF7", "CUL4A", "RBX1", "CUL5", "GLMN", "NEDD8", "UBA3", "NAE1", "UBE2M", "DCUN1D1", "CUL1", "CAND1")], \
# #                        "AQR,BUD31,CDC40,CDC5L,CRNKL1,PPIE,PPIL1,RBM22,SNRNP40,SYF2,XAB2": [("CRNKL1", "CDC5L", "SYF2", "PPIL1", "CDC40", "BUD31", "SNRNP40", "RBM22", "AQR", "PPIE", "XAB2")], \
# #                        "NDUFA10,NDUFA5,NDUFA7,NDUFS1,NDUFS2,NDUFS3,NDUFS4,NDUFS7,NDUFS8,NDUFV1,NDUFV3": [("NDUFA5", "NDUFS2", "NDUFA10", "NDUFS3", "NDUFV3", "NDUFS4", "NDUFV1", "NDUFS1", "NDUFA7", "NDUFS7", "NDUFS8")], \
# #                        "FKBP3,FLCN,FNIP2,MAPKAP1,MTOR,RICTOR,RPTOR,RRAGA,RRAGC,SLC38A9": [("FNIP2", "RRAGC", "SLC38A9", "FKBP3", "MTOR", "RICTOR", "MAPKAP1", "RPTOR", "RRAGA", "FLCN")], \
# #                        "CCT2,CCT3,CCT4,CCT5,CCT6A,CCT7,CCT8,PFDN2,PFDN4,TCP1,VBP1": [("CCT3", "TCP1", "CCT4", "CCT2", "CCT5", "CCT7", "CCT8", "CCT6A", "PFDN2", "VBP1", "PFDN4")], \
# #                        "PPP2CA,PPP2R1A,PPP2R1B,PPP2R2A,PPP2R5A,PPP2R5C,PPP2R5D,STRN3": [("PPP2R5D", "PPP2R1B", "PPP2R2A", "PPP2R5A", "STRN3", "PPP2CA", "PPP2R5C", "PPP2R1A")], \
# #                        "MZT1,REPS1,TUBGCP2,TUBGCP3,TUBGCP6,WWC2": [("TUBGCP6", "MZT1", "REPS1", "WWC2", "TUBGCP2", "TUBGCP3")], \
# #                        "CDC73,CTR9,LEO1,PAF1,RTF1,WDR61": [("WDR61", "CTR9", "PAF1", "LEO1", "RTF1", "CDC73")], \
# #                        "EIF2B1,EIF2B2,EIF2B4,EIF2B5,EIF2S1,EIF2S3": [("EIF2B1", "EIF2B4", "EIF2B2", "EIF2B5", "EIF2S3", "EIF2S1")], \
# #                        "MMP14,MMP2,MMP9,THBS1,TIMP1,TIMP2": [("TIMP2", "MMP2", "THBS1", "MMP9", "TIMP1", "MMP14")], \
# #                        "MCM2,MCM3,MCM4,MCM5,MCM6,MCM7": [("MCM3", "MCM5", "MCM2", "MCM6", "MCM4", "MCM7")], \
# #                        "CYC1,MT-CYB,UQCR10,UQCRC1,UQCRC2,UQCRFS1": [("CYC1", "MT-CYB", "UQCR10", "UQCRC1", "UQCRC2", "UQCRFS1")], \
# #                        "ARHGDIA,ARHGDIB,PAK4,RAC2,RAC3,RHOH": [("ARHGDIA", "RAC3", "PAK4", "RHOH", "ARHGDIB", "RAC2")], \
# #                        "ADM,ADM2,CALCRL,RAMP1,RAMP2,RAMP3": [("CALCRL", "RAMP1", "RAMP2", "ADM", "RAMP3", "ADM2")], \
# #                        "CHRAC1,DRAP1,NFYA,NFYB,NFYC,POLE3": [("NFYA", "NFYC", "CHRAC1", "POLE3", "DRAP1", "NFYB")], \
# #                        "CUL3,FBXL17,H2AFY,KEAP1,NFE2L2,SPOP": [("KEAP1", "FBXL17", "NFE2L2", "H2AFY", "SPOP", "CUL3")], \
# #                        "ZKSCAN4,ZNF213,ZNF446,ZNF496,ZSCAN12,ZSCAN21": [("ZKSCAN4", "ZNF496", "ZNF446", "ZSCAN12", "ZSCAN21", "ZNF213")], \
# #                        "CAB39,MST4,PDCD10,STK24,STK25,STRADA": [("STK24", "STRADA", "CAB39", "MST4", "STK25", "PDCD10")], \
# #                        "RIOK1,RIOK2,RPS14,RPS16,RPS18,RPS25,RPS28,RPS5": [("RIOK2", "RPS18", "RPS25", "RPS16", "RIOK1", "RPS14", "RPS28", "RPS5")], \
# #                        "MED11,MED14,MED17,MED18,MED22,MED27,MED29,MED6,MED8": [("MED17", "MED22", "MED11", "MED18", "MED8", "MED6", "MED14", "MED29", "MED27")], \
# #                        "GGA1,IGF1R,IGF2,IGF2R,INSR,IRS1,SORL1": [("IGF1R", "INSR", "IRS1", "SORL1", "GGA1", "IGF2R", "IGF2")], \
# #                        "FOXJ3,MED1,MED15,MED16,MED23,MED24,MED25,MED9,PUM2": [("MED1", "MED9", "MED15", "FOXJ3", "MED23", "PUM2", "MED25", "MED16", "MED24")], \
# #                        "EFNA1,EFNA5,EFNB2,EFNB3,EPHA2,EPHA3,EPHA4": [("EFNA5", "EPHA2", "EFNA1", "EFNB2", "EPHA4", "EFNB3", "EPHA3")], \
# #                        "E2F1,E2F3,E2F4,E2F5,E2F6,TFDP1,TFDP2": [("E2F1", "E2F4", "E2F5", "TFDP1", "E2F6", "E2F3", "TFDP2")], \
# #                        "DMTN,DYNC1H1,DYNC1I2,DYNC1LI2,DYNLRB1,PROCR": [("DYNC1I2", "DYNLRB1", "PROCR", "DMTN", "DYNC1LI2", "DYNC1H1")], \
# #                        "CTCF,MAU2,NIPBL,RAD21,SMC1A,SMC3,STAG1,STAG2": [("MAU2", "NIPBL", "RAD21", "STAG2", "CTCF", "STAG1", "SMC1A", "SMC3")], \
# #                        "CPSF3L,INTS1,INTS2,INTS4,INTS5,INTS6,INTS7,INTS8,INTS9": [("INTS1", "INTS6", "INTS2", "INTS9", "CPSF3L", "INTS4", "INTS7", "INTS5", "INTS8")], \
# #                        "CPSF1,CPSF2,CPSF3,CPSF4,CSTF3,FIP1L1,WDR33": [("CSTF3", "CPSF1", "CPSF2", "CPSF3", "CPSF4", "FIP1L1", "WDR33")], 
# #                        "CNR1,CXCR2,GNAI1,GNAO1,GNB1,RGS1,RGS16,SMO": [("CNR1", "CXCR2", "GNB1", "GNAO1", "RGS16", "SMO", "RGS1", "GNAI1")], \
# #                        "CNOT1,CNOT2,CNOT3,CNOT7,RQCD1,ZFP36": [("CNOT1", "CNOT2", "CNOT3", "CNOT7", "RQCD1", "ZFP36")], \
# #                        "CDK12,CRBN,DCAF15,DDA1,DDB1,DDB2,SAMHD1,VPRBP": [("DDB1", "DDB2", "VPRBP", "SAMHD1", "CDK12", "CRBN", "DCAF15", "DDA1")], \
# #                        "CASC3,CWC22,CWC27,EIF4A3,NHP2L1,PRPF6,SMU1,SNRNP200,TXNL4A": [("CASC3", "EIF4A3", "CWC22", "CWC27", "SNRNP200", "TXNL4A", "PRPF6", "NHP2L1", "SMU1")], \
# #                        "C3,CD46,CD55,CFB,CFH,CFI,CR2,VSIG4": [("C3", "CD55", "CFB", "CFH", "CFI", "CD46", "CR2", "VSIG4")], \
# #                        "BAK1,BAX,BCL2,BCL2L1,BCL2L11,BID,BOK,HUWE1,MCL1": [("BCL2", "BAX", "BCL2L11", "BCL2L1", "BID", "BAK1", "BOK", "HUWE1", "MCL1")], \
# #                        "AEBP2,EED,ELMSAN1,EZH1,EZH2,HDAC1,JARID2,MTA1,PHF6,RBBP4,SUZ12": [("EZH1", "EED", "EZH2", "ELMSAN1", "HDAC1", "MTA1", "PHF6", "RBBP4", "SUZ12", "AEBP2", "JARID2")],\
# #                        "ATP6V0D1,ATP6V1A,ATP6V1B2,ATP6V1D,ATP6V1E1,ATP6V1F,ATP6V1G1": [("ATP6V1A", "ATP6V1B2", "ATP6V1E1", "ATP6V1G1", "ATP6V0D1", "ATP6V1D", "ATP6V1F")]
#                       }
# module2name = {"AEBP2,EED,ELMSAN1,EZH1,EZH2,HDAC1,JARID2,MTA1,PHF6,RBBP4,SUZ12": [("PRC2", "top")], \
#                "ATP6V0D1,ATP6V1A,ATP6V1B2,ATP6V1D,ATP6V1E1,ATP6V1F,ATP6V1G1": [("V-ATPase", "top")], \
#                "BAK1,BAX,BCL2,BCL2L1,BCL2L11,BID,BOK,HUWE1,MCL1": [("BCL2 family", "top")], \
#                "C3,CD46,CD55,CFB,CFH,CFI,CR2,VSIG4": [("Complement system", "top")], \
#                "CASC3,CWC22,CWC27,EIF4A3,NHP2L1,PRPF6,SMU1,SNRNP200,TXNL4A": [("Spliceosome", "top")], \
#                "CDK12,CRBN,DCAF15,DDA1,DDB1,DDB2,SAMHD1,VPRBP": [("DCX E3 ligase", "top")], \
#                "CNOT1,CNOT2,CNOT3,CNOT7,RQCD1,ZFP36": [("CCR4-NOT complex", "top")], \
#                "CNR1,CXCR2,GNAI1,GNAO1,GNB1,RGS1,RGS16,SMO": [("G-protein signaling", "top")], \
#                "CPSF1,CPSF2,CPSF3,CPSF4,CSTF3,FIP1L1,WDR33": [("CPSF complex", "top")], \
#                "CPSF3L,INTS1,INTS2,INTS4,INTS5,INTS6,INTS7,INTS8,INTS9": [("Integrator complex", "top")], \
#                "CTCF,MAU2,NIPBL,RAD21,SMC1A,SMC3,STAG1,STAG2": [("Cohesin function", "top")], \
#                "DMTN,DYNC1H1,DYNC1I2,DYNC1LI2,DYNLRB1,PROCR": [("Dynein-1 complex", "top")], \
#                "E2F1,E2F3,E2F4,E2F5,E2F6,TFDP1,TFDP2": [("E2F-dependent\ntranscription", "top")], \
#                "EFNA1,EFNA5,EFNB2,EFNB3,EPHA2,EPHA3,EPHA4": [("Eph/ephrin signaling", "top")], \
#                "FOXJ3,MED1,MED15,MED16,MED23,MED24,MED25,MED9,PUM2": [("Mediator complex", "top")], \
#                "GGA1,IGF1R,IGF2,IGF2R,INSR,IRS1,SORL1": [("Insulin/IGF signaling", "top")], \
#                "MED11,MED14,MED17,MED18,MED22,MED27,MED29,MED6,MED8": [("Mediator complex", "top")], \
#                "RIOK1,RIOK2,RPS14,RPS16,RPS18,RPS25,RPS28,RPS5": [("40S ribosome\nbiogenesis", "top")], \
#                "CAB39,MST4,PDCD10,STK24,STK25,STRADA": [("Serine/threonine\nkinases", "top")], \
#                "ZKSCAN4,ZNF213,ZNF446,ZNF496,ZSCAN12,ZSCAN21": [("Zinc finger\nproteins", "top")], \
#                "CUL3,FBXL17,H2AFY,KEAP1,NFE2L2,SPOP": [("BTB E3 ligase", "top")], \
#                "CHRAC1,DRAP1,NFYA,NFYB,NFYC,POLE3": [("NF-Y", "top")], \
#                "ADM,ADM2,CALCRL,RAMP1,RAMP2,RAMP3": [("CALCRL/RAMPs complex", "top")], \
#                "ARHGDIA,ARHGDIB,PAK4,RAC2,RAC3,RHOH": [("Rho GTPase\nsignaling", "top")], \
#                "CYC1,MT-CYB,UQCR10,UQCRC1,UQCRC2,UQCRFS1": [("Complex III", "top")], \
#                "MCM2,MCM3,MCM4,MCM5,MCM6,MCM7": [("MCM2-7 complex", "top")], \
#                "MMP14,MMP2,MMP9,THBS1,TIMP1,TIMP2": [("Metalloproteinase\nactivity", "top")], \
#                "EIF2B1,EIF2B2,EIF2B4,EIF2B5,EIF2S1,EIF2S3": [("eIF2 activity", "top")], \
#                "CDC73,CTR9,LEO1,PAF1,RTF1,WDR61": [("PAF1 complex", "top")], \
#                "MZT1,REPS1,TUBGCP2,TUBGCP3,TUBGCP6,WWC2": [(r"$\gamma$" + "-tubulin complex", "top")], \
#                "PPP2CA,PPP2R1A,PPP2R1B,PPP2R2A,PPP2R5A,PPP2R5C,PPP2R5D,STRN3": [("PP2A", "top")], \
#                "CCT2,CCT3,CCT4,CCT5,CCT6A,CCT7,CCT8,PFDN2,PFDN4,TCP1,VBP1": [("TRiC", "top")], \
#                "FKBP3,FLCN,FNIP2,MAPKAP1,MTOR,RICTOR,RPTOR,RRAGA,RRAGC,SLC38A9": [("mTOR signaling", "top")], \
#                "NDUFA10,NDUFA5,NDUFA7,NDUFS1,NDUFS2,NDUFS3,NDUFS4,NDUFS7,NDUFS8,NDUFV1,NDUFV3": [("Complex I", "top")], \
#                "AQR,BUD31,CDC40,CDC5L,CRNKL1,PPIE,PPIL1,RBM22,SNRNP40,SYF2,XAB2": [("Spliceosome", "top")], \
#                "BCL6,BCOR,BCORL1,DOT1L,KDM2B,MLLT3,PCGF1,PCGF3,PCGF5,PCGF6,RING1,RNF2,RYBP": [("Non-canonical PRC1", "top")], \
#                "CREBBP,FOSB,FOSL2,IRF3,JUN,JUNB,JUND,KANSL1,KMT2A,MEN1,MTDH,MYB,NCOA1,NCOA3,NCOR1,NR3C1,PABPC1,PAIP1,PPARG,PSIP1,RXRB,SND1,WDR5": [("Transactivation", "top")], \
#                "ARID1A,ARNT,ARNTL,BRD4,CBFB,COPS2,COPS3,COPS4,COPS7B,DPF2,EGLN1,EPAS1,GPS1,HIF1A,NPAS2,RUNX1,SMARCA4,SMARCC2,SMARCE1,SOCS2,TCEB1,TCEB2,VHL,WHSC1L1": [("Ubiquitination pathway", "top")], \
#                "ICT1,MRP63,MRPL10,MRPL11,MRPL13,MRPL14,MRPL15,MRPL16,MRPL18,MRPL19,MRPL20,MRPL23,MRPL24,MRPL27,MRPL3,MRPL30,MRPL32,MRPL35,MRPL38,MRPL39,MRPL4,MRPL40,MRPL41,MRPL42,MRPL46,MRPL47,MRPL48,MRPL51,MRPL52,MRPL53,MRPL54,MRPS18A,MRPS30": [("Mitochondrial ribosome", "top")], \
#                "ALB,APC,B2M,BTRC,CD1A,CD1B,CD1C,CD1D,CD8A,CDH5,COX6A1,COX6B1,COX6C,CTNNB1,CTNNBIP1,FBXW11,FCGRT,HFE,HLA-A,HLA-B,HLA-C,HLA-DRA,HLA-DRB1,HLA-DRB5,HLA-E,HLA-F,HLA-G,HMHA1,IDH2,IRS2,JUP,KDR,LILRB2,LSP1,MR1,MT-CO1,MT-CO2,MT-CO3,NR5A2,PDIA3,PPP1CB,PPP1R16B,PPP1R3B,PYGB,PYGL,SYNM,TAPBP,TAPBPL,VEGFC": [("Antigen presentation", "top"), ("Complex IV", "top"), ("Wnt signaling", "top")], \
#                "ASF1A,CRCP,ERCC2,ERCC3,GTF2A1,GTF2A2,GTF2B,GTF2E1,GTF2E2,GTF2F1,GTF2F2,GTF2H1,GTF2H3,GTF2H4,HIRA,POLR1D,POLR2A,POLR2B,POLR2F,POLR2H,POLR2I,POLR2K,POLR2L,POLR3A,POLR3B,POLR3C,POLR3D,POLR3E,POLR3F,POLR3H,POLR3K,TAF1,TAF10,TAF11,TAF2,TAF3,TAF4,TAF5,TAF6,TAF7,TAF8,TAF9,TBP": [("Transcription initiation and elongation", "top")], \
#                "CAND1,CUL1,CUL4A,CUL5,DCUN1D1,GLMN,NAE1,NEDD8,RBX1,RNF7,UBA3,UBE2M": [("haha", "top")], \
#                "PSMA2,PSMA4,PSMA5,PSMA6,PSMA7,PSMB10,PSMB7,PSMC1,PSMC2,PSMC3,PSMC4,PSMC6,PSMD1,PSMD11,PSMD12,PSMD13,PSMD14,PSMD3,PSMD6,PSMD7,PSMD8,PSME4,SHFM1": [("Proteasome", "top")], \
#                "CASP1,CASP4,ELF4,GSDMD,IL18,IL18R1,PYCARD": [("Inflammatory response", "top")], \
#               }
# gene2secondary = {"HLA-A": {"LSP1", "IRS2", "SYNM", "NR5A2","KDR","VEGFC"}, 
#                   "TAPBP": {"PDIA3"}, \
#                   "HLA-B": {"IDH2","PPP1CB","PPP1R16B","PPP1R3B","PYGB","PYGL"}, \
#                   "B2M": {"ALB","FCGRT"}}
# module2figsize = {"AEBP2,EED,ELMSAN1,EZH1,EZH2,HDAC1,JARID2,MTA1,PHF6,RBBP4,SUZ12_1": (9, 10), \
#                   "ATP6V0D1,ATP6V1A,ATP6V1B2,ATP6V1D,ATP6V1E1,ATP6V1F,ATP6V1G1_1": (8, 10), \
#                   "BAK1,BAX,BCL2,BCL2L1,BCL2L11,BID,BOK,HUWE1,MCL1_1": (8, 10), \
#                   "C3,CD46,CD55,CFB,CFH,CFI,CR2,VSIG4_1": (8, 10), \
#                   "CASC3,CWC22,CWC27,EIF4A3,NHP2L1,PRPF6,SMU1,SNRNP200,TXNL4A_1": (8, 10), \
#                   "CDK12,CRBN,DCAF15,DDA1,DDB1,DDB2,SAMHD1,VPRBP_1": (8, 10), \
#                   "CNOT1,CNOT2,CNOT3,CNOT7,RQCD1,ZFP36_1": (8, 10), \
#                   "CNR1,CXCR2,GNAI1,GNAO1,GNB1,RGS1,RGS16,SMO_1": (8, 10),\
#                   "CPSF1,CPSF2,CPSF3,CPSF4,CSTF3,FIP1L1,WDR33_1": (8, 10), 
#                   "CPSF3L,INTS1,INTS2,INTS4,INTS5,INTS6,INTS7,INTS8,INTS9_1": (8, 10), \
#                   "CTCF,MAU2,NIPBL,RAD21,SMC1A,SMC3,STAG1,STAG2_1": (8, 10), \
#                   "DMTN,DYNC1H1,DYNC1I2,DYNC1LI2,DYNLRB1,PROCR_1": (8, 10), \
#                   "E2F1,E2F3,E2F4,E2F5,E2F6,TFDP1,TFDP2_1": (8, 10), \
#                   "EFNA1,EFNA5,EFNB2,EFNB3,EPHA2,EPHA3,EPHA4_1": (8, 10), \
#                   "FOXJ3,MED1,MED15,MED16,MED23,MED24,MED25,MED9,PUM2_1": (8, 10), \
#                   "GGA1,IGF1R,IGF2,IGF2R,INSR,IRS1,SORL1_1": (8, 10), \
#                   "MED11,MED14,MED17,MED18,MED22,MED27,MED29,MED6,MED8_1": (8, 10), \
#                   "RIOK1,RIOK2,RPS14,RPS16,RPS18,RPS25,RPS28,RPS5_1": (8, 10), \
#                   "CAB39,MST4,PDCD10,STK24,STK25,STRADA_1": (8, 10), \
#                   "ZKSCAN4,ZNF213,ZNF446,ZNF496,ZSCAN12,ZSCAN21_1": (8, 10), \
#                   "CUL3,FBXL17,H2AFY,KEAP1,NFE2L2,SPOP_1": (8, 10), \
#                   "CHRAC1,DRAP1,NFYA,NFYB,NFYC,POLE3_1": (8, 10), \
#                   "ADM,ADM2,CALCRL,RAMP1,RAMP2,RAMP3_1": (8, 10), \
#                   "ARHGDIA,ARHGDIB,PAK4,RAC2,RAC3,RHOH_1": (8, 10), \
#                   "CYC1,MT-CYB,UQCR10,UQCRC1,UQCRC2,UQCRFS1_1": (8, 10), \
#                   "MCM2,MCM3,MCM4,MCM5,MCM6,MCM7_1": (8, 10), \
#                   "MMP14,MMP2,MMP9,THBS1,TIMP1,TIMP2_1": (8, 10), \
#                   "EIF2B1,EIF2B2,EIF2B4,EIF2B5,EIF2S1,EIF2S3_1": (8, 10), \
#                   "CDC73,CTR9,LEO1,PAF1,RTF1,WDR61_1": (8, 10), \
#                   "MZT1,REPS1,TUBGCP2,TUBGCP3,TUBGCP6,WWC2_1": (8, 10), \
#                   "PPP2CA,PPP2R1A,PPP2R1B,PPP2R2A,PPP2R5A,PPP2R5C,PPP2R5D,STRN3_1": (8, 10), \
#                   "CCT2,CCT3,CCT4,CCT5,CCT6A,CCT7,CCT8,PFDN2,PFDN4,TCP1,VBP1_1": (9, 10), \
#                   "FKBP3,FLCN,FNIP2,MAPKAP1,MTOR,RICTOR,RPTOR,RRAGA,RRAGC,SLC38A9_1": (8, 10), \
#                   "NDUFA10,NDUFA5,NDUFA7,NDUFS1,NDUFS2,NDUFS3,NDUFS4,NDUFS7,NDUFS8,NDUFV1,NDUFV3_1": (8, 10), \
#                   "AQR,BUD31,CDC40,CDC5L,CRNKL1,PPIE,PPIL1,RBM22,SNRNP40,SYF2,XAB2_1": (8, 10), \
#                   "CAND1,CUL1,CUL4A,CUL5,DCUN1D1,GLMN,NAE1,NEDD8,RBX1,RNF7,UBA3,UBE2M_1": (10, 10), \
#                   "PSMA2,PSMA4,PSMA5,PSMA6,PSMA7,PSMB10,PSMB7,PSMC1,PSMC2,PSMC3,PSMC4,PSMC6,PSMD1,PSMD11,PSMD12,PSMD13,PSMD14,PSMD3,PSMD6,PSMD7,PSMD8,PSME4,SHFM1_1": (15, 15), \
#                   "CASP1,CASP4,ELF4,GSDMD,IL18,IL18R1,PYCARD_1": (8, 10), \
#                   "BCL6,BCOR,BCORL1,DOT1L,KDM2B,MLLT3,PCGF1,PCGF3,PCGF5,PCGF6,RING1,RNF2,RYBP_1": (10, 10), \
#                   "CREBBP,FOSB,FOSL2,IRF3,JUN,JUNB,JUND,KANSL1,KMT2A,MEN1,MTDH,MYB,NCOA1,NCOA3,NCOR1,NR3C1,PABPC1,PAIP1,PPARG,PSIP1,RXRB,SND1,WDR5_1": (16, 16), \
#                   "ARID1A,ARNT,ARNTL,BRD4,CBFB,COPS2,COPS3,COPS4,COPS7B,DPF2,EGLN1,EPAS1,GPS1,HIF1A,NPAS2,RUNX1,SMARCA4,SMARCC2,SMARCE1,SOCS2,TCEB1,TCEB2,VHL,WHSC1L1_1": (16, 16), \
#                   "ICT1,MRP63,MRPL10,MRPL11,MRPL13,MRPL14,MRPL15,MRPL16,MRPL18,MRPL19,MRPL20,MRPL23,MRPL24,MRPL27,MRPL3,MRPL30,MRPL32,MRPL35,MRPL38,MRPL39,MRPL4,MRPL40,MRPL41,MRPL42,MRPL46,MRPL47,MRPL48,MRPL51,MRPL52,MRPL53,MRPL54,MRPS18A,MRPS30_1": (20, 20), \
#                   "ASF1A,CRCP,ERCC2,ERCC3,GTF2A1,GTF2A2,GTF2B,GTF2E1,GTF2E2,GTF2F1,GTF2F2,GTF2H1,GTF2H3,GTF2H4,HIRA,POLR1D,POLR2A,POLR2B,POLR2F,POLR2H,POLR2I,POLR2K,POLR2L,POLR3A,POLR3B,POLR3C,POLR3D,POLR3E,POLR3F,POLR3H,POLR3K,TAF1,TAF10,TAF11,TAF2,TAF3,TAF4,TAF5,TAF6,TAF7,TAF8,TAF9,TBP_1": (28, 28), \
#                   "ALB,APC,B2M,BTRC,CD1A,CD1B,CD1C,CD1D,CD8A,CDH5,COX6A1,COX6B1,COX6C,CTNNB1,CTNNBIP1,FBXW11,FCGRT,HFE,HLA-A,HLA-B,HLA-C,HLA-DRA,HLA-DRB1,HLA-DRB5,HLA-E,HLA-F,HLA-G,HMHA1,IDH2,IRS2,JUP,KDR,LILRB2,LSP1,MR1,MT-CO1,MT-CO2,MT-CO3,NR5A2,PDIA3,PPP1CB,PPP1R16B,PPP1R3B,PYGB,PYGL,SYNM,TAPBP,TAPBPL,VEGFC_1": (18, 18), \
#                   "ALB,APC,B2M,BTRC,CD1A,CD1B,CD1C,CD1D,CD8A,CDH5,COX6A1,COX6B1,COX6C,CTNNB1,CTNNBIP1,FBXW11,FCGRT,HFE,HLA-A,HLA-B,HLA-C,HLA-DRA,HLA-DRB1,HLA-DRB5,HLA-E,HLA-F,HLA-G,HMHA1,IDH2,IRS2,JUP,KDR,LILRB2,LSP1,MR1,MT-CO1,MT-CO2,MT-CO3,NR5A2,PDIA3,PPP1CB,PPP1R16B,PPP1R3B,PYGB,PYGL,SYNM,TAPBP,TAPBPL,VEGFC_2": (8, 10), \
#                   "ALB,APC,B2M,BTRC,CD1A,CD1B,CD1C,CD1D,CD8A,CDH5,COX6A1,COX6B1,COX6C,CTNNB1,CTNNBIP1,FBXW11,FCGRT,HFE,HLA-A,HLA-B,HLA-C,HLA-DRA,HLA-DRB1,HLA-DRB5,HLA-E,HLA-F,HLA-G,HMHA1,IDH2,IRS2,JUP,KDR,LILRB2,LSP1,MR1,MT-CO1,MT-CO2,MT-CO3,NR5A2,PDIA3,PPP1CB,PPP1R16B,PPP1R3B,PYGB,PYGL,SYNM,TAPBP,TAPBPL,VEGFC_3": (8, 10), \
#                  }


# sys.path.append("/fs/cbsuhy01/storage/yz2296/cancer_hotspot_new/website/figures/figure_2")
# import figure2
# import importlib
# importlib.reload(figure2)

# output_path = "/local/storage/yz2296/cancer_hotspot_new/website/website/analysis/multiscale_map/"
# with open("/local/storage/yz2296/cancer_hotspot_new/website/website/database/metadata/uniprot2gene.json") as f:
#     uniprot2key = json.load(f)

# druggable_prots = set(pd.read_csv("/local/storage/yz2296/cancer_hotspot_new/website/website/analysis/cancer_drug_target_all_JZ.txt", sep = "\t")["UniProt ID"])
# druggable_genes = {uniprot2key[x] for x in druggable_prots if x in uniprot2key}
# cancer_genes = all_function_py3.get_cgc_genes(cgc_tier1 = "/local/storage/yz2296/cancer_hotspot_new/reference_data/CGC_0227_2022.tsv", cgc_tier2 = "/local/storage/yz2296/cancer_hotspot_new/reference_data/CGC_0227_2022.tsv")

# cluster_file = "/local/storage/yz2296/cancer_hotspot_new/website/website/pipeline/output/TCGA_19_signatures.txt"
# module_file = "/local/storage/yz2296/cancer_hotspot_new/website/website/pipeline/output/TCGA_19_subnetworks.txt"
# init_graph = "/local/storage/yz2296/cancer_hotspot_new/website/website/pipeline/output/TCGA_19/initial_state.graphml.gz"
# final_graph = "/local/storage/yz2296/cancer_hotspot_new/website/website/pipeline/output/TCGA_19/final_state.graphml.gz"
# cluster_icon = "/local/storage/yz2296/cancer_hotspot_new/website/figures/figure_3/cluster.png"
# lof_icon = "/local/storage/yz2296/cancer_hotspot_new/website/figures/figure_3/lof.png"

# G_init = nx.read_graphml(init_graph)
# G_init = nx.relabel_nodes(G_init, uniprot2key)
# G_final = nx.read_graphml(final_graph)
# G_final = nx.relabel_nodes(G_final, uniprot2key)
# for key in module2OrderedGenes:
#     for i,x in enumerate(module2OrderedGenes[key]):
#         G_init_sub = G_init.subgraph(set(x))
#         if key in module2name:
#             title, pos = module2name[key][i]
#         else:
#             title = None
#             pos = None
#         fig_width, fig_height = module2figsize[key + "_" + str(i + 1)]
#         fig, ax = plt.subplots(figsize = (fig_width, fig_height))
#         if key in {"ARID1A,ARNT,ARNTL,BRD4,CBFB,COPS2,COPS3,COPS4,COPS7B,DPF2,EGLN1,EPAS1,GPS1,HIF1A,NPAS2,RUNX1,SMARCA4,SMARCC2,SMARCE1,SOCS2,TCEB1,TCEB2,VHL,WHSC1L1", "ICT1,MRP63,MRPL10,MRPL11,MRPL13,MRPL14,MRPL15,MRPL16,MRPL18,MRPL19,MRPL20,MRPL23,MRPL24,MRPL27,MRPL3,MRPL30,MRPL32,MRPL35,MRPL38,MRPL39,MRPL4,MRPL40,MRPL41,MRPL42,MRPL46,MRPL47,MRPL48,MRPL51,MRPL52,MRPL53,MRPL54,MRPS18A,MRPS30", "ASF1A,CRCP,ERCC2,ERCC3,GTF2A1,GTF2A2,GTF2B,GTF2E1,GTF2E2,GTF2F1,GTF2F2,GTF2H1,GTF2H3,GTF2H4,HIRA,POLR1D,POLR2A,POLR2B,POLR2F,POLR2H,POLR2I,POLR2K,POLR2L,POLR3A,POLR3B,POLR3C,POLR3D,POLR3E,POLR3F,POLR3H,POLR3K,TAF1,TAF10,TAF11,TAF2,TAF3,TAF4,TAF5,TAF6,TAF7,TAF8,TAF9,TBP","ALB,APC,B2M,BTRC,CD1A,CD1B,CD1C,CD1D,CD8A,CDH5,COX6A1,COX6B1,COX6C,CTNNB1,CTNNBIP1,FBXW11,FCGRT,HFE,HLA-A,HLA-B,HLA-C,HLA-DRA,HLA-DRB1,HLA-DRB5,HLA-E,HLA-F,HLA-G,HMHA1,IDH2,IRS2,JUP,KDR,LILRB2,LSP1,MR1,MT-CO1,MT-CO2,MT-CO3,NR5A2,PDIA3,PPP1CB,PPP1R16B,PPP1R3B,PYGB,PYGL,SYNM,TAPBP,TAPBPL,VEGFC"}:
#             fontsize = 22
#         else:
#             fontsize = 28
#         figure2.shell_layout(cluster_file, G_init_sub, cluster_icon, lof_icon, ax, (-0.5 * fig_width, 0.5 * fig_width), \
#                              (-0.5 * fig_height, 0.5 * fig_height), ordered_genes = x, interval_ratio = 0.8, \
#                              edge_magnifier = 2, fontsize = fontsize, figure_padding = 1, \
#                              drug_target_genes = druggable_genes, known_cancer_genes = cancer_genes, \
#                              gene2secondary = gene2secondary, marker_zoom = 0.055, \
#                              title = title, title_pos = pos)
#         fig.savefig(output_path + key.split(",")[0] + "_" + str(i + 1) + ".png", transparent = True)
############################################################################################################################################


# fig_width = 95
# fig_height = 94
# painting_color = "#d9d9d9"
# text_color = "#252525"
# nucleus_color = "#f9fff3"
# cytosal_color = "#fffbf7"
# environment_color = "#f6faff"
# fig, ax = plt.subplots(figsize = (fig_width, fig_height))
# plt.xlim(0, fig_width)
# plt.ylim(- fig_height, 0)
# ax.fill([0,fig_width, fig_width, 0], [0, 0, -fig_height, -fig_height], environment_color) 
# ax.fill([0,80.5,80.5,0], [0,0,-79.5,-79.5], cytosal_color)
# ax.fill([0,71,71,0], [-79.5,-79.5,-89,-89], cytosal_color)
# circle = plt.Circle((71, -79.5), 10, edgecolor = "white", linewidth = 0, facecolor = cytosal_color, zorder = 1)
# ax.add_artist(circle)
# plt.text(17, -24.5, "Transcription regulation", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(5, -28.5, "Pre-mRNA\nprocessing", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(18.5, -28.5, "DNA replication\nand repair", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(46, -2.73, "Nuclear import\nof proteins", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(36.5, -34, "Cohesin\nturn over", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(51, -19.5, "Cell cycle\nregulation", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(80, -2, "Metastasis", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(67.5, -27, "Cytoskeletal\norganization", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(57.5, -35, "Membrane trafficking", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(33.5, -43.5, "Ubiquitin-Proteasome system", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(80, -34, "Immune responses", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(69, -73, "mRNA\ndegradation", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(11, -85.5, "Environment sensing", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(6, -46.5, "Circadian\nclock", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(54, -11.5, "Multitask", va = "center", ha = "center", fontsize = 72, fontweight = "heavy", \
#		  color = "#252525")
# plt.text(89.5, -18, "Cell adhesion", va = "center", ha = "center", fontsize = 50, fontweight = "normal", \
#		  color = "#252525")
# plt.text(29, -87, "Electron transport chain", va = "center", ha = "center", fontsize = 50, fontweight = "normal", \
#		  color = "#252525")
# plt.text(43, -86, "Mitochondrial\nprotein synthesis", va = "center", ha = "center", fontsize = 50, fontweight = "normal", \
#		  color = "#252525")
# plt.plot([-1, 29.5], [-26, -26], color = text_color, linestyle = "dashed", linewidth = 5)
# plt.plot([11, 11], [-41, -26], color = text_color, linestyle = "dashed", linewidth = 5)
# plt.plot([80.5, 80.5], [-0.7, 0], color = painting_color, linewidth = 40, zorder = 2)
# plt.plot([80.5, 80.5], [-3.8, -3.3], color = painting_color, linewidth = 40, zorder = 2)
# plt.plot([80.5, 80.5], [-16.5, -13.5], color = painting_color, linewidth = 40, zorder = 2)
# plt.plot([80.5, 80.5], [-32.6, -24.5], color = painting_color, linewidth = 40, zorder = 2)
# plt.plot([80.5, 80.5], [-79.5, -35.4], color = painting_color, linewidth = 40, zorder = 2)
# plt.plot([-1, 71], [-89, -89], color = painting_color, linewidth = 40, zorder = 2)
# ax = all_function_py3.draw_arc(ax, (80.5 - 10, -89 + 10), 10, (-90, 0), edge_color = painting_color, lw = 40)
# circle = plt.Circle((11, -16), 26.7, edgecolor = painting_color, linewidth = 40, facecolor = "white")
# ax.add_artist(circle)
# circle = plt.Circle((11, -16), 25.3, edgecolor = painting_color, linewidth = 40, facecolor = nucleus_color)
# ax.add_artist(circle)
# circle = plt.Circle((37, -5), 4, edgecolor = painting_color, linewidth = 0, facecolor = "white") # Nuclear import of proteins
# ax.add_artist(circle)
# circle = plt.Circle((34.7, -26.45), 5.2, edgecolor = painting_color, linewidth = 0, facecolor = "white") # Cohesin turn over
# ax.add_artist(circle)
# module2pos = {"ARID1A_DPF2_SMARCB1_SMARCC1_SMARCC2_1": (4, -7), \
#			   "CHRAC1_DR1_DRAP1_GTF2B_NFYA_NFYB_NFYC_POLE3_1": (12.5, -5.5), \
#			   "ASNA1_ASPH_BAP1_CCNC_CDK8_EIF1AD_F10_GLYR1_GPX7_KDM1B_MGST3_PRKAB1_1": (21.5, -5.5), \
#			   "CLIP2_MAX_MLX_MNT_MXD1_MXD4_PLIN3_1": (29, -9), \
#			   "ASH2L_KMT2A_MEN1_POGZ_PSIP1_RBBP5_SETD1A_WDR5_1": (5, -16), \
#			   "BCOR_BCORL1_KDM2B_PCGF1_PCGF2_PCGF3_PCGF5_PCGF6_RING1_RYBP_1": (15.5, -15), \
#			   "E2F1_E2F4_E2F6_RB1_TFDP1_1": (25, -16), \
#			   "CPSF1_CPSF4_CSTF3_FIP1L1_WDR33_1": (5, -34), \
#			   "POLD2_POLD3_POLH_POLK_REV1_SOX12_1": (18, -34), \
#			   "KPNA1_KPNA3_KPNA4_KPNA6_MAT2A_MAT2B_NUP50_1": (36, -5), \
#			   "ANP32A_APEX1_FAM32A_GZMA_GZMK_LIG4_NSFL1C_PHAX_PNKP_POLB_RPL22_SET_SRP14_TARDBP_1": (44.5, -11), \
#			   "CTCF_NIPBL_PDS5A_PDS5B_RAD21_SMC1A_SMC3_STAG1_STAG2_WAPAL_1": (34, -26), \
#			   "CCND3_CDC34_CDKN1B_CRABP2_DCTD_EIF3K_EIF3L_MALSU1_TIMM13_UBE2V2_ZNRF1_1": (50, -26), \
#			   "SEPT10_SEPT11_SEPT5_SEPT6_SEPT7_SEPT8_SEPT9_1": (62, -9.5), \
#			   "TUBG1_TUBGCP2_TUBGCP3_TUBGCP5_TUBGCP6_1": (62, -17.5), \
#			   "ABI1_BRK1_CYFIP1_NCKAP1_WASF1_1": (70, -9.5), \
#			   "ANXA2_MYH10_MYH9_MYL9_S100A4_S100P_1": (70, -17.3), \
#			   "CD47_ITGA2_ITGAV_ITGB1_ITGB3_ITGB6_ITGB8_SIRPA_1": (78, -8.5), \
#			   "HSD3B7_ICAM1_ITGAL_ITGAM_ITGAX_ITGB2_1": (78, -20), \
#			   "ECH1_ECI2_EEF1B2_EEF1D_LCN2_MMP14_MMP2_MMP9_TIMP1_TIMP2_TPT1_1": (88, -10.5), \
#			   "CASP1_CASP4_ELF4_IL18_IL1RL1_IL33_NGF_NGFR_PYCARD_RIPK2_1": (72.5, -40.5), \
#			   "C3_CD46_CFB_CFH_CFI_1": (84.5, -44.5), \
#			   "AHCY_ALB_B2M_CBR3_CD1D_CDC25B_DDX3X_FBXW11_FCGRT_GRHPR_HLA-A_HLA-B_HLA-DPA1_HPRT1_HSCB_IDH2_IDH3A_IDH3B_IDH3G_IRS2_ISCU_LILRB1_MR1_PLP1_PYGB_PYGL_RAB35_RUSC2_SDHA_SDHAF2_SDHB_SUCLA2_SUCLG1_SUCLG2_TAPBPL_1": (72.5, -51), \
#			   "B4GALT7_METTL1_PTAR1_RAB37_RAB3A_RABGGTA_RABGGTB_WDR4_1": (52, -40), \
#			   "NAPB_NELFE_SNAP23_STX2_STX3_VAMP3_1": (60.5, -40), \
#			   "COPS2_COPS3_COPS4_COPS7B_COPS8_GPS1_PHYH_1": (32, -48), \
#			   "CAND1_CUL1_CUL4A_DCUN1D1_RFWD2_1": (40.5, -48), \
#			   "ACP5_PDCD5_PDRG1_PFDN2_PFDN4_PFDN6_PRMT8_TNFRSF10B_TNFSF10_VBP1_VHL_1":(21, -49.3), \
#			   "ARNTL_ARNTL2_CLOCK_EPAS1_NPAS2_1": (5.5, -51.5), \
#			   "JAG1_NOTCH1_NOTCH2_POGLUT1_RBPJ_1":(6, -68), \
#			   "CRELD2_GGA1_IDE_IGF2_IGF2R_IGFBP3_NDUFB11_SORL1_1": (6, -76), \
#			   "LAMTOR1_LAMTOR4_RPTOR_RRAGA_RRAGC_SLC38A9_1": (13.5, -72), \
#			   "AHCY_ALB_B2M_CBR3_CD1D_CDC25B_DDX3X_FBXW11_FCGRT_GRHPR_HLA-A_HLA-B_HLA-DPA1_HPRT1_HSCB_IDH2_IDH3A_IDH3B_IDH3G_IRS2_ISCU_LILRB1_MR1_PLP1_PYGB_PYGL_RAB35_RUSC2_SDHA_SDHAF2_SDHB_SUCLA2_SUCLG1_SUCLG2_TAPBPL_2": (55, -61), \
#			   "FANCF_FANCG_FCER1A_OLFM2_PDHA1_PDHB_PDK3_RAB11B_RAB11FIP2_RAB7A_TOMM22_TOMM7_1": (41, -63), \
#			   "MRPL2_MRPL22_MRPS7_MRPS9_TUFM_1": (47, -72), \
#			   "CYC1_NDUFA9_NDUFAF3_NDUFB10_NDUFS2_NDUFS7_NDUFS8_1": (30, -64), \
#			   "ATP5F1_ATP5H_ATP5L_NDUFA12_NDUFB9_1": (26, -71), \
#			   "ICT1_MRPL10_MRPL12_MRPL37_MRPL40_MRPL47_MRPL50_1": (37, -74), \
#			   "CNOT1_CNOT2_CNOT3_CNOT6L_CNOT8_RQCD1_1": (68, -78)
#			  }
# for key in module2pos:
#	 x, y = module2pos[key]
#	 fig, ax = all_function_py3.display_image(output_path + key + ".png", ax, x, y, zoom = 1)
# plt.axis("off")
# plt.margins(0,0)
# fig.savefig("/local/storage/yz2296/cancer_hotspot_new/website/figures/figure_3/figure3.pdf", bbox_inches = "tight")
# fig.savefig("/local/storage/yz2296/cancer_hotspot_new/website/figures/figure_3/figure3.png", bbox_inches = "tight")


############################################################################################################################################
# figure 4b
############################################################################################################################################
# gene = "MAX"
# maf_file = "/local/storage/yz2296/cancer_hotspot_new/data/version_5/Mutation_data/Pan_cancer/Zhang_Mutations_obsolete.txt"
# cluster_file = "/local/storage/yz2296/cancer_hotspot_new/website/website/download/pancan_clusters.txt"
# output = "/local/storage/yz2296/cancer_hotspot_new/website/figures/figure_4/" + gene + "_LoF.png"
# ax = figure2.mutation_distribution_lof(gene,uniprot2len[key2uniprot[gene]], maf_file, cluster_file, \
#											gene_color = "#70d3ff", output = output, unit_length = 0.08, \
#											block_height = 1, marker_weight = 2, unitmut_height = 0.8, fontsize = 28, \
#											marker_size = 1500)

# resource_path = "/local/storage/yz2296/cancer_hotspot_new/website/figures/figure_4/"
# G = nx.read_graphml("/local/storage/yz2296/cancer_hotspot_new/website/website/database/pancan/initial_state.graphml.gz")
# G = nx.relabel_nodes(G, uniprot2key)
# edge_radius = 0.23
# node_radius = 0.8
# node2pos = {"MAX": (-5, 0), "MLX": (5, 0), "MXD1": (0, 2), "MXD4": (0, 0), "MNT": (0, -2)}
# for key in node2pos:
#	 node2pos[key] = (node2pos[key], resource_path + key + "_cartoon.png")
# node2pos["MAX"] = (node2pos["MAX"][0], resource_path + "MAX_1R05.png")
# node2pos["MXD1"] = (node2pos["MXD1"][0], resource_path + "MXD1_1NLW.png")
# node2color = {x: ["#252525", node_radius, ("#252525", "normal")] for x in node2pos}
# node2color["MAX"][1] = 1
# specific_edge = {x: (edge_radius, None) for x in [("MAX", "MNT"), ("MNT", "MLX"), ("MLX", "MXD1"), ("MXD1", "MAX")]}
# specific_edge[("MAX", "MAX")] = (0.9 * node_radius, "left")
# specific_edge[("MNT", "MNT")] = (0.9 * node_radius, "bottom")
# specific_edge[("MXD1", "MXD1")] = (0.9 * node_radius, "top")
# specific_edge[("MLX", "MLX")] = (0.9 * node_radius, "right")
# G_sub = nx.Graph(G.subgraph(node2pos.keys()))
# for u,v in list(G_sub.edges()):
#	 if u == v:   
#		 G_sub.remove_edge(u, v)
# ax = figure2.Graph_layout_FineGrained(G_sub, node2pos, node2color, edge2info = specific_edge, edge_radius = edge_radius, \
#						   edge_color = "#252525", edge_magnifier = 2, node_rim = 0.5, zoom = 20, fontsize = 160)
# ax = all_function_py3.display_image(resource_path + "MAX-MXD1_1NLW.png", ax, -2.7, 1.6, zoom = 1)
# ax = all_function_py3.display_image(resource_path + "MAX-MXD4_cartoon.png", ax, -2.35, -0.1, zoom = 1)
# ax = all_function_py3.display_image(resource_path + "MAX-MNT_cartoon.png", ax, -2.7, -1.8, zoom = 1)
# ax = all_function_py3.display_image(resource_path + "MLX-MXD4_cartoon.png", ax, 2.5, 0.2, zoom = 1)
# ax = all_function_py3.display_image(resource_path + "MLX-MNT_cartoon.png", ax, 2.85, -1.55, zoom = 1)
# ax = all_function_py3.display_image(resource_path + "MAX_LoF.png", ax, -5.3, -0.5, zoom = 1.5)
# # ax = all_function_py3.display_image(resource_path + "MAX-MAX_1R05.png", ax, \
# #									 -2.9 - 1.5 * node_radius - 0.9 * node_radius, 0.2, zoom = 1)
# fig = plt.gcf()
# fig.savefig(resource_path + "4b.pdf", bbox_inches = "tight")
# fig.savefig(resource_path + "4b.png", bbox_inches = "tight", transparent = True)
# # output = "/local/storage/yz2296/cancer_hotspot_new/website/figures/figure_4/4b.pdf"
# # for u,v,info in G_sub.edges(data = True):
# #	 if u == v:
# #		 print(u,v,info)





############################################################################################################################################
# figure 4c
############################################################################################################################################
# gene2color = {"MXD1": "#70d3ff", "MAX": "#fbba1d", "MXD4": "#70d3ff", "MLX": "#fbba1d", "MNT": "#70d3ff"}
# cluster_file = "/local/storage/yz2296/cancer_hotspot_new/website/website/download/pancan_clusters.txt"
# mutation_file = "/local/storage/yz2296/cancer_hotspot_new/data/version_5/Mutation_data/Pan_cancer/Zhang_Mutations.txt"
# output = "/local/storage/yz2296/cancer_hotspot_new/website/figures/figure_4/4c.pdf"
# figure2.module_layout(gene2color, key2uniprot, cluster_file, mutation_file, output = output)