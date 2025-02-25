import pandas as pd
import numpy as np
from multiprocessing import Pool
import os
import argparse
import sys
sys.path.append("/home/yz2296/bin")
from collections import defaultdict
import logging
from os import path
import scipy.stats

description = '''
Scripts for mapping elements in a bed-like file1 (--file_to_be_filtered) to the elements in another bed-like file2 (--reference_file). \
The script will add one column called "Overlapped_elements" to file1 (--file_to_be_filtered) and write the output to --output. \
Each value in "Overlapped_elements" is a list of elements in file2 (--reference_file) separated by coma. These elements overlap with the corresponding\
element in file1 (--file_to_be_filtered).
'''

def whether_an_element_in_reference(chromosome, start_pos, end_pos, reference_file):
	'''
	[NOTE]: All coordinates are 0-indexed. Start positions are included and end positions are not.
	chromosome, start_pos, end_pos: genomic region covered by a target element
	reference_file: tab-separated file with at least three columns: chr, start and end. Must have a header line.
	'''
	df_ref = pd.read_csv(reference_file, sep = "\t", dtype = {"start": int, "end": int})[["chr", "start", "end"]].drop_duplicates()
	df_ref = df_ref[(df_ref["start"] >= 0) & (df_ref["end"] > 0)]
	df_ref_overlap = df_ref[(df_ref["chr"] == chromosome) & (df_ref["start"] < end_pos) & (df_ref["end"] > start_pos)]
	if df_ref_overlap.shape[0] > 0:
		return [chromosome, start_pos, end_pos, ",".join(df_ref_overlap.apply(lambda row: row["chr"] + ":" + str(row["start"]) + "-" + \
			str(row["end"]), axis = 1).tolist())]
	else:
		return -1

def multi_run_wrapper(args):
	return whether_an_element_in_reference(*args)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = description)
	parser.add_argument('-o','--output', required=True, type=str, help='Path to the output file')
	parser.add_argument('-f','--file_to_be_filtered', required=True, type=str, help='A bed-like file with at least 3 columns: chr, \
		start(0-indexed, included), end(0-indexed, not included). Must have a header line. Please keep this file as small as possible.')
	parser.add_argument('-r','--reference_file', required=True, type=str, help='A bed-like file with at least 3 columns: chr, \
		start(0-indexed, included), end(0-indexed, not included). Must have a header line. Please keep this file as small as possible.')
	parser.add_argument('-t','--threads', type=int, default=35, help='[OPTIONAL] The number of threads to use (by default 35)')
	args = parser.parse_args()

	def get_function_input(row, reference_file):
		function_input.append([row["chr"], row["start"], row["end"], reference_file])
		return True

	df = pd.read_csv(args.file_to_be_filtered, sep = "\t", dtype = {"start": int, "end": int})
	df_tmp = df[["chr", "start", "end"]].drop_duplicates()
	df_tmp = df_tmp[(df_tmp["start"] >= 0) & (df_tmp["end"] > 0)]
	function_input = []
	df_tmp["status"] = df_tmp.apply(lambda row: get_function_input(row, args.reference_file), axis = 1)

	p = Pool(args.threads)
	output = p.map(multi_run_wrapper, function_input)
	p.close()

	chr_list = []
	start_list = []
	end_list = []
	overlap_element_list = []
	for element in output:
		if element != -1:
			chr_list.append(element[0])
			start_list.append(element[1])
			end_list.append(element[2])
			overlap_element_list.append(element[3])
	df_retain = pd.DataFrame({"chr": chr_list, "start": start_list, "end": end_list, "Overlapped_elements": overlap_element_list})
	pd.merge(df, df_retain, on = ["chr", "start", "end"]).to_csv(args.output, sep = "\t", header = True, index = None)

