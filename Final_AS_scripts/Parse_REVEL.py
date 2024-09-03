description = "Python Script to annotate REVEL score to mutations"

import os
import argparse
import pandas as pd
from os import path
import sys


file_path = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
	# user input
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-m','--Mutation_file', required = True, type = str, help = 'Path to the mutation file; Please include columns named ["chr", "pos", "ref", "alt"] with format{[12], [9077413], [A], [C]}')
    parser.add_argument('-c','--Human_reference_genome', required = True, type = str, help = 'Enter either hg19 or grch38')
    parser.add_argument('-o','--output', required = True, type = str, help = 'Path to the output file')
    parser.add_argument('-r','--REVEL_score', type = str, default = file_path + "/../REVEL/revel-v1.3_all_chromosomes.zip", help = '[OPTIONAL] Path to full REVEL score zipped file')
    args = parser.parse_args()

############################################################################################################################################
	# Check input files
	# -----------------
    if path.exists(args.Mutation_file) == False:
        sys.exit("No such file: " + args.Mutation_file)
    if args.Human_reference_genome not in {"hg19", "grch38"}:
        sys.exit("Please check --Human_reference_genome! Can only be set to hg19 or grch38.")
	############################################################################################################################################

    # load data
    genome = str(args.Human_reference_genome + "_pos")
    REVEL = pd.read_csv(args.REVEL_score, compression='zip', dtype = str)
    mutation_file = pd.read_csv(args.Mutation_file, sep = "\t", dtype = str)
    mutation_file = mutation_file.rename(columns = {"pos": genome})
    
    # user input might be redundant, extract_user_input = [[chr, pos, ref, alt]].drop_duplicates(), 
    extract_user_input = mutation_file[["chr", genome, "ref", "alt"]].drop_duplicates()
    df_merge = pd.merge(REVEL, extract_user_input, on = ["chr", genome, "ref", "alt"])
    df_final = {"chr": [], genome: [], "ref": [], "alt": [], "revel": []}
    for (chr, pos, ref, alt),df_pervariant in df_merge.groupby(["chr", genome, "ref", "alt"]):
        df_pervariant["info"] = df_pervariant.apply(lambda row: row["Ensembl_transcriptid"] + ":" + row["REVEL"], axis = 1)
        ALL_REVEL_SCORE = ",".join(df_pervariant["info"].tolist())
        df_final["chr"].append(chr)
        df_final[genome].append(pos)
        df_final["ref"].append(ref)
        df_final["alt"].append(alt)
        df_final["revel"].append(ALL_REVEL_SCORE)
    df_final = pd.DataFrame(df_final)
    final_merge = pd.merge(mutation_file, df_final, on =["chr", genome, "ref", "alt"], how = "left")
    final_merge.to_csv(args.output, sep = "\t", header = True, index = None)
