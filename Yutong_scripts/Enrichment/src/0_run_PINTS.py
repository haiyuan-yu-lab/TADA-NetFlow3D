import os
import argparse
import pandas as pd
    
def run_PINTS(job_name, bw_pl, bw_mn, alpha = None, output_path = ".", promoter_file = "/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Promoters/promoters_1kb_tss_centered.coding.v43.bed", disable_sticks = False):
	'''
	promoter_file can also be /local/storage/yz2296/Burden_test_for_Autism/PROcap/Supplement/promoters_1kb_tss_centered.prtlncR.v34.bed
	alpha: PINTS v1.0 used alpha = 1; PINTS v1.1 used alpha = 0.3; By default alpha = 0.3. 
	'''
	print(f"Running job: {job_name}")
	
	if os.path.exists(output_path) == False:
		os.makedirs(output_path, exist_ok=True)

	command = "pints_caller --bw-pl " + bw_pl + " --bw-mn " + bw_mn + " --save-to " + output_path + \
	" --file-prefix " + job_name + " --thread 50"
	if alpha != None:
		command = command + " --alpha " + str(alpha)
	if disable_sticks == True:
		# command = command + " --disable-small"
		command = command + " --remove-sticks"
	os.system(command)

	command = "pints_boundary_extender --bw-pl " + bw_pl + " --bw-mn " + bw_mn + " --divergent-files " + \
	output_path + "/" + job_name + "_*_divergent_peaks.bed --bidirectional-files " + output_path + "/" + job_name + \
	"_*_bidirectional_peaks.bed --unidirectional-files " + output_path + "/" + job_name + \
	"_*_unidirectional_peaks.bed --save-to " + output_path + " --div-ext-left 200 --div-ext-right 200 \
	--single-ext-left 200 --single-ext-right 200 --promoter-bed " + promoter_file
	os.system(command)

	for suffix in ["", "_e"]:
		df1 = pd.read_csv(output_path + "/" + job_name + "_1_unidirectional_peaks_element_200bp" + suffix + ".bed", sep = "\t", header = None, dtype = str)
		df2 = pd.read_csv(output_path + "/" + job_name + "_1_divergent_peaks_element_200bp" + suffix + ".bed", sep = "\t", header = None, dtype = str)
		df3 = pd.read_csv(output_path + "/" + job_name + "_1_bidirectional_peaks_element_200bp" + suffix + ".bed", sep = "\t", header = None, dtype = str)
		pd.concat([df1, df2, df3])[[0, 1, 2]].drop_duplicates().rename(columns = {0: "chr", 1: "start", \
			2: "end"}).to_csv(output_path + "/all_peaks_element_200bp" + suffix + ".bed", sep = "\t", header = True, index = None)
		os.system("rm " + output_path + "/" + job_name + "_1_*_peaks_element_200bp" + suffix + ".bed")
	
	df1 = pd.read_csv(output_path + "/all_peaks_element_200bp.bed", sep = "\t")
	df2 = pd.read_csv(output_path + "/all_peaks_element_200bp_e.bed", sep = "\t")
	df2["flag"] = "enhancer"
	df = pd.merge(df1, df2, on = ["chr", "start", "end"], how = "left")
	df[df["flag"].isna()].drop(columns = ["flag"]).to_csv(output_path + "/all_peaks_element_200bp_p.bed", sep = "\t", header = True, index = None)
	# os.remove(output_path + "/all_peaks_element_200bp.bed")
	# os.system("rm " + output_path + "/" + job_name + "_1_*_peaks.bed")
	return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the run_pints function.")
    parser.add_argument("--job_name", type=str, help="Job name")
    parser.add_argument("--bw_pl", type=str, help="BW PL")
    parser.add_argument("--bw_mn", type=str, help="BW MN")
    parser.add_argument("--alpha", type=float, default=None, help="Alpha value (optional)")
    parser.add_argument("--output_path", type=str, default=".", help="Output path (optional)")

    args = parser.parse_args()

    run_PINTS(args.job_name, args.bw_pl, args.bw_mn, args.alpha, args.output_path)
