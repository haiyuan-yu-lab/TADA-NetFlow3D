import argparse
import os
import logging

import pandas as pd
from helpers import consequence2varclass_coding, transform_variant, Exclude_Sib, Exclude_Mul, annotate_conservation_to_element, parse_annotated_Rev


def set_parser(parser):
    parser.add_argument("--dnv_vep_path", 
                        help="Path to DNV VEP output.",
                        type=str,
                        required=True,
                        )
    
    parser.add_argument("--trans_vep_path",
                        help="Path to transmitted VEP output.",
                        type=str,
                        required=True,
                        )
    
    parser.add_argument("--raw_dir",
                        help="Path to raw UCSF data.",
                        type=str,
                        required=True,
                        )
    
    parser.add_argument("--family",
                        help="Path to ASD family information.",
                        type=str,
                        required=True,
                        )
    
    parser.add_argument("--metadata",
                        help="Path to ASD metadata.",
                        type=str,
                        required=True,
                        )
    
    parser.add_argument("--out_path", 
                        help="Path to save output.",
                        type=str,
                        required=True,
                        )


### Preprocess de novo variants
def DNV_preprocess(ucsf_vep_path, ucsf_raw, asd_family, asd_metadata, outdir, coding_flag = "damaging"):
    """
    Preprocess de novo variants from annotated VEP output and raw data.
    
    Parameters:
        ucsf_vep_path: str
            Path to VEP output.
        ucsf_raw: str
            Path to raw UCSF data.
        asd_family: str
            Path to ASD family information.
        asd_metadata: str
            Path to ASD metadata.
        outdir: str
            Path to save output.
        coding_flag: 1. "damaging"; 2. "strict"; 3. "loose"; 4. None.
            1. "damaging": exclude LOF & damaging missense variants
            2. "strict": exclude variants with ANY coding consequences
            3. "loose": exclude variants with ONLY coding consequences
            4. None: do not exclude any coding variants
    """
    # Load and process VEP output
    UCSF_vep = pd.read_csv(ucsf_vep_path, sep = "\t", dtype = str)
    UCSF_vep["Variant_Classification"] = UCSF_vep["Consequence"].apply(consequence2varclass_coding)

    # load and process raw data
    UCSF_hg38 = pd.read_csv(ucsf_raw+"UCSF_20211001_dnms_WGS_hg38.tsv", sep = "\t", dtype = str)
    UCSF_hg38["CHROM"] = UCSF_hg38["CHROM"].apply(lambda x: x[3:])
    UCSF_hg38 = UCSF_hg38.rename(
        columns = {"CHROM": "chr", 
                   "POS":"pos", 
                   "REF":"ref", 
                   "ALT":"alt", 
                   "STATUS": "Pheno", 
                   "SAMPLE": "SampleID"}
                   )

    # extra processing to make pos/ref/alt the same format as VEP files
    UCSF_hg38["#Uploaded_variation"] = UCSF_hg38["ID"].apply(transform_variant)

    # Annotate raw data
    UCSF_merged = pd.merge(UCSF_hg38, UCSF_vep, on = "#Uploaded_variation")
    # print(UCSF_merged)

    # Exclude siblings with ASD risk
    roles = pd.read_csv(asd_family, index_col = None)
    UCSF_merged = Exclude_Sib(UCSF_merged, roles)

    # Exclude multiplex families
    UCSF_merged = Exclude_Mul(UCSF_merged, roles)

    # Exclude coding variants
    if coding_flag == "strict":
        coding_var = set(UCSF_merged[~UCSF_merged["Variant_Classification"].isna()]["#Uploaded_variation"])
        UCSF_merged = UCSF_merged[~UCSF_merged["#Uploaded_variation"].isin(coding_var)]
    elif coding_flag == "loose":
        UCSF_merged = UCSF_merged[UCSF_merged["Variant_Classification"].isna()] # any variants with noncoding consequences will be included.
    elif coding_flag == "damaging": 
        # rest of the tags are: missense_mutation, silent, in_frame_del, in_frame_ins; 
        # silent: synonymous mutations, in_frame_del: reading frame of the gene is preserved and not disrupted
        lof = {"Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Nonstop_Mutation","Translation_Start_Site","Splice_Site"}
        UCSF_lof = UCSF_merged[UCSF_merged["Variant_Classification"].isin(lof)]
        
        UCSF_lof.to_csv(outdir+"/UCSF_DNV_preprocessed_LOF.txt", sep = "\t", header = True, index = None)
        UCSF_merged[UCSF_merged["Variant_Classification"] == "Missense_Mutation"].to_csv(outdir+"/mis_tmp.txt", sep = "\t", header = True, index = None)
        
        # REVEL: rare missense variant pathogenicity prediction
        command = "python /fs/cbsuhyfs1/storage/resources/variant_effects/tools/Parse_REVEL.py -m "+outdir+"mis_tmp.txt -c grch38 -o " + outdir + "/mis_tmp_annotated.txt"
        os.system(command)

        UCSF_mis = pd.read_csv(outdir + "/mis_tmp_annotated.txt", sep = "\t").astype(str)
        UCSF_mis["parsed_revel"] = UCSF_mis["revel"].apply(parse_annotated_Rev)
        UCSF_Dmis = UCSF_mis[UCSF_mis["parsed_revel"] >= 0.5]
        UCSF_Dmis.to_csv(outdir + "/UCSF_DNV_preprocessed_Dmis.txt", sep = "\t", header = True, index = None)
        exclude = set(pd.concat([UCSF_Dmis, UCSF_lof])["#Uploaded_variation"])
        UCSF_merged = UCSF_merged[~UCSF_merged["#Uploaded_variation"].isin(exclude)]
        os.remove(outdir+"/mis_tmp.txt")
        os.remove(outdir+"/mis_tmp_annotated.txt")

    # Annotate conservation
    UCSF_merged["chr"] = UCSF_merged["chr"].apply(lambda x: "chr" + str(x))
    UCSF_merged["start"] = UCSF_merged["pos"].apply(lambda x: int(x) - 1)
    UCSF_merged["end"] = UCSF_merged.apply(lambda row: int(row["pos"]) + len(row["ref"]) - 1, axis = 1)
    UCSF_merged = annotate_conservation_to_element(UCSF_merged, outdir, None, chr_col = "chr", start_col = "start", end_col = "end", start_col_index = 0, flag = "max")

    # Annotate sex
    sex = pd.read_csv(asd_metadata, sep = "\t")
    UCSF_merged = pd.merge(UCSF_merged, sex[["SampleID", "sex"]].drop_duplicates(), on = "SampleID", how = "left")

    # Annotate gnomADg_AF
    output = outdir + "/UCSF_DNV_ExcludeCoding_" + coding_flag + ".txt"
    UCSF_merged.to_csv(output, sep = "\t", header = True, index = None)

    # Annotate gnomADg_AF (gnomAD version: r3.1)
    command = "python /fs/cbsuhyfs1/storage/resources/gnomad/tools/gnomAD_grab.py -m " + output + " -g hg38 -v r3.1 -s genome -a " + output
    os.system(command)
    return output


def Transmitted_preprocess(ucsf_vep_path, ucsf_raw, asd_family, asd_metadata, outdir, coding_flag = "damaging"):
    """
    Preprocess transmitted variants from raw data.

    Parameters:
        ucsf_vep_path: str
            Path to VEP output.
        ucsf_raw: str
            Path to raw UCSF data.
        asd_family: str
            Path to ASD family information.
        asd_metadata: str
            Path to ASD metadata.
        outdir: str
            Path to save output.
        coding_flag: 1. "damaging"; 2. "strict"; 3. "loose"; 4. None.
            1. "damaging": exclude LOF variants
            2. "strict": exclude variants with ANY coding consequences
            3. "loose": exclude variants with ONLY coding consequences
            3. None: do not exclude any coding variants
    """
    # load and process raw data
    Inherited = []
    for file in ["KOR-STAR.rare.inherited.tsv.gz", "SSC_WGS9209s.rare.inherited.tsv.gz", "SPARK_WGS_1.rare.inherited.tsv.gz", "SPARK_WGS_2.rare.inherited.tsv.gz", "SPARK_WGS_3.rare.inherited.tsv.gz"]:
        Inherited.append(pd.read_csv(ucsf_raw + file, sep = "\t", dtype = str))
    Inherited = pd.concat(Inherited).rename(columns = {"sample_id": "SampleID"}).drop(columns = ["gnomad_AF"]) # gnomAD version: r2.1.1 annotation in vep files
    Inherited[["chr", "pos"]] = pd.DataFrame(Inherited["locus"].apply(lambda x: x.split(":")).tolist(), index = Inherited.index)
    
    # Inherited files also need extra processing - vep files are of format (chr_pos_-/T), whereas
    # raw data in format (chr_pos_ref/alt, e.g. 1_123456_A/AT)
    Inherited["ID"] = Inherited["chr"] + ":" + Inherited["pos"] + ":" + Inherited["ref"] + ":" + Inherited["alt"]
    Inherited["#Uploaded_variation"] = Inherited["ID"].apply(transform_variant)

    Inherited["start"] = Inherited["pos"].apply(lambda x: int(x) - 1)
    Inherited["end"] = Inherited.apply(lambda x: int(x["pos"]) + len(x["ref"]) - 1, axis = 1)
    logging.warning("transmitted line original: " + str(Inherited[["#Uploaded_variation", "SampleID"]].drop_duplicates().shape[0]))

    # load variant annotation
    var_annot = pd.read_csv(ucsf_vep_path, sep = "\t", dtype = str)
    Inherited = pd.merge(Inherited, var_annot, on = "#Uploaded_variation")
    logging.warning("transmitted line after mapping: " + str(Inherited[["#Uploaded_variation", "SampleID"]].drop_duplicates().shape[0]))

    # Exclude siblings with ASD risk
    roles = pd.read_csv(asd_family, index_col = None)
    Inherited = Exclude_Sib(Inherited, roles)

    # Exclude coding variants
    Inherited["Variant_Classification"] = Inherited["Consequence"].apply(consequence2varclass_coding)
    if coding_flag == "strict":
        coding_var = set(Inherited[~Inherited["Variant_Classification"].isna()]["#Uploaded_variation"])
        Inherited = Inherited[~Inherited["#Uploaded_variation"].isin(coding_var)]
    elif coding_flag == "loose":
        Inherited = Inherited[Inherited["Variant_Classification"].isna()]
    elif coding_flag == "damaging":
        lof = {"Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Nonstop_Mutation","Translation_Start_Site","Splice_Site"}
        Inherited_lof = Inherited[Inherited["Variant_Classification"].isin(lof)]
        Inherited_lof.to_csv(outdir + "/UCSF_Transmitted_preprocessed_LOF.txt", sep = "\t", header = True, index = None)
        Inherited = Inherited[~Inherited["#Uploaded_variation"].isin(set(Inherited_lof["#Uploaded_variation"]))]

    # Annotate sex and phenotype
    sex = pd.read_csv(asd_metadata, sep = "\t")
    Inherited = pd.merge(Inherited, sex[["SampleID", "sex", "Pheno"]].drop_duplicates(), on = "SampleID", how = "left")

    # Write output
    output = outdir + "/UCSF_Transmitted_ExcludeCoding_" + coding_flag + ".txt"
    Inherited.to_csv(output, sep = "\t", header = True, index = None)

    # Annotate with gnomAD version: r3.1
    command = "python /fs/cbsuhyfs1/storage/resources/gnomad/tools/gnomAD_grab.py -m " + output + " -g hg38 -v r3.1 -s genome -a " + output
    os.system(command)

    return output


def main(args):
    # DNV_preprocess(args.dnv_vep_path, args.raw_dir, args.family, args.metadata, args.out_path)
    Transmitted_preprocess(args.trans_vep_path, args.raw_dir, args.family, args.metadata, args.out_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    set_parser(parser)
    args = parser.parse_args()
    main(args)
