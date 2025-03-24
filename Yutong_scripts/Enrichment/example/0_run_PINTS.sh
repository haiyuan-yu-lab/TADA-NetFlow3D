#!/bin/bash

# DATA_DIR=/fs/cbsuhy01/storage/NetFlow_ASD/Yutong_output/Enrichment/data/enhancer
DATA_DIR=/fs/cbsuhy02/storage/yz2676/ASD/enhancer
mkdir -p ${DATA_DIR}

cell_types=("SC3" "CHTN22" "HN21" "CHTN13" "SC2" "SC8" "Microglia"
            "Neuron_4R" "LUHMES_Li" "ReNcellCX")

SC3_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/SC3/alignmentmerge_SC3a_and_SC3b_pl.bw
SC3_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/SC3/alignmentmerge_SC3a_and_SC3b_mn.bw

CHTN22_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/CHTN22/alignmentmerge_CHTN22a_and_CHTN22b_pl.bw
CHTN22_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/CHTN22/alignmentmerge_CHTN22a_and_CHTN22b_mn.bw

HN21_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/HN21/alignmentmerge_HN21a_and_HN21b_pl.bw
HN21_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/HN21/alignmentmerge_HN21a_and_HN21b_mn.bw

CHTN13_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/CHTN13/alignmentmerge_CHTN13a_and_CHTN13b_pl.bw
CHTN13_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/CHTN13/alignmentmerge_CHTN13a_and_CHTN13b_mn.bw

SC2_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/SC2/alignmentmerge_SC2a_and_SC2b_pl.bw
SC2_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/SC2/alignmentmerge_SC2a_and_SC2b_mn.bw

SC8_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/SC8/alignmentmerge_SC8a_and_SC8b_pl.bw
SC8_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/SC8/alignmentmerge_SC8a_and_SC8b_mn.bw

LUHMES_Li_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/LUHMES.Li/PROcap_luhmes_hg38_merged_alignment_pl.bw
LUHMES_Li_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/LUHMES.Li/PROcap_luhmes_hg38_merged_alignment_mn.bw

ReNcellCX_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/ReNcellCX/alignmentmerge_procap_rencellcx_bat1_rep1_and_procap_rencellcx_bat1_rep2_pl.bw
ReNcellCX_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/ReNcellCX/alignmentmerge_procap_rencellcx_bat1_rep1_and_procap_rencellcx_bat1_rep2_mn.bw

Neuron_4R_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/PROcap/Xiuqi/batch_3/Neuron_4R.yingying.all.PROcap_5pl.bw
Neuron_4R_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/PROcap/Xiuqi/batch_3/Neuron_4R.yingying.all.PROcap_5mn.bw

Microglia_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/PROcap/Xiuqi/batch_3/MG_WT.PROcap_5pl.bw
Microglia_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/PROcap/Xiuqi/batch_3/MG_WT.PROcap_5mn.bw

for cell_type in "${cell_types[@]}"
do
    OUTPATH=${DATA_DIR}/${cell_type}
    mkdir -p ${OUTPATH}
    PL_VAR="${cell_type}_PL_PATH"
    MN_VAR="${cell_type}_MN_PATH"
    python src/0_run_PINTS.py \
        --job_name ${cell_type} \
        --bw_pl ${!PL_VAR} \
        --bw_mn ${!MN_VAR} \
        --alpha 0.3 \
        --output_path ${OUTPATH}
done