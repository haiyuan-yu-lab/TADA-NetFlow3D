#!/bin/env bash
out=manifests/procap_data.tsv
mkdir -p "$(dirname "$out")"

cell_types=("SC3" "CHTN22" "HN21" "CHTN13" "SC2" "SC8" "MG_WT"
            "Neuron_4R" "LUHMES_Li" "ReNcellCX" "C2" "C1" "SC1")

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

MG_WT_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/PROcap/Xiuqi/batch_3/MG_WT.PROcap_5pl.bw
MG_WT_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/PROcap/Xiuqi/batch_3/MG_WT.PROcap_5mn.bw

C2_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/HEK293/alignmentmerge_C2a_and_C2b_pl.bw
C2_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/HEK293/alignmentmerge_C2a_and_C2b_mn.bw

C1_PL_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/K562/alignmentmerge_C1a_and_C1b_pl.bw
C1_MN_PATH=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/PROcap_hg38/K562/alignmentmerge_C1a_and_C1b_mn.bw

SC1_PL_PATH=/fs/cbsuhy02/storage/yz2676/data/PROcap/SC1/SC1_PROcap.merged.5pl.bw
SC1_MN_PATH=/fs/cbsuhy02/storage/yz2676/data/PROcap/SC1/SC1_PROcap.merged.5mn.bw

{
  printf "cell_type\tbw_pl\tbw_mn\n"
  for ct in "${cell_types[@]}"; do
    pl_var="${ct}_PL_PATH"
    mn_var="${ct}_MN_PATH"
    pl="${!pl_var}"
    mn="${!mn_var}"
    if [[ -z "$pl" || -z "$mn" ]]; then
      echo "[Warning] missing path for $ct" >&2
      continue
    fi
    printf "%s\t%s\t%s\n" "$ct" "$pl" "$mn"
  done
} > "$out"
echo "wrote $out"
