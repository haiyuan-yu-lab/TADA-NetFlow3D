#!/bin/bash

fu_raw=/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/TADA_input_coding_FuRaw.txt
fu_scaled=/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/TADA_input_coding_FuScaled.txt

params=/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/sample_UCSF.txt
load_params() {
    local f="$1"
    while IFS=$'\t' read -r key value; do
        export "$key"="$value"
    done < "$f"
}
load_params "$params"

awk -v OFS="\t" -v FS="\t" \
  -v dn_case_scale="$(awk -v a="$n_dn_case_noncoding" -v b="$n_dn_case_coding" 'BEGIN{print a/b}')" \
  -v dn_con_scale="$(awk -v a="$n_dn_control_noncoding" -v b="$n_dn_control_coding" 'BEGIN{print a/b}')" \
  -v in_t_scale="$(awk -v a="$n_t_case_noncoding" -v b="$n_t_case_coding" 'BEGIN{print a/b}')" \
  -v in_u_scale="$(awk -v a="$n_t_control_noncoding" -v b="$n_t_control_coding" 'BEGIN{print a/b}')" \
    'NR==1{
        for (i=1; i<=NF; i++) h[i]=$i
        print $0
        next
    }
    {
        for (i=1; i<=NF; i++) {
            if (i<=16) {
                out[i]=$i
            } else{
                name=h[i]; val=$i
                if (name ~ /^dn/ && name !~ /\.sib$/) {
                    out[i] = int(val * dn_case_scale + 0.5)
                } else if (name ~ /^dn/ && name ~ /\.sib$/) {
                    out[i] = int(val * dn_con_scale + 0.5)
                } else if (name ~ /^in\.t/) {
                    out[i] = int(val * in_t_scale + 0.5)
                } else if (name ~ /^in\.u/) {
                    out[i] = int(val * in_u_scale + 0.5)
                } else {
                    printf("ERROR: unexpected column header after col16: col=%d header=%s\n", i, h) > "/dev/stderr"
                    exit 2
                }
            }
            printf "%s%s", out[i], (i<NF?OFS:ORS)
        }
    }' "$fu_raw" > "$fu_scaled"