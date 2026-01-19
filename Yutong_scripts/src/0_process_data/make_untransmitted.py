import argparse
import subprocess
import sys
import pandas as pd
from typing import Iterable, List, Set, Tuple

def set_parser():
    parser = argparse.ArgumentParser(description="Extract untransmitted variants from parents in trios.")
    parser.add_argument("--trios", 
                        required=True, 
                        help="TSV: child\\tfather\\tmother\\tchild_pheno")
    
    parser.add_argument("--transmitted-tsv", 
                        required=True, 
                        help="tabix-indexed TSV.gz keyed by child sample id")
    
    parser.add_argument("--parent-tsv", 
                        required=True, 
                        help="tabix-indexed TSV.gz keyed by parent sample id")
    
    parser.add_argument("--out-body", 
                        required=True, 
                        help="Output VCF body (no header) path")
    
    parser.add_argument("--only-child", 
                        default="", 
                        help="Optional: process only this child")

    return parser.parse_args()

def parse_trios(path: str) -> List[Tuple[str, str, str, str]]:
    trios = pd.read_csv(path, sep="\t", header=None, 
                        names=["child", "father", "mother", "pheno"])
    trios_list = []
    for _, row in trios.iterrows():
        trios_list.append((row["child"], row["father"], row["mother"], row["pheno"]))
    return trios_list

def run_tabix(tsv_bgz: str, chrom: str) -> Iterable[str]:
    """
    Yields lines from a tabix-indexed TSV.bgz for a given 'chrom' (here: sample id).
    Assumes file was indexed with -s 1 -b 2 -e 3 and col1 is sample id.
    """
    # tabix will exit nonzero if chrom not present; treat as empty.
    cmd = ["tabix", tsv_bgz, chrom]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        return []
    return p.stdout.splitlines()

def varkey_from_vcf_fields(fields: List[str]) -> str:
    # fields: CHROM POS ID REF ALT QUAL FILTER INFO ...
    # Use CHROM:POS:REF:ALT.
    return f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"

def extract_transmitted_keys(transmitted_tsv_bgz: str, child: str) -> Tuple[Set[str], bool, int]:
    """
    Fetch transmitted variants for this child (from inherited TSV).
    Returns (keys, present, n_lines), where present=False means child_id not found in TSV.
    """
    lines = list(run_tabix(transmitted_tsv_bgz, child))
    if not lines:
        return set(), False, 0

    keys: Set[str] = set()
    for line in lines:
        parts = line.split("\t")
        vcf_fields = parts[3:]  # original VCF columns start here
        if len(vcf_fields) < 8:
            raise ValueError(f"Unexpected VCF line with <8 fields: {line}")
        keys.add(varkey_from_vcf_fields(vcf_fields))
    return keys, True, len(lines)

def emit_untransmitted_for_parent(
    parent_lines: List[str],
    child_id: str,
    child_pheno: str,
    transmitted: Set[str],
    out
    ) -> int:
    """
    Given tabix-returned lines for ONE parent, write untransmitted ones.
    Appends INFO fields: OFFSPRING, STATUS.
    Returns number of records written.
    """
    written = 0
    for line in parent_lines:
        parts = line.split("\t")
        vcf_fields = parts[3:]  # original VCF columns
        if len(vcf_fields) < 8:
            raise ValueError(f"Unexpected VCF line with <8 fields: {line}")

        k = varkey_from_vcf_fields(vcf_fields)
        if k in transmitted:
            continue

        info = vcf_fields[7]
        info += f";OFFSPRING={child_id};STATUS={child_pheno}"
        vcf_fields[7] = info
        out.write("\t".join(vcf_fields) + "\n")
        written += 1
    return written

def main():
    args = set_parser()

    trios = parse_trios(args.trios)

    skipped_child_absent_in_tx = 0
    skipped_parent_missing = 0   # strict: either father missing OR mother missing -> skip
    processed_trios = 0
    tx_lines_total = 0
    tx_keys_total = 0
    out_records_father = 0
    out_records_mother = 0
    out_records_total = 0

    children_processed = set()

    with open(args.out_body, "w") as out:
        for (child, father, mother, pheno) in trios:
            if args.only_child and child != args.only_child:
                continue

            transmitted, child_present, n_tx_lines = extract_transmitted_keys(args.transmitted_tsv, child)
            if not child_present:
                skipped_child_absent_in_tx += 1
                continue

            # father / mother variants not transmitted
            father_lines = list(run_tabix(args.parent_tsv, father))
            mother_lines = list(run_tabix(args.parent_tsv, mother))
            if not father_lines or not mother_lines:
                skipped_parent_missing += 1
                continue

            # Passed gates -> process
            processed_trios += 1
            children_processed.add(child)

            tx_lines_total += n_tx_lines
            tx_keys_total += len(transmitted)

            wf = emit_untransmitted_for_parent(father_lines, child, pheno, transmitted, out)
            wm = emit_untransmitted_for_parent(mother_lines, child, pheno, transmitted, out)

            out_records_father += wf
            out_records_mother += wm
            out_records_total += (wf + wm)

    # ---- summary output to stderr ----
    print("\n=== Untransmitted (STRICT TRIO) summary ===", file=sys.stderr)

    print(f"Skipped (child absent in transmitted TSV): {skipped_child_absent_in_tx}", file=sys.stderr)
    print(f"Skipped (missing father or mother in parent TSV): {skipped_parent_missing}", file=sys.stderr)

    print(f"Processed trios (passed both gates): {processed_trios}", file=sys.stderr)
    print(f"Unique children processed: {len(children_processed)}", file=sys.stderr)

    if processed_trios > 0:
        print(f"Transmitted TSV child-lines fetched (sum over processed trios): {tx_lines_total}", file=sys.stderr)
        print(f"Transmitted key count (sum over processed trios): {tx_keys_total}", file=sys.stderr)
        print(f"Avg transmitted keys per processed trio: {tx_keys_total / processed_trios:.2f}", file=sys.stderr)

    print(f"Output untransmitted records written: {out_records_total}", file=sys.stderr)
    print(f"  - from fathers: {out_records_father}", file=sys.stderr)
    print(f"  - from mothers: {out_records_mother}", file=sys.stderr)
    print("=== End summary ===\n", file=sys.stderr)

    
if __name__ == "__main__":
    main()
