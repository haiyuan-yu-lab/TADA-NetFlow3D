"""Assemble the manuscript TADA table and export reported network modules."""

from __future__ import annotations

from pathlib import Path
import csv
import gzip

import numpy as np
import pandas as pd

from combine_cohort_bf import q2p

FU_CODING_LOG10_BF_COLUMNS = {
    "l10BF_ASD_PTV": "ptv",
    "l10BF_ASD_misB": "misb",
    "l10BF_ASD_misA": "misa",
}
FU_CNV_LOG10_BF_COLUMNS = {"l10BF_ASD_DEL": "del", "l10BF_ASD_DUP": "dup"}

MODEL_COLUMNS = {
    "exc": "l10BF_EM_exc",
    "mic": "l10BF_EM_mic",
    "ast": "l10BF_EM_ast",
    "cortex": "l10BF_EM_cortex",
    "other_cortex": "l10BF_EM_other_cortex",
    "combined": "l10_BF_EM",
}


def bayesian_fdr(bayes_factors: pd.Series | np.ndarray, pi0: float = 0.95) -> np.ndarray:
    """Calculate FDR in input order using the manuscript posterior calculation."""
    bf = np.asarray(bayes_factors, dtype=float)
    order = np.argsort(-bf)
    bf_sorted = bf[order]
    pi_alt = 1.0 - pi0
    posterior_alt = pi_alt * bf_sorted / (1.0 - pi_alt + pi_alt * bf_sorted)
    posterior_null = 1.0 - posterior_alt
    fdr_sorted = np.cumsum(posterior_null) / np.arange(1, len(bf_sorted) + 1)
    fdr = np.empty_like(fdr_sorted)
    fdr[order] = fdr_sorted
    return fdr


def noncoding_log_bf(noncoding: pd.DataFrame) -> pd.Series:
    """Read noncoding evidence after flooring within each dataset and component."""
    if "bf_floor_separate" not in noncoding.columns:
        raise ValueError("Noncoding table lacks a bf_floor_separate column")
    values = pd.to_numeric(noncoding["bf_floor_separate"], errors="raise")
    if not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("Noncoding bf_floor_separate must be finite and positive")
    return np.log(values)


def combine_coding_and_noncoding(coding: pd.DataFrame, noncoding: pd.DataFrame) -> pd.DataFrame:
    """Add prepared coding and noncoding log BFs and floor the final heat at zero."""
    noncoding = noncoding.copy()
    noncoding["log_bf.em_separate"] = noncoding_log_bf(noncoding)
    merged = coding.merge(noncoding[["gene_id", "log_bf.em_separate"]], on="gene_id")
    merged["log_bf"] = merged["log_bf.coding"] + merged["log_bf.em_separate"]
    merged["bf"] = np.exp(merged["log_bf"])
    merged["fdr_addEM"] = bayesian_fdr(merged["bf"])
    merged["heat_score"] = merged["log_bf"].clip(lower=0)
    return merged


def assemble_tada_table(fu: pd.DataFrame, noncoding: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Apply the Figure 2 evidence combination and Table 8 column schema."""
    output = fu[["gene_id", "gene"]].rename(columns={"gene": "Hugo_symbol"}).copy()
    for model, column in MODEL_COLUMNS.items():
        table = noncoding[model].set_index("gene_id")
        bf = output["gene_id"].map(table["bf_floor_separate"])
        output[column] = bf.clip(lower=1).apply(np.log10)
    coding_columns = list(FU_CODING_LOG10_BF_COLUMNS) + list(FU_CNV_LOG10_BF_COLUMNS)
    output["l10_BF_coding"] = fu[coding_columns].sum(axis=1)
    output["l10_BF_coding_noncoding"] = output["l10_BF_coding"] + output["l10_BF_EM"]
    # Preserve the notebook's product of component BFs for FDR calculation.
    components = fu[coding_columns].copy()
    components["l10BF_ASD_EM"] = output["l10_BF_EM"]
    bf = (10 ** components).prod(axis=1, skipna=True)
    output["FDR_TADA_ASD_addEM"] = bayesian_fdr(bf, pi0=0.95)
    output["p_TADA_ASD_addEM"] = q2p(output["FDR_TADA_ASD_addEM"], output["gene_id"])
    return output


def check_numeric_table(actual: pd.DataFrame, expected_file: Path, key: str) -> None:
    """Compare identifiers and every reference column, allowing TSV float precision."""
    expected = pd.read_csv(expected_file, sep="\t")
    if actual[key].duplicated().any() or set(actual[key]) != set(expected[key]):
        raise ValueError(f"Gene identifiers differ from {expected_file.name}")
    actual = actual.set_index(key).loc[expected[key]]
    expected = expected.set_index(key)
    for column in expected.columns:
        if column not in actual:
            raise ValueError(f"Missing reference column: {column}")
        if pd.api.types.is_numeric_dtype(expected[column]):
            np.testing.assert_allclose(actual[column], expected[column], rtol=1e-8, atol=1e-12,
                                       equal_nan=True, err_msg=f"{expected_file.name}: {column}")
        elif actual[column].fillna("").tolist() != expected[column].fillna("").tolist():
            raise ValueError(f"{expected_file.name}: {column} differs")


def protein_sets(table: pd.DataFrame, column: str) -> list[tuple[str, ...]]:
    return sorted(tuple(sorted(str(s).split(","))) for s in table[column])


def export_modules(subnetwork_file: Path, id_mapping: Path, output_dir: Path,
                   expected_dir: Path, size_cutoff: int = 5) -> tuple[int, int]:
    """Select reported modules and map proteins as in the Figure 4 notebook."""
    all_modules = pd.read_csv(subnetwork_file, sep="\t")
    expected = pd.read_csv(expected_dir / "netflow3d_subnetworks.tsv", sep="\t")
    if protein_sets(all_modules, "Subnetwork_genes") != protein_sets(expected, "Subnetwork_genes"):
        raise ValueError("NetFlow3D subnetworks differ from Supplementary Table 10")
    modules = all_modules[all_modules["Subnetwork_size"] > size_cutoff].copy()
    names = {}
    opener = gzip.open if id_mapping.suffix == ".gz" else open
    with opener(id_mapping, "rt") as handle:
        for protein, kind, identifier in csv.reader(handle, delimiter="\t"):
            if kind == "Gene_Name":
                names[protein.split("-")[0]] = identifier.split(".")[0]
    records = []
    for row in modules.itertuples(index=False):
        proteins = row.Subnetwork_genes.split(",")
        missing = set(proteins) - names.keys()
        if missing:
            raise ValueError(f"Gene_Name mappings missing for: {', '.join(sorted(missing))}")
        records.append({"Proteins": ",".join(sorted(proteins)),
                        "Genes": ",".join(sorted({names[p] for p in proteins}))})
    records.sort(key=lambda r: (-len(r["Proteins"].split(",")), r["Proteins"]))
    reported = pd.DataFrame(records, columns=["Proteins", "Genes"])
    reference = pd.read_csv(expected_dir / "reported_modules.tsv", sep="\t")
    actual_pairs = sorted(zip(protein_sets_by_row(reported, "Proteins"),
                              protein_sets_by_row(reported, "Genes")))
    expected_pairs = sorted(zip(protein_sets_by_row(reference, "Proteins"),
                                protein_sets_by_row(reference, "Genes")))
    if actual_pairs != expected_pairs:
        raise ValueError("Reported module proteins or genes differ from Supplementary Table 11")
    genes = sorted({g for s in reported["Genes"] for g in s.split(",")})
    reported.to_csv(output_dir / "reported_modules.tsv", sep="\t", index=False)
    (output_dir / "reported_genes.txt").write_text("\n".join(genes) + "\n")
    return len(reported), len(genes)


def protein_sets_by_row(table: pd.DataFrame, column: str) -> list[tuple[str, ...]]:
    return [tuple(sorted(str(s).split(","))) for s in table[column]]
