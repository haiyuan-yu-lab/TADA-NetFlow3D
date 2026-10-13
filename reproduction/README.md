# Manuscript inputs and reference results

Use the [manuscript workflow](../docs/reproducing_results.md) to reproduce
Supplementary Tables 8 and 10. The configuration connects the required
inputs; no manual assembly of tables is needed.

## Inputs

| Location | Contents |
|---|---|
| `../examples/input/tada/` | Combined noncoding counts, gene information, priors and sample sizes |
| `input/tada/SSC_KOR_STAR/` and `input/tada/SPARK/` | Noncoding counts for Table 8 cell and tissue analyses |
| `input/tada/fu_coding_bayes_factors.tsv` | Published coding evidence used in Table 8 |
| `input/netflow3d/coding_bayes_factors.tsv` | Prepared coding evidence for NetFlow3D score integration |
| `input/netflow3d/expressed_genes.txt` | Fetal frontal cortex expression list for network filtering |
| `input/netflow3d/disrupted_interactions.tsv` | Normalized protein pairs from Supplementary Table 9 |

The two coding evidence tables serve different analyses: the Fu table is used
for Table 8, while the prepared NetFlow3D table supplies one combined coding
value per gene in `log_bf.coding` (natural-log BF). The workflow combines this
value with noncoding evidence to prepare the network gene scores.

The noncoding mutation-rate column is `mut.em`. The tables contain aggregate
counts and prepared identifiers, without participant identifiers. The complete
HINT HQ interactome and UniProt mapping are provided in `../resources/netflow3d/`;
their paths are set in the manuscript configuration.

## Reference results

`expected_output/` contains Table 8, integrated NetFlow3D scores, all Table 10
subnetworks and Table 11 module memberships. Generated results are saved under
`results/manuscript/`.
