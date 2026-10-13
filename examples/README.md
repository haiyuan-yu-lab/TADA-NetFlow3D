# Examples

To reproduce the manuscript results, use the
[manuscript workflow](../docs/reproducing_results.md). The examples below show
how to run the individual methods.

Run TADA for SSC/KOR/STAR and SPARK, or run NetFlow3D with a prepared real-data
subset. All commands below run from the repository root and use
[`config/example_config.yaml`](../config/example_config.yaml).

## TADA

Run both cohorts, or select one:

```bash
bash workflows/run_tada.sh config/example_config.yaml
# Alternatively, run SPARK only:
bash workflows/run_tada.sh config/example_config.yaml --cohort SPARK
```

The configuration enables de novo and inherited coding and noncoding evidence.
Each cohort is fitted separately.

| Input file | Contents |
|---|---|
| `input/tada/shared/gene_information.tsv` | Gene information |
| `input/tada/shared/priors.tsv` | Model priors |
| `input/tada/SSC_KOR_STAR/coding_counts.tsv` | SSC/KOR/STAR coding counts |
| `input/tada/SPARK/coding_counts.tsv` | SPARK coding counts |
| `input/tada/SSC_KOR_STAR/noncoding_counts.tsv` | SSC/KOR/STAR noncoding counts |
| `input/tada/SPARK/noncoding_counts.tsv` | SPARK noncoding counts |
| `input/tada/SSC_KOR_STAR/sample_sizes.tsv` | SSC/KOR/STAR sample sizes |
| `input/tada/SPARK/sample_sizes.tsv` | SPARK sample sizes |

Results are written to `results/examples/tada/<cohort>/tada_results.tsv`.

## NetFlow3D

Run with the bundled network inputs and prepared gene scores:

```bash
bash workflows/run_netflow3d.sh config/example_config.yaml
```

The network is a randomly sampled subset of HINT HQ. The configuration uses
1,000 background trials. Gene scores are supplied directly, so no BF conversion
is needed before running this example.

| Input file | Contents |
|---|---|
| `input/netflow3d/gene_scores.tsv` | Prepared gene scores (`Gene` and `Score`) |
| `input/netflow3d/interactome.tsv` | HINT HQ binary interaction subset, as UniProt pairs |
| `input/netflow3d/id_mapping.tsv` | UniProt mappings to gene and transcript IDs |
| `input/netflow3d/disrupted_interactions.tsv` | Disrupted interactions |
| `input/netflow3d/expressed_genes.txt` | Expressed transcript IDs (ENST) |

Results are written under `results/examples/netflow/output/example/`.

The example's subnetworks should not be interpreted as biological findings
from the study. Sampling removes many proteins and interactions, changing
heat diffusion and the background threshold, and the subset does
not represent the full network's interaction structure.

For your own analyses, see [input formats](../docs/input_formats.md) and
[reference data](../docs/reference_data.md).

## Output directories

`results/` contains generated runs and is ignored by Git.
`examples/expected_output/` contains saved reference results for the examples.
