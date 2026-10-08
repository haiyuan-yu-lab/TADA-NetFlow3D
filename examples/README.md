# Complete real TADA examples

These examples fit SSC/KOR/STAR and SPARK independently with full prepared
aggregate input tables. They use FuRaw coding counts and Combined Simple
noncoding counts with control-derived rare-damaging-rate scaling (`rd_control`).
They are not downsampled or combined across cohorts.

## Input files

| Example file | Contents | Data rows |
|---|---|---:|
| `input/tada/shared/gene_information.tsv` | Gene information | 18,128 |
| `input/tada/shared/priors.tsv` | Model priors | 5 |
| `input/tada/SSC_KOR_STAR/coding_counts.tsv` | SSC/KOR/STAR coding counts | 18,128 |
| `input/tada/SPARK/coding_counts.tsv` | SPARK coding counts | 18,128 |
| `input/tada/SSC_KOR_STAR/noncoding_counts.tsv` | SSC/KOR/STAR noncoding counts | 33,468 |
| `input/tada/SPARK/noncoding_counts.tsv` | SPARK noncoding counts | 33,468 |
| `input/tada/SSC_KOR_STAR/sample_sizes.tsv` | SSC/KOR/STAR sample sizes | 8 |
| `input/tada/SPARK/sample_sizes.tsv` | SPARK sample sizes | 8 |

## Run and compare

From the repository root:

```bash
bash workflows/run_examples.sh
# Or select a cohort:
bash workflows/run_tada.sh config/example_config.yaml --cohort SPARK
```

Shared settings are in `config/example_config.yaml`: all de novo and inherited
coding/enhancer evidence enabled, `mut-col: mut.em`, and `in-floor: false`.
Priors and sample sizes come from the bundled tables; `pi0` is 0.95.

Results are written to `results/examples/tada/<cohort>/tada_results.tsv`.
Corresponding expected outputs are under `examples/expected_output/tada/<cohort>/`.
Expected outputs are generated with R 4.0.5 using the supplied YAML settings.
Each output represents a separate cohort fit; it is not a combined-cohort result.

No synthetic dataset or NetFlow3D network subset is bundled. NetFlow3D examples
will be added separately; the overall example script currently runs TADA only.
