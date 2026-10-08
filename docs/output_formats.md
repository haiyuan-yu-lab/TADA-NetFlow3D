# Output formats

## TADA

The TADA result contains `gene_id`, `Hugo_Symbol`, `LOEUF`, `pLI`, component Bayes
factors for selected evidence, total `bf`, Bayesian `fdr`, and a back-transformed
`pval`. It is a tab-separated file with a header, no row names, and no quoting
of text fields. Numeric calculations and row order follow the TADA runner.

## TADA-to-NetFlow3D conversion

The converter writes `Gene` and `Score`. By default `Score = ln(bf)`.

## NetFlow3D

All files are written under `<output-dir>/<job-name>/`. Start with the subnetwork
results table; the other files support network inspection and background reuse.

| File | Purpose |
|---|---|
| `<job-name>_subnetworks.txt` | Main results: `Subnetwork_genes` lists comma-separated UniProt IDs and `Subnetwork_size` gives the number of proteins, sorted by descending size |
| `initial_state.graphml.gz` | Network before diffusion, with protein heat scores and interaction weights |
| `bg_delta.txt` | Background calibration results; reuse with `--background-delta-file` to skip new background trials |
| `final_state.graphml.gz` | Directed heat-transfer network after thresholding, for inspection or visualization |

Rerunning the same output directory and job name replaces generated files;
unrelated files are left in place. Use a different job name to retain both runs.
