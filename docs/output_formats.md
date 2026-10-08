# Output formats

## TADA

The TADA result contains `gene_id`, `Hugo_Symbol`, `LOEUF`, `pLI`, component Bayes
factors for selected evidence, total `bf`, Bayesian `fdr`, and a back-transformed
`pval`. It is a tab-separated file with a header, no row names, and no quoting
of text fields. Numeric calculations and row order follow the TADA runner.

## TADA-to-NetFlow3D conversion

The converter writes `Gene` and `Score`. By default `Score = ln(bf)`.

## NetFlow3D

- `<job>/subnetworks.txt`: comma-separated UniProt members and module size.
- `<job>/initial_state.graphml.gz`: weighted input network and node heat scores.
- `<job>/choose_delta.txt`: permutation results used to select delta.
- `<job>/final_state.graphml.gz`: directed heat-flow graph after thresholding.

Downstream biological summaries are outside this repository’s scope.
