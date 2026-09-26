# Output formats

## TADA

The TADA result contains `gene_id`, `gene`, `LOEUF`, `pLI`, component Bayes
factors, total `bf`, Bayesian `fdr`, and a back-transformed `pval`.

## TADA-to-NetFlow3D conversion

The converter writes `Gene` and `Score`. By default `Score = ln(bf)`.

## NetFlow3D

- `<job>/subnetworks.txt`: comma-separated UniProt members and module size.
- `<job>/initial_state.graphml.gz`: weighted input network and node heat scores.
- `<job>/choose_delta.txt`: permutation results used to select delta.
- `<job>/final_state.graphml.gz`: directed heat-flow graph after thresholding.

Output details remain provisional until the scientific implementation is
finalized. Downstream biological summaries are intentionally out of scope.
