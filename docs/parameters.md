# Parameter reference

## TADA flags

- `--dn-em`: include de novo noncoding enhancer evidence.
- `--in-em`: include inherited noncoding enhancer evidence.
- `--in-floor`: floor inherited Bayes factors at 1.

Coding PTV, MisB, and MisA categories are enabled in the provisional script.

## NetFlow3D parameters

- `--restart_prob`: fraction of heat retained during diffusion; default `0.2`.
- `--weighted_edge`: multiplier for qualifying disrupted edges; default `5.0`.
- `--resolution`: `low` selects the minimum permuted delta; `high` selects the
  median.
- `--max_subnetwork_size`: target maximum size used during delta selection.
- `--delta_trial`: number of randomized networks used to select delta.
- `--random_seed`: first seed in the deterministic sequence.
- `--threads`: worker-process count.

Defaults are provisional and must be reconciled with the finalized manuscript
analysis before release.
