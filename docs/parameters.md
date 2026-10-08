# Parameter reference

## TADA flags

All long options use hyphens, without underscore aliases.

| Options | Meaning / default |
|---|---|
| `--info-file`, `--prior-file`, `--sample-file`, `--output-file` | Required paths |
| `--coding-file` | Required if any coding component is selected |
| `--noncoding-file`, `--mut-col` | Required if any noncoding component is selected |
| `--dn-coding`, `--in-coding` | Select all three coding classes for that component |
| `--dn-noncoding`, `--in-noncoding` | Select the corresponding enhancer component |
| `--dn-ptv`, `--in-ptv` | Select individual PTV components |
| `--dn-misb`, `--in-misb` | Select individual MisB components |
| `--dn-misa`, `--in-misa` | Select individual MisA components |
| `--dn-em`, `--in-em` | Select individual enhancer components |
| `--in-floor` | Optional inherited BF flooring; default off |
| `--help`, `-h` | Show usage |

Every evidence switch defaults to off. No selected evidence is an error.
Group switches retain the local implementation's behavior; flags have no value
argument. Gene metadata and statistical parameters come from the input tables.

## NetFlow3D parameters

- `--restart_prob`: fraction of heat retained during diffusion; default `0.2`.
- `--weighted_edge`: multiplier for qualifying disrupted edges; default `5.0`.
- `--resolution`: `low` selects the minimum permuted delta; `high` selects the
  median.
- `--max_subnetwork_size`: target maximum size used during delta selection.
- `--delta_trial`: number of randomized networks used to select delta.
- `--random_seed`: first seed in the deterministic sequence.
- `--threads`: worker-process count.

These are the existing NetFlow3D runner defaults. Set explicit parameters for
your analysis. They are not automatically inferred from the TADA YAML.
