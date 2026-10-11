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

Long options use hyphens. Required arguments are job name, gene scores, and the
output parent directory. Input paths are relative to the caller's working
directory or absolute. Network and ID mapping defaults use the complete files
under `resources/netflow3d/`; supply other paths to use your own references.
Optional analysis settings can be omitted to use the
defaults; `--help` displays the runner's default seed.

| Option | Meaning / default |
|---|---|
| `-j`, `--job-name` | Required single directory name; also prefixes the subnetwork file |
| `-g`, `--gene-score` | Required `Gene`/`Score` TSV |
| `-o`, `--output-dir` | Required parent directory; results go under `<output-dir>/<job-name>/` |
| `-b`, `--binary-interactome` | Headerless PPI table; bundled HINT HQ network by default |
| `-p`, `--id-mapping` | UniProt mapping table; bundled complete mapping by default |
| `-i`, `--disrupted-interaction` | Optional `P1`/`P2` TSV; omitted by default |
| `-X`, `--expressed-genes` | Optional expressed-gene list; omitted by default |
| `-r`, `--restart-prob` | Restart probability; `0.5` |
| `-w`, `--weighted-edge` | Weight assigned to qualifying disrupted edges; `5` |
| `-R`, `--resolution` | `low` (minimum finite background delta) or `high` (median); `low` |
| `-x`, `--max-subnetwork-size` | Background calibration size cutoff, not a final module-size limit; `5` |
| `-D`, `--delta-trial` | Background trial count; `1000` |
| `-S`, `--random-seed` | Random seed for reproducible background trials; trials use consecutive seeds |
| `--swaps-per-edge` | Swaps per between-protein edge for rewire/both; `10` |
| `--background-null` | `rewire`, `heat`, or `both`; `both`. See below for details. |
| `--background-weights` | `unit`, `shuffle`, or `keep`, subject to the combinations below; `unit` |
| `--heat-shuffle` | `uniform` or `binned`, used by heat/both; `uniform`. See below for details. |
| `--background-delta-file` | Optional existing trial table; omitted by default |
| `--delta-decimals` | Nonnegative decimal places for half-up rounding with ties passing; `3` |
| `-t`, `--n-workers` | Background worker processes; `8` |
| `-h`, `--help` | Show usage |

| Background mode | Randomization | Allowed background weights |
|---|---|---|
| `rewire` | Degree-preserving edge swaps; heat remains in place | `unit`, `shuffle` |
| `heat` | Heat permutation; connections remain in place | `unit`, `keep` |
| `both` | Edge swaps followed by heat permutation | `unit`, `shuffle` |

`unit` sets background weights to one. `shuffle` redistributes between-protein
weights over rewired edges while preserving self-loop weights. `keep` leaves
weights in place. `uniform` permutes heat across all proteins; `binned` permutes
within ten groups ranked by partner count. These settings affect the background,
not the observed network's heat or disrupted-edge weights.

The bash wrapper sets OMP, MKL, OpenBLAS, and NumExpr to one thread per worker.
Only low/high resolution is exposed, and the diffusion implementation is fixed.
The same settings can be provided under `netflow3d` in YAML, using the long
option names without `--`. Omitted settings use these defaults.
