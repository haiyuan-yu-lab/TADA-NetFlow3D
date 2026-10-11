# Reproduce the manuscript results

Run the study's gene association and network analyses using the bundled inputs.

## Run the workflow

Set up and activate the environment described in the [README](../README.md),
then run from the repository root:

```bash
bash workflows/reproduce_manuscript.sh config/manuscript_config.yaml
```

The configuration already points to the study inputs and complete network
references. The workflow runs TADA, combines the prepared coding evidence with
noncoding evidence, and runs NetFlow3D. Progress appears in the terminal;
individual run logs are saved with the results.

## Outputs

Files are saved under `results/manuscript/`:

| File or directory | Contents |
|---|---|
| `tada_results.tsv` | Gene association results, FDR and p-values |
| `netflow3d_gene_scores.tsv` | Integrated gene scores used by NetFlow3D |
| `reported_modules.tsv` | Protein and gene memberships for the reported modules |
| `tada/` | TADA results and logs |
| `netflow3d/manuscript/` | All subnetworks, graph states and background results |
| `netflow3d.log` | NetFlow3D progress and diagnostics |

Repeating a run replaces its outputs. Change `output-dir` in the YAML to keep
a separate run.

## Run in stages

To run the gene association analysis alone:

```bash
bash workflows/reproduce_manuscript.sh config/manuscript_config.yaml --stage tada
```

Afterward, prepare the gene scores and run NetFlow3D:

```bash
bash workflows/reproduce_manuscript.sh config/manuscript_config.yaml --stage netflow
```

Use the same configuration and output directory for both stages.

## Reuse network backgrounds

To reuse background trials, set `netflow3d.background-delta-file` in the YAML to
an existing `bg_delta.txt` from the same inputs and settings. The default `null`
calculates new trials.

## Runtime

Allow several hours for the complete workflow. The original network job took
approximately 1 hour 45 minutes with 16 workers and a recorded peak memory usage
of 38.3 GiB. The current default is 8 workers; runtime and memory usage depend
on the available computing resources and worker count. For worker settings and
other options, see the
[parameter reference](parameters.md#netflow3d-parameters).

For input descriptions, see [the study input inventory](../reproduction/README.md)
and [reference data](reference_data.md). For analyses with other data, see
[the custom-data guide](custom_data.md).
