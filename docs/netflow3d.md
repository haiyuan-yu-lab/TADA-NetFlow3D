# NetFlow3D guide

NetFlow3D places TADA-derived heat on an optionally expression-filtered protein
network, upweights selected disrupted interactions, performs network diffusion,
selects a delta threshold using randomized backgrounds, and reports strongly
connected hot subnetworks. Prepared gene scores from other sources can also be
used.

## Running an analysis

See [README: Run NetFlow3D](../README.md#3-run-netflow3d) for a complete command.
The bash entry point runs in your active environment:

```bash
bash workflows/run_netflow3d.sh --help
```

Supply a prepared `Gene`/`Score` table with `--gene-score`, a protein network
with `--binary-interactome`, and a UniProt mapping with `--id-mapping`.
Use `--job-name` to name the run and `--output-dir` to choose its output parent
directory. Input paths are relative to your working directory or absolute.
See [input formats](input_formats.md) for the table columns.

## Optional arguments

- `--disrupted-interaction` supplies interactions to upweight. Omit it to use
  the network without disrupted-edge weighting.
- `--expressed-genes` restricts the network to proteins mapped from an expression
  list. Omit it to use the network without expression filtering.
- `--n-workers` controls background worker processes. Reduce it when fewer CPUs
  or less memory are available.

The default analysis settings can be used without listing them in the command.
See the [parameter reference](parameters.md#netflow3d-parameters) for defaults,
background modes, and optional overrides. All long options use hyphens.

## YAML configuration

Run the bundled subset with the shared configuration:

```bash
bash workflows/run_netflow3d.sh config/example_config.yaml
```

The `netflow3d` section uses the same option names without `--`. Omitted options
use the runner defaults; optional file paths may be `null`. Paths inside YAML
are relative to the repository root or absolute. The configuration filename is
relative to your current directory or absolute. The TADA section is ignored
when running NetFlow3D. Edit a copy to use your own data.

## Reusing backgrounds

Pass `--background-delta-file` to reuse a saved background trial table instead
of calculating new trials. Use a table generated with the same inputs and
calibration settings. The runner checks recorded background settings but does
not verify that the input datasets match.

## Saving results

Results go into `<output-dir>/<job-name>/`. Repeating that combination replaces
its generated files. Change the job name to keep another run. Unrelated files
in the directory are preserved. See [output formats](output_formats.md#netflow3d)
for the generated files.
