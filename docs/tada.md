# TADA guide

TADA can combine de novo and inherited coding evidence for PTV, MisB, and MisA
classes. It can also include de novo and inherited regulatory-element evidence
and allocate that evidence to target genes using inverse LOEUF weighting.
Inputs are prepared count tables, mutation rates, gene information, priors,
and sample sizes.

## Direct R invocation

See [README: Option A](../README.md#option-a-explicit-arguments) for the complete
SPARK command with explicit input paths and analysis options.

See [input formats](input_formats.md) for the contents of each input table.
`--mut-col` selects the noncoding mutation-rate column; `--output-file` sets
the results path.

All evidence categories default to disabled. Group switches select all coding
classes or the corresponding enhancer component; individual class switches are
also available. For noncoding-only analysis, omit the coding file, disable coding
switches, and use the information table as the gene backbone. For coding-only
analysis, omit the noncoding file and disable its switches.

Run `Rscript src/tada/run_tada.R --help` for all options. Long options use hyphens
only. Direct R input paths are relative to the caller's directory; the helper
file is resolved relative to the runner.

## YAML configuration

Use [the example YAML](../config/example_config.yaml) to save input paths and
analysis options. Run all configured entries from the repository root with:

```bash
bash workflows/run_tada.sh config/example_config.yaml
```

Add `--cohort SPARK` to run only SPARK. Each entry produces a separate result;
the launcher does not merge inputs or combine results. `bash workflows/run_examples.sh`
runs both bundled TADA examples.

For a single dataset, keep just one entry under `tada.cohorts` and give it a
name of your choice. 

Option names match the R arguments without `--`; use `true` or `false` for
analysis switches. Settings within a cohort override shared settings under
`tada`. Give each entry its own `output-file`.

Input and output paths in YAML are relative to the repository root or absolute.
The YAML filename on the command line is relative to your working directory or
absolute. Copy the example configuration and edit its paths and settings for
your own data.

## Calculations and outputs

The selected mutation column is converted to numeric; `.` and missing rates
are excluded. Enhancer log evidence is distributed using inverse LOEUF weights.
`--in-floor` applies the existing inherited coding floor and the inherited
noncoding floor after gene allocation. The output includes gene metadata,
component and total BFs, FDR, and p-values in an unquoted TSV.

Dependencies: R with GenomicRanges, stringr, and dplyr; Python with PyYAML for
the configuration launcher. See the [example README](../examples/README.md)
for bundled files, settings, and expected outputs.
