# Analyses with your own data

TADA and NetFlow3D can be run separately or connected using prepared gene scores.
For the manuscript results, use the [reproduction workflow](reproducing_results.md).

## Run TADA and NetFlow3D

Run the following commands from the repository root.

### 1. Run TADA

Choose either Option A or Option B to run the SPARK example with the same
inputs and analysis settings.

#### Option A: Explicit arguments

```bash
mkdir -p results/examples/tada/SPARK

Rscript src/tada/run_tada.R \
  --coding-file examples/input/tada/SPARK/coding_counts.tsv \
  --noncoding-file examples/input/tada/SPARK/noncoding_counts.tsv \
  --info-file examples/input/tada/shared/gene_information.tsv \
  --prior-file examples/input/tada/shared/priors.tsv \
  --sample-file examples/input/tada/SPARK/sample_sizes.tsv \
  --dn-coding \
  --in-coding \
  --dn-noncoding \
  --in-noncoding \
  --mut-col mut.em \
  --output-file results/examples/tada/SPARK/tada_results.tsv
```

| Input argument | Contents |
|---|---|
| `--coding-file` | Gene-level coding variant counts and mutation rates |
| `--noncoding-file` | Noncoding variant counts and mutation rates |
| `--info-file` | Gene information, including LOEUF weights |
| `--prior-file` | Model priors |
| `--sample-file` | Sample sizes for each evidence category |

The four analysis switches enable de novo (`dn`) and inherited (`in`) coding
and noncoding evidence. All evidence categories are disabled by default.
`--mut-col mut.em` selects the noncoding mutation-rate column. For coding-only
or noncoding-only analysis, supply the corresponding input and enable only
its switches; keep the information, prior, and sample-size tables.

Results are written to the file specified by `--output-file`. See the
[TADA guide](tada.md) for individual category switches and other options.

#### Option B: YAML configuration

Read the input paths and analysis options from the bundled
[configuration](../config/example_config.yaml):

```bash
bash workflows/run_tada.sh config/example_config.yaml --cohort SPARK
```

Edit a copy of the configuration to use your own input paths and settings.
See the [TADA guide](tada.md#yaml-configuration) for running both cohorts
and for configuration rules.

### 2. Prepare gene scores

The NetFlow3D example uses the prepared real-data subset in
`examples/input/netflow3d/gene_scores.tsv`; no conversion is needed to run it.
For your own TADA results, the optional converter writes a `Gene`/`Score` table:

```bash
python scripts/tada_to_netflow3d.py \
  --tada-results results/examples/tada/SPARK/tada_results.tsv \
  --output results/examples/netflow/input/netflow3d_gene_scores.tsv
```

The converter defaults to `Score = max(0, ln(BF))`, using the total TADA Bayes
factor. You can also supply your own prepared `Gene`/`Score` table.

### 3. Run NetFlow3D

Run the bundled real-data subset with explicit input paths. Disrupted
interactions and an expressed-gene list are optional for your own analyses.

```bash
bash workflows/run_netflow3d.sh \
  --job-name example \
  --gene-score examples/input/netflow3d/gene_scores.tsv \
  --binary-interactome examples/input/netflow3d/interactome.tsv \
  --id-mapping examples/input/netflow3d/id_mapping.tsv \
  --disrupted-interaction examples/input/netflow3d/disrupted_interactions.tsv \
  --expressed-genes examples/input/netflow3d/expressed_genes.txt \
  --output-dir results/examples/netflow/output
```

The command uses the default analysis settings, including 1,000 background
trials. You can run the same example from the shared configuration instead:

```bash
bash workflows/run_netflow3d.sh config/example_config.yaml
```

Use `bash workflows/run_examples.sh` to run both TADA cohorts and the prepared
NetFlow3D subset. The examples use their own bundled inputs; the NetFlow3D
subset does not consume the TADA example outputs. Optional settings are listed
in the [parameter reference](parameters.md#netflow3d-parameters).

Results are written to `results/examples/netflow/output/example/`. Repeating
the same job name and output directory replaces that job's generated files; a different job name
keeps the runs separate. See the [NetFlow3D guide](netflow3d.md) for the
method.


See [input formats](input_formats.md) for your own tables and
[reference data](reference_data.md) for full network inputs.
