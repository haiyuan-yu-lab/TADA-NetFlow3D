# TADA-NetFlow3D

Combine gene-level Bayesian association evidence with network diffusion.
TADA fits selected coding and noncoding evidence; NetFlow3D consumes prepared
gene scores, a protein interaction network, and optional disrupted edges.
Manuscript figures, enrichment analyses, and experimental validation are outside
this repository's scope.

## Workflow

```text
coding and/or noncoding counts
                    |
                    v
                  TADA
                    |
             gene BF and FDR
                    |
                    v
        TADA-to-NetFlow3D conversion
                    |
                    v
                NetFlow3D
                    |
          subnetworks + graph states
```

## Setup

```bash
conda env create -f environment.yml
conda activate tada-netflow3d
```

Python utilities can also be installed with `python -m pip install -e .`.
TADA requires R with GenomicRanges, stringr, and dplyr.

## Quick start

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
[TADA guide](docs/tada.md) for individual category switches and other options.

#### Option B: YAML configuration

Read the input paths and analysis options from the bundled
[configuration](config/example_config.yaml):

```bash
bash workflows/run_tada.sh config/example_config.yaml --cohort SPARK
```

Edit a copy of the configuration to use your own input paths and settings.
See the [TADA guide](docs/tada.md#yaml-configuration) for running both cohorts
and for configuration rules.

### 2. Convert TADA output to gene scores

```bash
python scripts/tada_to_netflow3d.py \
  --tada-results results/example/tada_results.tsv \
  --output results/example/netflow3d_gene_scores.tsv
```

The default NetFlow3D score is the natural logarithm of the total TADA Bayes
factor.

### 3. Run NetFlow3D

```bash
python src/netflow3d/run_netflow3d.py \
  --job_name example \
  --gene_score results/example/netflow3d_gene_scores.tsv \
  --binary_interactome examples/input/interactome.tsv \
  --id_mapping examples/input/id_mapping.tsv \
  --disrupted_interaction examples/input/disrupted_interactions.tsv \
  --threads 1 \
  --delta_trial 2 \
  --max_subnetwork_size 10 \
  --output_path results/example/netflow3d
```

## Example data

Example TADA inputs are provided for SSC/KOR/STAR and SPARK. They include coding
and noncoding variant counts, gene information, priors, and sample sizes. Each
run produces a gene-level results table containing Bayes factors, FDR, and
p-values. See [examples/README.md](examples/README.md) for the bundled files and
expected outputs, and [input formats](docs/input_formats.md) for preparing your
own data.

## Layout and documentation

- `src/`: scientific TADA and NetFlow3D implementations.
- `workflows/`: overall example, TADA, and NetFlow3D bash entry points.
- `config/`: shared example YAML.
- `scripts/`: YAML launcher, optional table-header checks, and BF conversion.
- `examples/`: complete prepared TADA inputs and expected results.
- `resources/`: reference-data manifest and documentation.

[TADA](docs/tada.md) · [NetFlow3D](docs/netflow3d.md) ·
[Input formats](docs/input_formats.md) · [Output formats](docs/output_formats.md) ·
[Parameters](docs/parameters.md) · [Reference data](docs/reference_data.md) ·
[Troubleshooting](docs/troubleshooting.md) · [Developer guide](docs/developer_guide.md)

## Reproducibility and testing

```bash
bash tests/run_checks.sh
```

The checks cover basic code and input-header validation. They do not execute
TADA or NetFlow3D analyses.

## Data policy

The bundled TADA examples contain real aggregate count tables. No
individual-level participant data are included. Users are responsible for
confirming that their own data may be processed and shared under the relevant
consent and governance requirements.

## Citation and license

See `CITATION.cff` and `LICENSE.md`. Publication metadata and the software license
remain subject to finalization.

## Support

Please use the GitHub issue tracker for reproducible bug reports. Include the
command, software versions, input column names, and complete error message; do
not attach controlled participant data.
