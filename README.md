# TADA-NetFlow3D

Combine gene-level Bayesian association evidence with network diffusion.
TADA fits selected coding and noncoding evidence; NetFlow3D consumes prepared
gene scores, a protein interaction network, and optional disrupted edges.
Manuscript figures, enrichment analyses, and experimental validation are outside
this repository's scope.

## Reproduce the manuscript results

### 1. Set up the environment

```bash
conda env create -f environment.yml
conda activate tada-netflow3d
```

Run the following commands from the repository root.

### 2. Run the workflow

Input and output paths are provided in
[config/manuscript_config.yaml](config/manuscript_config.yaml).

```bash
bash workflows/reproduce_manuscript.sh config/manuscript_config.yaml
```

This command runs the study's gene association and network analyses and saves
the gene results and protein interaction modules under `results/manuscript/`.

### 3. Inspect the results

| Output under `results/manuscript/` | Contents |
|---|---|
| `tada_results.tsv` | Gene association results |
| `netflow3d_gene_scores.tsv` | Integrated scores used by NetFlow3D |
| `reported_modules.tsv` | The 13 reported modules, with protein and gene identifiers |
| `netflow3d/manuscript/` | All subnetworks, network states and background results |

The workflow uses prepared coding Bayes factors and runs network background
trials, which can take several hours. See the [reproduction guide](docs/reproducing_results.md)
for staged runs, saved intermediate results and runtime information.

To reproduce Table 8 alone, using only bundled inputs:

```bash
bash workflows/reproduce_manuscript.sh config/manuscript_config.yaml --stage tada
```

## Examples and custom analyses

For a smaller network example and the individual TADA and NetFlow3D commands,
see [examples](examples/README.md) and [analyses with your own data](docs/custom_data.md).
These examples introduce the interfaces; the manuscript workflow above uses
the complete study inputs.

[TADA](docs/tada.md) · [NetFlow3D](docs/netflow3d.md) ·
[Input formats](docs/input_formats.md) · [Output formats](docs/output_formats.md) ·
[Parameters](docs/parameters.md) · [Reference data](docs/reference_data.md)

## Repository layout

- `reproduction/`: study inputs and reference results for the manuscript workflow.
- `src/`: TADA, NetFlow3D and manuscript evidence-integration code.
- `workflows/`: bash entry points for reproduction and individual methods.
- `config/`: manuscript and example configurations.
- `examples/`: individual-method inputs and expected results.
- `scripts/`: configuration launchers and table utilities.
- `resources/`: complete HINT HQ interactome, UniProt mapping and their documentation.

## Development

```bash
bash tests/run_checks.sh
```

These checks cover code syntax, repository references and input headers. They
run converter checks but do not execute the full scientific workflow.
See the [developer guide](docs/developer_guide.md).

## Data, citation and support

Bundled study tables contain aggregate counts, prepared inputs and results;
no participant identifiers are included. The complete HINT HQ interactome and
UniProt mapping are provided under `resources/netflow3d/`.

See `CITATION.cff` and `LICENSE.md`. Publication metadata and the software license
remain subject to finalization. For bug reports, use the GitHub issue tracker
and include the command and error message.
