# TADA-NetFlow3D

TADA-NetFlow3D is an integrative computational framework for systematically characterizing the effects of rare variation, from coding and noncoding variant classes to individual genes and biological modules.

The framework combines two complementary components. Extended TADA integrates rare variants across coding and noncoding classes to quantify gene-level association, enabling noncoding variant burden observed at the element level to be translated into evidence for specific target genes.

NetFlow3D extends this gene-level evidence to the systems level by propagating association signals through a protein–protein interaction network informed by three-dimensional protein structural information. It identifies significantly interconnected gene modules in which genetic association converges on shared biological processes and pathways.

Together, TADA-NetFlow3D provides a multiscale framework connecting rare genetic variation → gene-level association → molecular networks and biological processes, allowing coding and noncoding genetic evidence to be jointly analyzed while preserving the distinct information contributed by each.

## Workflow

```text
coding counts + optional noncoding counts
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

This repository ends at core TADA and NetFlow3D results. Manuscript figures,
enrichment analyses, experimental validation, and other downstream analyses are
intentionally excluded.

## Installation

```bash
conda env create -f environment.yml
conda activate tada-netflow3d
```

Python-only installation is also available:

```bash
python -m pip install -e .
```

## Quick start

Run the bundled synthetic example from the repository root:

```bash
bash workflows/run_example.sh
```

Results are written under `results/example/`. The example uses a deliberately
small synthetic network and must not be interpreted biologically.

### 1. Run TADA

```bash
Rscript src/tada/run_tada.R \
  --coding-file examples/input/coding_counts.tsv \
  --noncoding-file examples/input/noncoding_counts.tsv \
  --info-file examples/input/gene_information.tsv \
  --prior-file examples/input/priors.tsv \
  --sample-file examples/input/sample_sizes.tsv \
  --dn-em --in-em \
  --output-file results/example/tada_results.tsv
```

Coding PTV, MisB, and MisA categories are enabled by the current provisional
implementation. `--dn-em` and `--in-em` enable de novo and inherited noncoding
enhancer evidence. `--in-floor` optionally floors inherited Bayes factors at 1.

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

## Documentation

- [TADA guide](docs/tada.md)
- [NetFlow3D guide](docs/netflow3d.md)
- [Input formats](docs/input_formats.md)
- [Output formats](docs/output_formats.md)
- [Parameters](docs/parameters.md)
- [Reference data](docs/reference_data.md)
- [Troubleshooting](docs/troubleshooting.md)
- [Developer replacement guide](docs/developer_guide.md)

## Reproducibility and testing

```bash
bash tests/run_checks.sh
```

The checks validate repository references and input schemas and run the small
Python NetFlow3D example when dependencies are available. The full example also
requires R.

## Data policy

No individual-level participant data are included. All bundled example records
are synthetic. Users are responsible for confirming that their own data may be
processed and shared under the relevant consent and governance requirements.

## Citation and license

Citation metadata are recorded in `CITATION.cff`. Publication-specific fields
and the final software license must be reviewed before a public release. Until a
license is selected, no reuse permission is granted beyond applicable law.

## Support

Please use the GitHub issue tracker for reproducible bug reports. Include the
command, software versions, input column names, and complete error message; do
not attach controlled participant data.
