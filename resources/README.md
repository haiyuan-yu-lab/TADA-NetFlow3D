# Reference resources

Production NetFlow3D runs require a binary protein interactome and a UniProt ID
mapping table. These resources are not committed here because they are large and
have independent versioning and licensing requirements.

The synthetic files in `examples/input/` are sufficient for testing. Before a
public release, maintainers should add approved download URLs, release versions,
checksums, citations, and licenses to `manifest.tsv` and implement those entries
in `scripts/download_reference_data.py`.
