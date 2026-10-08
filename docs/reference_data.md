# Reference data

Large third-party files are intentionally not stored in Git. The production
workflow is expected to use a versioned HINT binary interactome, UniProt ID
mapping, and optionally PIONEER protein-interface annotations.

Before public release, `resources/manifest.tsv` must contain an approved source
URL, exact version, SHA-256 checksum, citation, and license for every resource.
Download resources from their documented sources and compare their checksums
with the manifest. No automatic downloader is provided.

Do not copy reference files from private cluster paths into a public repository
without checking redistribution terms.

The bundled [NetFlow3D example](../examples/README.md) uses a small subset and
does not require full reference downloads. For your own analysis, supply the
external files through `--binary-interactome` and `--id-mapping`, or set their
paths in your YAML configuration.
