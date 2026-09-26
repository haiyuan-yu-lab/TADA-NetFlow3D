# Reference data

Large third-party files are intentionally not stored in Git. The production
workflow is expected to use a versioned HINT binary interactome, UniProt ID
mapping, and optionally PIONEER protein-interface annotations.

Before public release, `resources/manifest.tsv` must contain an approved source
URL, exact version, SHA-256 checksum, citation, and license for every resource.
The downloader should refuse files whose checksums do not match.

Do not copy reference files from private cluster paths into a public repository
without checking redistribution terms.
