# Reference data

NetFlow3D uses two reference files, both provided under `resources/netflow3d/`:

| Argument | Bundled file | Contents |
|---|---|---|
| `--binary-interactome` | `HomoSapiens_binary_HINThq.txt` | Complete HINT HQ binary protein interaction network; headerless UniProt pairs |
| `--id-mapping` | `HUMAN_9606_idmapping.dat.gz` | Complete human UniProt identifier mapping; accession, identifier type and identifier |

The [manuscript configuration](../config/manuscript_config.yaml) uses these files.
The NetFlow3D runner also uses them by default. The mapping remains gzip-compressed
and can be read directly without unpacking it.

For your own references, pass other paths through these arguments or your YAML.
See [input formats](input_formats.md) for the expected columns. Changing the
reference network or mapping can change the resulting subnetworks.

The [small NetFlow3D example](../examples/README.md) explicitly uses its own
subset network and mapping instead of the complete references.
