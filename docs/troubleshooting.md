# Troubleshooting

## A file cannot be found

Run commands from the repository root or provide absolute paths. Public scripts
must not rely on laboratory-specific filesystem locations.

## Gene scores do not map to the network

Confirm that the `Gene` values match identifiers in the third column of the
UniProt mapping table. Version suffixes are removed from that mapping column;
input `Gene` values must match the resulting identifiers. Unmatched score rows
do not contribute heat. A run with no matching identifiers can complete without
reporting a mapping error, so check the mapping if the output has no heat or
subnetworks. The runner does not comprehensively validate score values; provide
finite numeric scores.

## NetFlow3D uses too much memory

Diffusion can use substantial memory within connected components. Reduce
`--n-workers` to lower the number of simultaneous background solves. The bash
wrapper limits each numerical library to one thread per worker.

## A run directory already exists

NetFlow3D writes under `<output-dir>/<job-name>/` and replaces its generated
files when that job is repeated. Use a new `--job-name` to keep separate results.
