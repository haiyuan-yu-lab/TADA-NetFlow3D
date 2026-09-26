# Troubleshooting

## A file cannot be found

Run commands from the repository root or provide absolute paths. Public scripts
must not rely on laboratory-specific filesystem locations.

## Gene scores do not map to the network

Confirm the gene namespace and mapping release. Ensembl version suffixes are
stripped, but incompatible IDs cannot be inferred automatically.

## NetFlow3D uses too much memory

The provisional implementation uses dense matrices within connected components.
Use an appropriate compute node and record memory requirements. Optimization is
a planned scientific-code polishing task.

## A run directory already exists

Choose a new job name or move the previous result. The program refuses to mix
new output with an existing nonempty run.
