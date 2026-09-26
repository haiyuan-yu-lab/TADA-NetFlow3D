# NetFlow3D guide

`src/netflow3d/run_netflow3d.py` places TADA-derived heat on proteins, optionally
filters the network by expression, upweights selected disrupted interactions,
performs network diffusion, selects a delta threshold by permutation, and
reports connected hot subnetworks.

Production runs require approved full reference resources. The bundled network
is synthetic and intended only for interface testing.

The current implementation calculates dense diffusion matrices independently
for connected components and may require substantial memory on large networks.
Record interactome version, mapping version, parameters, random seed, and
software version for every analysis.
