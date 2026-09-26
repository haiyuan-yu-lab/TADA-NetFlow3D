# TADA guide

`src/tada/run_tada.R` is the supported entry point. The present implementation
combines coding de novo and inherited evidence for PTV, MisB, and MisA classes.
It can optionally add de novo and inherited regulatory-element evidence and
allocate that evidence to target genes using inverse LOEUF weighting.

Run `Rscript src/tada/run_tada.R --help` for the current interface. Input column
definitions are in `input_formats.md`.

Scientific assumptions, priors, mutation-rate derivation, count-table creation,
and cohort QC must be documented by the analyst. Raw-variant preprocessing is
not part of this repository.

The implementation is provisional. A finalized replacement should preserve the
documented flags and output columns or include an explicit migration note.
