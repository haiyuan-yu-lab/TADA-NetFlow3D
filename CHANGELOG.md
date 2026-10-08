# Changelog

## Unreleased

### Added

- NetFlow3D supports network rewiring, heat shuffling, or both to estimate
  background thresholds, and can reuse a saved background table.
- NetFlow3D method documentation explains heat diffusion and threshold rounding.

- TADA can analyze coding evidence, noncoding evidence, or both.
- Real TADA example data and expected results are included for SSC/KOR/STAR
  and SPARK.
- TADA runs can be saved in a YAML configuration and repeated for one or
  multiple datasets. Each dataset produces its own result.

### Usage changes

- NetFlow3D uses hyphenated option names and requires an output directory.
- NetFlow3D defaults to 1,000 background trials combining edge rewiring with
  uniform heat shuffling, and three-decimal threshold rounding.
- NetFlow3D results are grouped by job name; repeating a job replaces its
  generated files. The trial table is `bg_delta.txt` and subnetwork files are
  named `<job-name>_subnetworks.txt`.

- TADA now requires you to select which evidence to analyze; coding evidence
  is no longer enabled automatically.
- TADA option names use hyphens, such as `--coding-file`.
- Examples run with bash. Use `bash workflows/run_examples.sh` to run the
  bundled TADA examples.

## 0.1.0-dev — 2026-09-26

- Established the core TADA and NetFlow3D repository structure.
