# Changelog

## Unreleased

### Added

- TADA can analyze coding evidence, noncoding evidence, or both.
- Real TADA example data and expected results are included for SSC/KOR/STAR
  and SPARK.
- TADA runs can be saved in a YAML configuration and repeated for one or
  multiple datasets. Each dataset produces its own result.

### Usage changes

- TADA now requires you to select which evidence to analyze; coding evidence
  is no longer enabled automatically.
- TADA option names use hyphens, such as `--coding-file`. 
- Examples run with bash. Use `bash workflows/run_examples.sh` to run the
  bundled TADA examples.

## 0.1.0-dev — 2026-09-26

- Established the core TADA and NetFlow3D repository structure.
