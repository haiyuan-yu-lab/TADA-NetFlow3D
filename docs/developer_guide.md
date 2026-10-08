# Developer guide

Scientific implementations live in `src/`; launchers and format utilities live
in `scripts/`; bash workflows live in `workflows/`. TADA YAML keys follow the
runner's long-option names. Keep configuration, help, and documentation aligned.

For code changes:

1. Keep scientific changes separate from interface or formatting changes.
2. Document required inputs, selected evidence, output schemas, and dependencies.
3. Keep output formats stable. Document any incompatible changes in `CHANGELOG.md`.
4. Run `bash tests/run_checks.sh` and relevant scientific comparisons.

Do not add manuscript notebooks or participant-level data to the core software.
