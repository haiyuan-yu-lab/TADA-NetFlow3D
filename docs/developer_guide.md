# Developer replacement guide

The current source files are provisional. Finalized scripts may replace them as
long as the repository remains self-contained and the documented interface is
preserved.

For each replacement:

1. Keep the entry-point location and existing flags where practical.
2. Remove personal absolute paths and undocumented imports.
3. Validate required columns and parameter ranges.
4. Preserve or deliberately version output schemas.
5. Update the environment and reference manifest.
6. Run `bash tests/run_checks.sh` and the full synthetic example.
7. Record user-visible changes in `CHANGELOG.md`.

Do not add figure notebooks, participant-level data, manuscript-only analyses,
or generated results to the core repository.
