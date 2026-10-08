#!/usr/bin/env python3
"""Run independent TADA cohorts from one YAML configuration."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess

import yaml

ROOT = Path(__file__).resolve().parents[1]
PATH_OPTIONS = (
    "coding-file", "noncoding-file", "info-file", "prior-file",
    "sample-file", "output-file",
)
VALUE_OPTIONS = (*PATH_OPTIONS, "mut-col")
FLAG_OPTIONS = (
    "dn-coding", "in-coding", "dn-noncoding", "in-noncoding",
    "dn-ptv", "in-ptv", "dn-misb", "in-misb", "dn-misa", "in-misa",
    "dn-em", "in-em", "in-floor",
)
OPTIONS = set(VALUE_OPTIONS + FLAG_OPTIONS)


def check_options(options: dict, label: str) -> None:
    unknown = set(options) - OPTIONS
    if unknown:
        raise ValueError(f"{label}: unknown option(s): {', '.join(sorted(map(str, unknown)))}")
    for key, value in options.items():
        if key in FLAG_OPTIONS:
            if type(value) is not bool:
                raise ValueError(f"{label}.{key}: expected a YAML boolean (true or false)")
        elif value is not None and (not isinstance(value, str) or not value.strip()):
            raise ValueError(f"{label}.{key}: expected a nonempty string or null")


def build_commands(config_file: Path, cohort: str | None = None) -> list[tuple[str, list[str], Path]]:
    with config_file.open() as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict) or "tada" not in config:
        raise ValueError("Configuration must contain a top-level 'tada' mapping")
    unknown_sections = set(config) - {"tada", "netflow3d"}
    if unknown_sections:
        raise ValueError(f"Unknown configuration section(s): {', '.join(sorted(map(str, unknown_sections)))}")
    # A shared configuration may include NetFlow3D settings; this launcher uses TADA only.
    tada = config["tada"]
    if not isinstance(tada, dict):
        raise ValueError("tada must be a mapping")
    cohorts = tada.get("cohorts")
    if not isinstance(cohorts, dict) or not cohorts:
        raise ValueError("tada.cohorts must be a nonempty mapping")
    shared = {key: value for key, value in tada.items() if key != "cohorts"}
    check_options(shared, "tada")
    if cohort is not None and cohort not in cohorts:
        raise ValueError(f"Unknown cohort {cohort!r}; available: {', '.join(map(str, cohorts))}")
    commands = []
    outputs = set()
    for name, overrides in cohorts.items():
        if not isinstance(name, str) or not name:
            raise ValueError("Cohort names must be nonempty strings")
        if not isinstance(overrides, dict):
            raise ValueError(f"tada.cohorts.{name} must be a mapping")
        check_options(overrides, f"tada.cohorts.{name}")
        if cohort is not None and name != cohort:
            continue
        options = {**shared, **overrides}
        for key in PATH_OPTIONS:
            value = options.get(key)
            if value is None:
                continue
            path = Path(value)
            if not path.is_absolute():
                path = ROOT / path
            path = path.resolve()
            options[key] = str(path)
        if options.get("output-file") is None:
            raise ValueError(f"{name}: output-file is required to prepare the output directory")
        output = Path(options["output-file"])
        if output in outputs:
            raise ValueError(f"Multiple cohorts use the same output-file: {output}")
        outputs.add(output)
        args = ["Rscript", str(ROOT / "src/tada/run_tada.R")]
        for key in VALUE_OPTIONS:
            if options.get(key) is not None:
                args.extend(("--" + key, options[key]))
        args.extend("--" + key for key in FLAG_OPTIONS if options.get(key, False))
        commands.append((name, args, output))
    return commands


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="YAML path relative to the current directory, or absolute")
    parser.add_argument("--cohort", help="Run one configured cohort; default: all, in YAML order")
    args = parser.parse_args()
    try:
        commands = build_commands(args.config, args.cohort)
    except (OSError, ValueError, yaml.YAMLError) as error:
        parser.error(str(error))
    for name, command, output in commands:
        output.parent.mkdir(parents=True, exist_ok=True)
        print(f"Running TADA: {name} -> {output}", flush=True)
        subprocess.run(command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
