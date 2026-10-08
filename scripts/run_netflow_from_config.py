#!/usr/bin/env python3
"""Run NetFlow3D using the netflow3d section of a YAML configuration."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[1]
PATH_OPTIONS = {
    "gene-score", "output-dir", "disrupted-interaction", "binary-interactome",
    "id-mapping", "expressed-genes", "background-delta-file",
}
OPTIONS = PATH_OPTIONS | {
    "job-name", "resolution", "restart-prob", "weighted-edge",
    "max-subnetwork-size", "delta-trial", "swaps-per-edge", "background-null",
    "background-weights", "heat-shuffle", "random-seed", "n-workers",
    "delta-decimals",
}


def build_command(config_file: Path) -> list[str]:
    with config_file.open() as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict) or not isinstance(config.get("netflow3d"), dict):
        raise ValueError("Configuration must contain a top-level 'netflow3d' mapping")
    unknown_sections = set(config) - {"tada", "netflow3d"}
    if unknown_sections:
        raise ValueError(f"Unknown configuration section(s): {', '.join(sorted(map(str, unknown_sections)))}")
    options = config["netflow3d"]
    unknown = set(options) - OPTIONS
    if unknown:
        raise ValueError(f"netflow3d: unknown option(s): {', '.join(sorted(map(str, unknown)))}")
    command = [sys.executable, str(ROOT / "src/netflow3d/run_netflow3d.py")]
    for key, value in options.items():
        if value is None:
            continue
        if type(value) not in (str, int, float) or (isinstance(value, str) and not value.strip()):
            raise ValueError(f"netflow3d.{key}: expected a nonempty string, number, or null")
        if key in PATH_OPTIONS:
            if not isinstance(value, str):
                raise ValueError(f"netflow3d.{key}: expected a path string or null")
            path = Path(value)
            if not path.is_absolute():
                path = ROOT / path
            value = str(path.resolve())
        command.extend(("--" + key, str(value)))
    return command


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="YAML path relative to the current directory, or absolute")
    args = parser.parse_args()
    try:
        command = build_command(args.config)
    except (OSError, ValueError, yaml.YAMLError) as error:
        parser.error(str(error))
    subprocess.run(command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
