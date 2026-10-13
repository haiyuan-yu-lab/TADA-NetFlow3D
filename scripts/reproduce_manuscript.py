#!/usr/bin/env python3
"""Reproduce the manuscript TADA results and NetFlow3D modules."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import subprocess
import sys

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src/reproduction"))
from combine_cohort_bf import combine_tables
from manuscript_results import (
    MODEL_COLUMNS, assemble_tada_table,
    combine_coding_and_noncoding, export_modules,
)


def path(value: str) -> Path:
    p = Path(value)
    return p if p.is_absolute() else ROOT / p


def run_command(command: list[str], log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("w") as log:
        result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"Command failed; see {log_file}")


def run_jobs(jobs: list[tuple[list[str], Path]], workers: int, label: str) -> None:
    if workers < 1:
        raise ValueError("Worker counts must be positive")
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run_command, command, log_file) for command, log_file in jobs]
        for index, future in enumerate(as_completed(futures), 1):
            future.result()
            if index == len(jobs) or index % 50 == 0:
                print(f"{label}: {index}/{len(jobs)} runs completed", flush=True)


def tada_command(config: dict, cohort: str, output: Path) -> list[str]:
    shared = config["inputs"]
    command = ["Rscript", str(ROOT / "src/tada/run_tada.R")]
    for key in ("info-file", "prior-file"):
        command.extend(("--" + key, str(path(shared[key]))))
    command.extend(("--sample-file", str(path(shared["cohorts"][cohort]["sample-file"])),
                    "--output-file", str(output)))
    return command


def run_tada(config: dict, output_dir: Path) -> dict[str, pd.DataFrame]:
    destination = output_dir / "tada"
    jobs = []
    for cohort, inputs in config["inputs"]["cohorts"].items():
        for model in MODEL_COLUMNS:
            output = destination / "runs" / cohort / f"{model}.tsv"
            output.parent.mkdir(parents=True, exist_ok=True)
            command = tada_command(config, cohort, output)
            command.extend(("--noncoding-file", str(path(inputs["noncoding-files"][model])),
                            "--mut-col", "mut.em", "--dn-noncoding", "--in-noncoding"))
            jobs.append((command, output.with_suffix(".log")))
    print("Running TADA analyses…", flush=True)
    run_jobs(jobs, config.get("execution", {}).get("tada-workers", 4), "TADA")
    combined = {}
    cohorts = list(config["inputs"]["cohorts"])
    if cohorts != ["SSC_KOR_STAR", "SPARK"]:
        raise ValueError("The manuscript configuration requires SSC_KOR_STAR followed by SPARK")
    for model in MODEL_COLUMNS:
        tables = [pd.read_csv(destination / "runs" / c / f"{model}.tsv", sep="\t") for c in cohorts]
        combined[model] = combine_tables(*tables)
        combined[model].to_csv(destination / f"{model}_combined.tsv", sep="\t", index=False)
    fu = pd.read_csv(path(config["inputs"]["fu-coding-file"]), sep="\t")
    table = assemble_tada_table(fu, combined)
    table.to_csv(output_dir / "tada_results.tsv", sep="\t", index=False)
    print(f"TADA results saved to {output_dir / 'tada_results.tsv'}", flush=True)
    return combined


def require_reference_files(config: dict) -> None:
    for option in ("binary-interactome", "id-mapping"):
        filename = path(config["netflow3d"][option])
        if not filename.is_file():
            raise ValueError(f"Reference file missing: {filename}. Set {option} in your YAML; "
                             "see docs/reproducing_results.md. To run TADA only, use --stage tada.")


def run_netflow(config: dict, output_dir: Path, combined: dict[str, pd.DataFrame]) -> None:
    coding = pd.read_csv(path(config["inputs"]["netflow-coding-bf-file"]), sep="\t")
    integrated = combine_coding_and_noncoding(coding, combined["combined"])
    scores = integrated[["gene_id", "heat_score"]].rename(columns={"gene_id": "Gene", "heat_score": "Score"})
    expected_dir = path(config["expected-output-dir"])
    scores.to_csv(output_dir / "netflow3d_gene_scores.tsv", sep="\t", index=False)
    print(f"Gene scores saved to {output_dir / 'netflow3d_gene_scores.tsv'}", flush=True)
    command = [sys.executable, "-u", str(ROOT / "src/netflow3d/run_netflow3d.py"),
               "--gene-score", str(output_dir / "netflow3d_gene_scores.tsv"),
               "--output-dir", str(output_dir / "netflow3d")]
    path_options = {"binary-interactome", "id-mapping", "expressed-genes",
                    "disrupted-interaction", "background-delta-file"}
    options = config["netflow3d"]
    for key, value in options.items():
        if value is not None:
            command.extend(("--" + key, str(path(value) if key in path_options else value)))
    if options.get("background-delta-file") is not None:
        print("Running NetFlow3D with saved background trials…", flush=True)
    else:
        print("Running NetFlow3D; background calibration can take several hours…", flush=True)
    run_command(command, output_dir / "netflow3d.log")
    subnetworks = output_dir / "netflow3d" / options["job-name"] / f"{options['job-name']}_subnetworks.txt"
    export_modules(subnetworks, path(options["id-mapping"]), output_dir, expected_dir)
    print(f"Modules saved to {output_dir / 'reported_modules.tsv'}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--stage", choices=("all", "tada", "netflow"), default="all",
                        help="Run the complete workflow, TADA only, or integration and NetFlow3D after TADA")
    args = parser.parse_args()
    try:
        with args.config.open() as handle:
            config = yaml.safe_load(handle)
        output_dir = path(config["output-dir"])
        if args.stage in {"all", "netflow"}:
            require_reference_files(config)
        output_dir.mkdir(parents=True, exist_ok=True)
        if args.stage in {"all", "tada"}:
            combined = run_tada(config, output_dir)
        else:
            combined = {model: pd.read_csv(output_dir / "tada" / f"{model}_combined.tsv", sep="\t")
                        for model in MODEL_COLUMNS}
        if args.stage in {"all", "netflow"}:
            run_netflow(config, output_dir, combined)
        print(f"Completed. Results: {output_dir}", flush=True)
    except (OSError, ValueError, KeyError, RuntimeError, AssertionError) as error:
        parser.exit(1, f"Reproduction stopped: {error}\n")


if __name__ == "__main__":
    main()
