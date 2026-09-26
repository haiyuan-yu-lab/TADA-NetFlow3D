#!/usr/bin/env python3
"""Reference-data downloader placeholder for the provisional scaffold."""

from pathlib import Path
import csv


def main() -> None:
    manifest = Path(__file__).resolve().parents[1] / "resources" / "manifest.tsv"
    with manifest.open() as handle:
        pending = [row["resource"] for row in csv.DictReader(handle, delimiter="\t") if row["url"] == "PENDING"]
    if pending:
        raise SystemExit(
            "Approved download locations are still pending for: " + ", ".join(pending)
            + ". See docs/reference_data.md."
        )


if __name__ == "__main__":
    main()
