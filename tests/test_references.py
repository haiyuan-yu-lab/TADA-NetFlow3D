#!/usr/bin/env python3
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = re.compile(r"/(Genomics|fs|local|home)/")


def main() -> None:
    errors = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or ".git" in path.parts or "source_material" in path.parts:
            continue
        if path.suffix.lower() not in {".py", ".r", ".sh", ".md", ".yaml", ".yml"}:
            continue
        text = path.read_text(errors="replace")
        if FORBIDDEN.search(text):
            errors.append(f"private absolute path: {path.relative_to(ROOT)}")
    if errors:
        raise SystemExit("\n".join(errors))


if __name__ == "__main__":
    main()
