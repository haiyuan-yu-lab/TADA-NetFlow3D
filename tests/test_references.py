#!/usr/bin/env python3
from pathlib import Path
import re
import subprocess


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = re.compile(r"/(Genomics|fs|local|home)/")


def main() -> None:
    errors = []
    # Check public repository files, including new files, while honoring Git ignores.
    names = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"], cwd=ROOT,
    ).decode().split("\0")
    for name in filter(None, names):
        path = ROOT / name
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
