#!/usr/bin/env python3
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tada = Path(tmp) / "tada.tsv"
        output = Path(tmp) / "score.tsv"
        pd.DataFrame({"gene_id": ["A", "B", "C"], "bf": [0.1, 1.0, 10.0]}).to_csv(tada, sep="\t", index=False)
        for flags, expected in [
            ([], [0.0, 0.0, np.log(10.0)]),
            (["--transform", "log10"], [0.0, 0.0, 1.0]),
            (["--transform", "none"], [0.1, 1.0, 10.0]),
        ]:
            subprocess.run(
                ["python3", str(ROOT / "scripts/tada_to_netflow3d.py"), "--tada-results", str(tada), "--output", str(output), *flags],
                check=True,
            )
            result = pd.read_csv(output, sep="\t")
            assert result.columns.tolist() == ["Gene", "Score"]
            assert result["Gene"].tolist() == ["A", "B", "C"]
            assert np.allclose(result["Score"], expected)


if __name__ == "__main__":
    main()
