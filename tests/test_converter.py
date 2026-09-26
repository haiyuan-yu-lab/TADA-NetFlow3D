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
        pd.DataFrame({"gene_id": ["A", "B"], "bf": [1.0, np.e]}).to_csv(tada, sep="\t", index=False)
        subprocess.run(
            ["python3", str(ROOT / "scripts/tada_to_netflow3d.py"), "--tada-results", str(tada), "--output", str(output)],
            check=True,
        )
        result = pd.read_csv(output, sep="\t")
        assert result.columns.tolist() == ["Gene", "Score"]
        assert np.allclose(result["Score"], [0.0, 1.0])


if __name__ == "__main__":
    main()
