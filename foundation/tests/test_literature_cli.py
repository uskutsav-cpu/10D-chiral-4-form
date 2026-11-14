"""Prevent the literature command from replacing a certified map incorrectly."""
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("degree,target,legacy", [(8, 10, False), (10, 8, False),
                                                 (8, 8, True), (10, 10, True)])
def test_incompatible_output_cannot_overwrite_canonical_map(degree, target, legacy):
    path = ROOT / f"results/order{target}_change_of_basis.json"
    original = path.read_bytes()
    command = [sys.executable, str(ROOT / "scripts/map_literature_basis.py"),
               "--degree", str(degree), "--out", str(path)]
    if legacy:
        command.append("--legacy")
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 2
    assert "cannot overwrite a canonical map" in result.stderr
    assert path.read_bytes() == original
