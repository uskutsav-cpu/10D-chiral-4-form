#!/usr/bin/env python3
"""Build a deterministic, graph-only reproducibility archive and file manifest."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]

from sdinv.catalog import atomic_write_json
from verify_rank81 import verify


def build():
    basis = json.loads((ROOT / "results/rank81_basis.json").read_text())
    certificate = json.loads((ROOT / "results/rank81_certificate.json").read_text())
    verify(basis, certificate)
    for path in ("results/rank81_lorentz.json", "results/order8_change_of_basis.json",
                 "results/order10_change_of_basis.json", "results/classification_validation.json",
                 "paper/manuscript.pdf"):
        if not (ROOT / path).is_file():
            raise ValueError(f"missing release deliverable {path}")
    validation = json.loads((ROOT / "results/classification_validation.json").read_text())
    if not validation["passed"]:
        raise ValueError("release validation did not pass")
    # This bundle's pytest entry point is deliberately its portable delivery
    # suite. The full historical tensor suite remains in the Git repository.
    files = list((ROOT / "src/sdinv").glob("*.py"))
    files += [ROOT / "scripts" / f"{name}.py" for name in (
        "run_6d", "run_10d", "run_degree", "degree10_pipeline", "degree12_pipeline",
        "graph_to_latex", "search_rank81", "verify_rank81", "map_literature_basis",
        "validate_rank81_lorentz", "verify_independent_audit", "build_classification_release")]
    files += [ROOT / p for p in (
        "INDEPENDENT_AUDIT.md", "tests/test_roadmap.py", "tests/test_literature_cli.py",
        "tests/test_published_degree10_source_reading.py", "requirements.txt", "requirements-lock.txt", "pytest.ini",
        "docs/CLASSIFICATION_ROADMAP.md", "docs/degree10.md", "docs/degree12.md",
        "docs/PUBLISHED_DEGREE10_INDEX_AUDIT.md",
        "results/10d_order8.json", "results/10d_order10.json", "results/10d_order12.json",
        "results/10d_graph_catalog.json", "results/degree10/B10_coordinates_per_prime.json",
        "results/rank81_basis.json", "results/rank81_certificate.json", "results/rank81_lorentz.json",
        "results/order8_change_of_basis.json", "results/order10_change_of_basis.json",
        "results/classification_validation.json", "paper/manuscript.tex", "paper/manuscript.pdf",
        "paper/tables/certificate_cells.tex", "paper/tables/rank81_formulas.tex")]
    files += [p for p in (ROOT / "results/audit").rglob("*")
              if p.is_file() and p.suffix in {".py", ".json", ".md"}
              and "__pycache__" not in p.parts]
    files = sorted(set(files))
    entries = [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size,
                "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
    directory = ROOT / "release/classification"
    directory.mkdir(parents=True, exist_ok=True)
    manifest = {"schema": 1, "files": entries, "scope": "Portable graph-only certificate, delivery tests, formulas and literature-map evidence. Full historical test suite is available in the Git repository."}
    atomic_write_json(directory / "manifest.json", manifest)
    readme = """# Exact graph classification package

Create a Python environment, install requirements.txt, then run:

    python -m pytest tests -o addopts='' -q
    python scripts/verify_rank81.py
    python scripts/verify_rank81.py --recompute
    python scripts/verify_independent_audit.py
    python scripts/map_literature_basis.py --degree 8
    python scripts/map_literature_basis.py --degree 10
    python scripts/validate_rank81_lorentz.py

The first verification command checks saved matrices; --recompute repeats
every graph evaluation. All 81 functions and all 126 coordinates are explicit.
See docs/CLASSIFICATION_ROADMAP.md for proof boundaries and literature maps.
manifest.json records SHA-256 hashes of every payload file.
The independent checker verifies saved minors, exact CRT evaluations and
rational maps without recomputing the tensor contractions. Source polynomial
identities rely on the explicitly cited Hilbert upper bounds. The frozen
literature transcription failed and was repaired; its historical matrices
remain under results/audit/frozen-literature/.

This portable delivery suite is distinct from the full historical tensor
suite in the repository. No external spinor archive is required.
"""
    (directory / "README.md").write_text(readme)
    archive = directory / "rank81-graph-classification.zip"
    payload = {entry["path"]: (ROOT / entry["path"]).read_bytes() for entry in entries}
    payload["manifest.json"] = (directory / "manifest.json").read_bytes()
    payload["README.md"] = readme.encode()
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as stream:
        for name, content in sorted(payload.items()):
            info = zipfile.ZipInfo(name, date_time=(2026, 9, 7, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            stream.writestr(info, content)
    checksum = hashlib.sha256(archive.read_bytes()).hexdigest()
    (directory / "SHA256SUMS").write_text(f"{checksum}  {archive.name}\n")
    print(f"Built {archive.name}: {len(files)} payload files; SHA256 {checksum}")


if __name__ == "__main__":
    build()
