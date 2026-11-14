#!/usr/bin/env python3
"""Check saved exact witnesses; --recompute independently reruns every graph."""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sdinv.certificate import digest, integral_directions, verify_witness


def verify(manifest, certificate, recompute=False):
    directions, indices = integral_directions()
    if certificate.get("schema") != 1 or certificate.get("complete") is not True:
        raise ValueError("certificate is incomplete or unsupported")
    if certificate["basis_sha256"] != digest(manifest):
        raise ValueError("certificate basis hash mismatch")
    if certificate["integral_directions_sha256"] != digest(directions.tolist()) or certificate["coordinate_component_indices"] != indices:
        raise ValueError("coordinate convention mismatch")
    cells = [[w["prime"], w["seed"]] for w in certificate["witnesses"]]
    if cells != certificate["planned_cells"] or len({tuple(c) for c in cells}) != len(cells):
        raise ValueError("missing, reordered or duplicate witness cells")
    if len({p for p, _ in cells}) < 2 or len({s for _, s in cells}) < 2:
        raise ValueError("release gate requires at least two primes and two seeds")
    return [verify_witness(manifest, witness, recompute) for witness in certificate["witnesses"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis", default=str(ROOT / "results/rank81_basis.json"))
    parser.add_argument("--certificate", default=str(ROOT / "results/rank81_certificate.json"))
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(Path(args.basis).read_text())
    certificate = json.loads(Path(args.certificate).read_text())
    results = verify(manifest, certificate, args.recompute)
    print(json.dumps({"passed": True, "mode": "graph_recomputation" if args.recompute else "saved_matrix_verification", "witnesses": results}, indent=2))


if __name__ == "__main__":
    main()
