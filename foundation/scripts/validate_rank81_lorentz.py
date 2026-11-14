#!/usr/bin/env python3
"""Check every fixed graph under an exact spatial rotation and a genuine boost."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
from sdinv.catalog import atomic_write_json
from sdinv.certificate import context, digest, validate_basis
from sdinv.contract import planned_value
from sdinv.serialize import from_record
from sdinv.published_degree8_invariants import evaluate_selected_basis


def transform(form, matrix, prime):
    for axis in range(5):
        # A sum of ten reduced products fits exactly in int64 for our primes.
        form = np.moveaxis(np.tensordot(matrix, form, axes=(1, axis)), 0, axis) % prime
    return form


def main():
    manifest = json.loads((ROOT / "results/rank81_basis.json").read_text())
    certificate = json.loads((ROOT / "results/rank81_certificate.json").read_text())
    if not certificate["complete"]:
        raise ValueError("finish the rank certificate first")
    items = validate_basis(manifest)
    witness = certificate["witnesses"][0]
    prime = witness["prime"]
    form, _, _ = context(prime, witness["coordinates"])
    eta = np.diag([-1] + [1] * 9)
    rotations = {}
    rotation = np.eye(10, dtype=np.int64)
    t = 5
    c, s = (1 - t*t) * pow(1 + t*t, -1, prime) % prime, 2*t * pow(1 + t*t, -1, prime) % prime
    rotation[2, 2], rotation[3, 3], rotation[2, 3], rotation[3, 2] = c, c, s, -s % prime
    rotations["spatial_rotation"] = rotation
    boost = np.eye(10, dtype=np.int64)
    t = 7
    c, s = (1 + t*t) * pow(1 - t*t, -1, prime) % prime, 2*t * pow(1 - t*t, -1, prime) % prime
    boost[0, 0], boost[1, 1], boost[0, 1], boost[1, 0] = c, c, s, s
    rotations["lorentz_boost"] = boost
    tensor_values = evaluate_selected_basis(form, prime)
    result = {"schema": 1, "prime": prime, "seed": witness["seed"],
              "basis_sha256": digest(manifest), "checks": {}, "passed": False}
    for name, matrix in rotations.items():
        if not np.array_equal(matrix.T @ eta @ matrix % prime, eta % prime):
            raise ValueError("test matrix is not Lorentz")
        transformed = transform(form, matrix, prime)
        matched = []
        for item, expected in zip(items, witness["values"]):
            scalar = planned_value(from_record(item["graph"]), transformed, 10, 5, True, prime,
                                   max_memory_bytes=2 * 1024**3)
            if scalar != expected:
                raise ValueError(f"{name} fails for {item['id']}")
            matched.append(item["id"])
        if evaluate_selected_basis(transformed, prime) != tensor_values:
            raise ValueError(f"{name} fails for selected literature basis")
        result["checks"][name] = {"matrix": matrix.tolist(), "matched_invariants": matched,
                                   "literature_order8_basis_matched": True}
        print(f"PASS: {name}, all {len(matched)} graphs and six literature tensors", flush=True)
    result["passed"] = True
    atomic_write_json(ROOT / "results/rank81_lorentz.json", result)


if __name__ == "__main__":
    main()
