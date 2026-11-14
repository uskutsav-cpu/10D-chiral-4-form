#!/usr/bin/env python3
"""Rebuild a graph-only rank-81 witness from the existing discovered basis.

Discovery is already complete: degree10_pipeline.py and degree12_pipeline.py
retain their resumable catalog searches. This final stage evaluates exactly
the same 81 graphs in each field. Checkpoints retain every accepted raw row.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path
import sys

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
from sdinv.catalog import atomic_write_json
from sdinv.certificate import (digest, evaluate_item, integral_directions,
                               make_witness, require_prime, validate_basis, verify_witness)
from sdinv.checkpoint import load_checkpoint, write_checkpoint
from sdinv.modp import RankSieve


def _job(args):
    return evaluate_item(*args)


def evaluate_cell(manifest, prime, seed, checkpoint_dir, workers):
    require_prime(prime)
    coordinates = np.random.default_rng(seed).integers(0, prime, 126).tolist()
    engine = {str(path.relative_to(ROOT)): path.read_bytes().hex() for path in
              [ROOT / "src/sdinv" / f"{name}.py" for name in
               ("certificate", "contract", "forms", "modp", "serialize", "graphs")]}
    identity = {"schema": 1, "basis_sha256": digest(manifest), "prime": prime,
                "seed": seed, "coordinates": coordinates, "engine_sha256": digest(engine)}
    checkpoint = checkpoint_dir / f"p{prime}_s{seed}.json"
    state = load_checkpoint(checkpoint, identity) if checkpoint.exists() else {"evaluated": [], "pivots": [], "rank": 0}
    items = validate_basis(manifest)
    evaluated = state["evaluated"]
    if [r["id"] for r in evaluated] != [x["id"] for x in items[:len(evaluated)]]:
        raise ValueError("checkpoint is not a prefix of the fixed basis")
    sieve = RankSieve(126, prime)
    for row in evaluated:
        if not sieve.add(row["jacobian_row"]):
            raise ValueError("checkpoint contains dependent rows")
    if state["rank"] != sieve.rank or state["pivots"] != sieve.pivots:
        raise ValueError("checkpoint pivot metadata mismatch")
    # Batches bound both the submission queue and completed rows awaiting the
    # coordinator. Worker count is explicit because contractions can use GiB.
    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    try:
        for start in range(len(evaluated), len(items), workers):
            jobs = [(item, prime, coordinates) for item in items[start:start + workers]]
            output = pool.map(_job, jobs) if pool else map(_job, jobs)
            for result in output:
                if not sieve.add(result["jacobian_row"]):
                    raise ValueError(f"fixed basis lost rank at {result['id']}")
                evaluated.append(result)
                state = {"evaluated": evaluated, "pivots": sieve.pivots, "rank": sieve.rank}
                write_checkpoint(checkpoint, identity, state)
                if sieve.rank % 10 == 0 or sieve.rank == 81:
                    print(f"p={prime} seed={seed}: {sieve.rank}/81", flush=True)
    finally:
        if pool:
            pool.shutdown(wait=True, cancel_futures=True)
    witness = make_witness(manifest, prime, seed, coordinates, evaluated)
    verify_witness(manifest, witness)
    return witness


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis", default=str(ROOT / "results/rank81_basis.json"))
    parser.add_argument("--out", default=str(ROOT / "results/rank81_certificate.json"))
    parser.add_argument("--checkpoint-dir", type=Path, default=ROOT / "work/rank81")
    parser.add_argument("--primes", type=int, nargs="+", default=[32749, 32719])
    parser.add_argument("--seeds", type=int, nargs="+", default=[20260907, 20260908])
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1 or len(set(args.primes)) != len(args.primes) or len(set(args.seeds)) != len(args.seeds):
        parser.error("workers must be positive and primes/seeds must be distinct")
    manifest = json.loads(Path(args.basis).read_text())
    validate_basis(manifest)
    directions, indices = integral_directions()
    certificate = {"schema": 1, "basis_sha256": digest(manifest),
                   "coordinate_convention": "A_I = F_I for sorted 5-tuples I containing 0; F = sum_I A_I (e_I + *e_I).",
                   "coordinate_component_indices": indices,
                   "integral_directions_sha256": digest(directions.tolist()),
                   "proof": "Each graph is an integer polynomial in A. A nonzero 81x81 Jacobian minor modulo a prime proves generic characteristic-zero rank at least 81. The known quotient dimension supplies the upper bound; this is not a full invariant-ring presentation or a global coordinate chart.",
                   "upper_bound_source": "https://arxiv.org/abs/2509.14350v2",
                   "planned_cells": [[p, s] for p in args.primes for s in args.seeds],
                   "complete": False, "witnesses": []}
    for prime, seed in certificate["planned_cells"]:
        certificate["witnesses"].append(evaluate_cell(manifest, prime, seed, args.checkpoint_dir, args.workers))
        atomic_write_json(args.out, certificate)
    certificate["complete"] = True
    atomic_write_json(args.out, certificate)
    print(f"PASS: fixed 81 graphs, {len(certificate['witnesses'])} exact nonzero minors")


if __name__ == "__main__":
    main()
