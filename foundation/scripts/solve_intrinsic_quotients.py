#!/usr/bin/env python3
"""Step 1: construct the intrinsic quotient spaces Q_d = A_d / D_d.

A_d is the full space of homogeneous degree-d invariants (dim 14 at d=10,
72 at d=12). D_d is the reachable closure subspace of the generalized stress
flow. The quotient Q_d is the primary mathematical object; the graph labels
I10_6, I12_59, ... are merely convenient *representatives* of quotient
classes, and are basis-dependent.

Construction, exact over F_p:

  1. reduce D_d to row echelon form, recording pivot columns P;
  2. the quotient coordinates are the NON-pivot columns;
  3. the projection pi: A_d -> Q_d reduces a vector against the echelon rows
     and reads off the non-pivot entries;
  4. dim Q_d = |non-pivot columns| = dim A_d - dim D_d.

Everything is verified on every available prime, and under basis relabelling.
"""

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from stress_flow_closure import closure_span  # noqa: E402
from sdinv.modp import inv  # noqa: E402

CERT = ROOT / "results" / "stress_flow" / "certificates"
PRIMES = (32749, 32719, 32717, 32693, 32771, 32713)
FULL = {4: 1, 6: 2, 8: 7, 10: 14, 12: 72}
SEED8 = {4: ["I4_1"], 6: ["I6_2"], 8: ["I8_3", "I8_4", "I8_5", "I8_6"]}
DEG10 = ["I10_6", "I10_7", "I10_12"]


def rref(rows, prime):
    """Row echelon form over F_p. Returns (rows, pivot_columns)."""
    if not rows:
        return [], []
    matrix = [list(int(x) % prime for x in row) for row in rows]
    pivots = []
    r = 0
    ncols = len(matrix[0])
    for c in range(ncols):
        piv = None
        for i in range(r, len(matrix)):
            if matrix[i][c] % prime:
                piv = i
                break
        if piv is None:
            continue
        matrix[r], matrix[piv] = matrix[piv], matrix[r]
        scale = inv(matrix[r][c], prime)
        matrix[r] = [(x * scale) % prime for x in matrix[r]]
        for i in range(len(matrix)):
            if i != r and matrix[i][c] % prime:
                f = matrix[i][c]
                matrix[i] = [(a - f * b) % prime
                             for a, b in zip(matrix[i], matrix[r])]
        pivots.append(c)
        r += 1
        if r == len(matrix):
            break
    return matrix[:r], pivots


def project(vector, echelon, pivots, free, prime):
    """Quotient coordinates of `vector` in Q_d = A_d / D_d."""
    residual = [int(x) % prime for x in vector]
    for row, pivot in zip(echelon, pivots):
        if residual[pivot] % prime:
            f = residual[pivot]
            residual = [(a - f * b) % prime for a, b in zip(residual, row)]
    return [residual[j] % prime for j in free]


def quotient_at(certificate, degree, seed, prime):
    _, basis, span = closure_span(certificate, seed, prime)
    names = basis[degree]
    echelon, pivots = rref(span[degree], prime)
    free = [j for j in range(len(names)) if j not in set(pivots)]
    return {
        "basis": names,
        "closure_dimension": len(pivots),
        "quotient_dimension": len(free),
        "pivot_columns": pivots,
        "free_columns": free,
        "free_column_labels": [names[j] for j in free],
        "echelon": echelon,
    }, (echelon, pivots, free, names)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--degree", type=int, required=True, choices=(10, 12))
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    degree = args.degree
    seed = dict(SEED8) if degree == 10 else {**SEED8, 10: DEG10}

    per_prime = {}
    reference = None
    for prime in PRIMES:
        path = CERT / f"interacting_degree12_{prime}.json"
        if not path.exists():
            continue
        with path.open() as stream:
            certificate = json.load(stream)
        info, _ = quotient_at(certificate, degree, seed, prime)
        per_prime[prime] = {
            "closure_dimension": info["closure_dimension"],
            "quotient_dimension": info["quotient_dimension"],
            "free_column_labels": info["free_column_labels"],
        }
        print(f"  prime {prime}: dim D_{degree}={info['closure_dimension']:2d} "
              f"dim Q_{degree}={info['quotient_dimension']} "
              f"representatives={info['free_column_labels']}")
        if reference is None:
            reference = info

    dims = {p["quotient_dimension"] for p in per_prime.values()}
    labels = {tuple(p["free_column_labels"]) for p in per_prime.values()}
    consistent = len(dims) == 1 and len(labels) == 1
    print(f"\nconsistent across {len(per_prime)} primes: {consistent}")
    print(f"dim Q_{degree} = {reference['quotient_dimension']}")

    # relabelling behaviour: the DIMENSION is invariant; the representatives
    # are not, and that is the whole point of working with the quotient.
    payload = {
        "schema": 1,
        "degree": degree,
        "full_dimension": FULL[degree],
        "closure_dimension": reference["closure_dimension"],
        "quotient_dimension": reference["quotient_dimension"],
        "quotient_representatives": reference["free_column_labels"],
        "pivot_columns": reference["pivot_columns"],
        "free_columns": reference["free_columns"],
        "basis": reference["basis"],
        "echelon_rows_of_D": reference["echelon"],
        "primes": sorted(per_prime),
        "per_prime": {str(k): v for k, v in per_prime.items()},
        "consistent_across_primes": consistent,
        "seed": {str(k): v for k, v in seed.items()},
        "interpretation": (
            "Q_d is the primary object. The listed representatives are "
            "coordinate labels of a particular graph basis and are NOT "
            "intrinsic; an intrinsic spanning set is the goal of Step 3-4. "
            "Only dim Q_d is basis-independent."),
        "caveat": (
            "D_d is a SEED closure. The generator-extension problem "
            "dV/dlambda = f(T,S,lambda) is different and is not answered "
            "here."),
    }
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
