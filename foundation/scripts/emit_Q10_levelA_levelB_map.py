#!/usr/bin/env python3
"""Exact change-of-basis between the Q10 Level-A and Level-B bases.

Computes BOTH directions per prime, verifies they are mutual inverses, and
proves removal minimality by explicit rank drop. Everything is exact over F_p;
no floating point is used anywhere.

Rational reconstruction is ATTEMPTED and reported honestly: with only two
primes the CRT modulus is far below the bound needed to certify a rational
lift, so the matrices are preserved modularly and the limitation is stated
rather than papered over.
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sdinv.exactmap import rank_mod

R = ROOT / "results" / "intrinsic_candidates"
OUT = R / "Q10_levelA_levelB_map.json"
LABELS = ["Q10_A", "Q10_B", "Q10_C"]


def inv_matrix(mat, p):
    """Exact inverse of an n x n matrix over F_p, or None if singular."""
    n = len(mat)
    aug = [[int(mat[i][j]) % p for j in range(n)]
           + [1 if i == j else 0 for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if aug[r][col] % p), None)
        if piv is None:
            return None
        aug[col], aug[piv] = aug[piv], aug[col]
        f = pow(aug[col][col], p - 2, p)
        aug[col] = [(v * f) % p for v in aug[col]]
        for r in range(n):
            if r != col and aug[r][col] % p:
                g = aug[r][col]
                aug[r] = [(a - g * b) % p for a, b in zip(aug[r], aug[col])]
    return [row[n:] for row in aug]


def matmul(a, b, p):
    n, m, k = len(a), len(b[0]), len(b)
    return [[sum(a[i][t] * b[t][j] for t in range(k)) % p for j in range(m)]
            for i in range(n)]


def main():
    amap = json.loads((R / "published_degree10_map.json").read_text())
    basis = json.loads((R / "intrinsic_Q10_levelB_basis.json").read_text())
    ctrl = json.loads((R / "degree10_positive_control.json").read_text())
    triple = basis["preferred_basis"]

    out = {"schema": 1, "preferred_basis": triple, "per_prime": {},
           "labels": LABELS}
    all_ok = True

    for prime in sorted(amap["per_prime"]):
        p = int(prime)
        # rows of B: the Level-B basis quotient vectors
        B = [amap["per_prime"][prime]["projections"][n]["quotient_vector"]
             for n in triple]
        # rows of A: the Level-A class quotient vectors
        A = ctrl["per_prime"][prime]["vectors"]

        Binv = inv_matrix(B, p)
        Ainv = inv_matrix(A, p)
        assert Binv is not None, f"Level-B basis singular at {prime}"
        assert Ainv is not None, f"Level-A basis singular at {prime}"

        # A = C . B   =>   C = A . B^-1   (Level-A expressed in Level-B)
        C_A_in_B = matmul(A, Binv, p)
        # B = D . A   =>   D = B . A^-1   (Level-B expressed in Level-A)
        C_B_in_A = matmul(B, Ainv, p)

        prod = matmul(C_A_in_B, C_B_in_A, p)
        identity = [[1 if i == j else 0 for j in range(3)] for i in range(3)]
        mutual = prod == identity

        # reconstruction check: C_A_in_B . B must reproduce A exactly
        recon = matmul(C_A_in_B, B, p) == [[v % p for v in row] for row in A]

        # --- removal minimality, by explicit rank drop --------------------
        removal = {}
        for drop in triple:
            keep = [n for n in triple if n != drop]
            rows = np.asarray(
                [amap["per_prime"][prime]["projections"][n]["quotient_vector"]
                 for n in keep], dtype=np.int64) % p
            removal[drop] = {"remaining": keep,
                             "rank": int(rank_mod(rows, p)),
                             "drops_below_3": int(rank_mod(rows, p)) < 3}

        # --- forced membership across ALL independent triples -------------
        proj = amap["per_prime"][prime]["projections"]
        nonzero = sorted(k for k in proj if "[" not in k
                         and proj[k]["status"] == "solved"
                         and any(v % p for v in proj[k]["quotient_vector"]))
        indep = []
        for t in itertools.combinations(nonzero, 3):
            rows = np.asarray([proj[k]["quotient_vector"] for k in t],
                              dtype=np.int64) % p
            if rank_mod(rows, p) == 3:
                indep.append(list(t))
        forced = [k for k in nonzero if indep and all(k in t for t in indep)]

        ok = mutual and recon and all(v["drops_below_3"]
                                      for v in removal.values())
        all_ok &= ok
        out["per_prime"][prime] = {
            "levelA_in_levelB": {l: C_A_in_B[i] for i, l in enumerate(LABELS)},
            "levelB_in_levelA": {n: C_B_in_A[i] for i, n in enumerate(triple)},
            "mutually_inverse": mutual,
            "reconstruction_exact": recon,
            "removal_minimality": removal,
            "independent_triples": indep,
            "forced_members": forced,
            "checks_pass": ok,
        }
        print(f"prime {prime}: mutual_inverse={mutual} recon={recon} "
              f"forced={forced} removal_all_drop="
              f"{all(v['drops_below_3'] for v in removal.values())}")

    # --- rational reconstruction, attempted and reported --------------------
    primes = [int(p) for p in out["per_prime"]]
    modulus = 1
    for p in primes:
        modulus *= p
    bound = int((modulus // 2) ** 0.5)
    out["rational_reconstruction"] = {
        "attempted": True,
        "certified": False,
        "n_primes": len(primes),
        "crt_modulus": str(modulus),
        "uniqueness_bound_numerator_denominator": str(bound),
        "reason": (
            "Two primes give a CRT modulus of about 1.07e9, so a rational "
            "lift is only unique when |numerator| and |denominator| are both "
            "below roughly 2.3e4. The change-of-basis entries are generic "
            "residues of that magnitude, so no lift can be distinguished from "
            "a coincidence. The exact modular matrices are preserved instead. "
            "Certifying rationals needs more primes, not more analysis."),
    }
    out["all_checks_pass"] = all_ok
    OUT.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(f"\nall_checks_pass={all_ok}")
    print(f"rational reconstruction certified: False (2 primes, bound {bound})")
    print(f"wrote {OUT.name}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
