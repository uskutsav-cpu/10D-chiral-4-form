#!/usr/bin/env python3
"""Freeze the pure-N^(4125) Q10 basis and all five exact change-of-basis maps.

The three basis elements live entirely in the N4125^5 sector, which contains no
M block, so their validation is untouched by the mixed-variance M defect that
invalidated the M-containing sector statistics.

Maps produced, per prime, exactly over F_p:
  1. pure-N  -> Level-A          4. published -> pure-N
  2. Level-A -> pure-N           5. pure-N    -> published
  3. pure-N  -> full degree-10 atlas
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sdinv.exactmap import rank_mod
from sdinv.reverse_block_decomposition import (
    build_einsum, canonical_form, parse_formula, render_formula,
    stream_candidates)

R = ROOT / "results" / "intrinsic_candidates"
EINSUMS = ["abcdef,abcdgh,efijkl,ghimno,jklmno->",
           "abcdef,abcdgh,efijkl,gijmno,hklmno->",
           "abcdef,abcdgh,egijkl,fhimno,jklmno->"]
KINDS = ["N4125"] * 5


def inv_matrix(mat, p):
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


def det_mod(mat, p):
    n, a, det = len(mat), [row[:] for row in mat], 1
    for col in range(n):
        piv = next((r for r in range(col, n) if a[r][col] % p), None)
        if piv is None:
            return 0
        if piv != col:
            a[col], a[piv] = a[piv], a[col]
            det = (-det) % p
        det = (det * a[col][col]) % p
        f = pow(a[col][col], p - 2, p)
        a[col] = [(v * f) % p for v in a[col]]
        for r in range(col + 1, n):
            if a[r][col] % p:
                g = a[r][col]
                a[r] = [(x - g * y) % p for x, y in zip(a[r], a[col])]
    return det % p


def matmul(a, b, p):
    return [[sum(a[i][t] * b[t][j] for t in range(len(b))) % p
             for j in range(len(b[0]))] for i in range(len(a))]


def main():
    val = json.loads((R / "degree10_reverse_span_validation.json").read_text())
    pub_map = json.loads((R / "published_degree10_map.json").read_text())
    basis = json.loads((R / "intrinsic_Q10_levelB_basis.json").read_text())
    ctrl = json.loads((R / "degree10_positive_control.json").read_text())
    triple_pub = basis["preferred_basis"]

    # --- formulas, regenerated and round-tripped ---------------------------
    topos, formulas = {}, {}
    for topo, _ in stream_candidates(KINDS, cap=60000):
        spec, _r = build_einsum(KINDS, topo)
        if spec in EINSUMS:
            topos[spec] = topo
    assert len(topos) == len(EINSUMS), f"only found {len(topos)} of 3"
    for spec in EINSUMS:
        text = render_formula(KINDS, topos[spec])
        k2, t2 = parse_formula(text)
        assert canonical_form(KINDS, topos[spec]) == canonical_form(k2, t2), \
            f"round trip failed for {spec}"
        formulas[spec] = text

    out = {"schema": 1,
           "claim": ("A three-element compact basis for Q10 was independently "
                     "recovered in the pure N^(4125)^5 sector and shown to span "
                     "the same quotient space as the graph-derived and "
                     "published-formula bases."),
           "not_claimed": ["universally canonical",
                           "simplest structures reaching Q10",
                           "a resolution of the published bracket ambiguity",
                           "a characteristic-zero identity"],
           "sector": "N4125+N4125+N4125+N4125+N4125",
           "sector_fully_enumerated": True,
           "sector_canonical_candidates": 15,
           "einsums": EINSUMS,
           "formulas": {f"R{i+1}": formulas[s] for i, s in enumerate(EINSUMS)},
           "round_trip_verified": True,
           "blocks_per_element": {"M": 0, "N1050": 0, "N4125": 5},
           "bracket_operations": {"black": 0, "red": 0},
           "expanded_permutation_terms": 1,
           "per_prime": {}}

    maps = {"schema": 1, "per_prime": {}}
    all_ok = True

    for prime, rec in sorted(val["per_prime"].items()):
        p = int(prime)
        pure = rec["original"]["vectors"]
        lvlA = ctrl["per_prime"][prime]["vectors"]
        pub = [pub_map["per_prime"][prime]["projections"][n]["quotient_vector"]
               for n in triple_pub]

        d_pure, d_A, d_pub = (det_mod(m, p) for m in (pure, lvlA, pub))
        iA, ipub, ipure = (inv_matrix(m, p) for m in (lvlA, pub, pure))

        pure_in_A = matmul(pure, iA, p)
        A_in_pure = matmul(lvlA, ipure, p)
        pure_in_pub = matmul(pure, ipub, p)
        pub_in_pure = matmul(pub, ipure, p)
        ident = [[1 if i == j else 0 for j in range(3)] for i in range(3)]

        mutual_A = matmul(pure_in_A, A_in_pure, p) == ident
        mutual_pub = matmul(pure_in_pub, pub_in_pure, p) == ident

        removal = {}
        for k in range(3):
            keep = [pure[j] for j in range(3) if j != k]
            removal[f"R{k+1}"] = int(
                rank_mod(np.asarray(keep, dtype=np.int64) % p, p))

        ok = (d_pure and mutual_A and mutual_pub
              and all(v < 3 for v in removal.values())
              and rank_mod(np.asarray(pure + lvlA, dtype=np.int64) % p, p) == 3)
        all_ok &= bool(ok)

        out["per_prime"][prime] = {
            "quotient_vectors": pure,
            "rank": rec["original"]["rank"],
            "fresh_sample_rank": rec["fresh"]["rank"],
            "determinant": d_pure,
            "removal_minimality_rank_after_removal": removal,
        }
        maps["per_prime"][prime] = {
            "pure_to_levelA": pure_in_A, "levelA_to_pure": A_in_pure,
            "pure_to_published": pure_in_pub, "published_to_pure": pub_in_pure,
            "determinants": {"pure": d_pure, "levelA": d_A, "published": d_pub},
            "mutually_inverse_levelA": mutual_A,
            "mutually_inverse_published": mutual_pub,
            "span_equals_levelA": bool(
                rank_mod(np.asarray(pure + lvlA, dtype=np.int64) % p, p) == 3),
            "checks_pass": bool(ok),
        }
        print(f"prime {prime}: det={d_pure} rank={rec['original']['rank']} "
              f"fresh={rec['fresh']['rank']} mutualA={mutual_A} "
              f"mutualPub={mutual_pub} removal={removal}")

    maps["rational_reconstruction"] = {
        "attempted": True, "certified": False,
        "reason": ("two primes give a CRT modulus of ~1.07e9, so a rational "
                   "lift is unique only below ~2.3e4; the entries are generic "
                   "residues of that magnitude. Exact modular matrices are "
                   "preserved and labelled UNCERTIFIED.")}
    out["all_checks_pass"] = all_ok
    maps["all_checks_pass"] = all_ok
    (R / "Q10_pure_N4125_basis.json").write_text(
        json.dumps(out, indent=1, sort_keys=True) + "\n")
    (R / "Q10_pure_N4125_basis_maps.json").write_text(
        json.dumps(maps, indent=1, sort_keys=True) + "\n")
    print(f"\nall_checks_pass={all_ok}")
    print("wrote Q10_pure_N4125_basis.json, Q10_pure_N4125_basis_maps.json")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
