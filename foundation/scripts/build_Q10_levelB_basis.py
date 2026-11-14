#!/usr/bin/env python3
"""Select a preferred compact Level-B basis for Q10 from the twelve candidates.

Reads the committed projection artifact, enumerates every independent
three-element subset, scores them under a documented simplicity criterion, and
proves removal minimality.

The claim this supports is deliberately narrow: the selected triple is
*preferred under the criterion recorded here* and *minimal among the twelve
published candidates*. It is NOT claimed to be canonical in any absolute sense
-- that would need a proof over a search space that has not been enumerated.
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sdinv.exactmap import rank_mod
from sdinv.published_degree10_invariants import PUBLISHED_DEGREE10

MAP = ROOT / "results" / "intrinsic_candidates" / "published_degree10_map.json"
OUT_BASIS = ROOT / "results" / "intrinsic_candidates" / "intrinsic_Q10_levelB_basis.json"
OUT_SEARCH = ROOT / "results" / "intrinsic_candidates" / "Q10_basis_search.json"

# Measured single-evaluation wall time, seconds, at prime 32749 (seed 5).
RUNTIME = {"P10_01": 0.0, "P10_02": 7.3, "P10_03": 5.8, "P10_04": 5.3,
           "P10_05": 4.9, "P10_06": 10.5, "P10_07": 5.0, "P10_08": 5.3,
           "P10_09": 5.1, "P10_10": 5.4, "P10_11": 7.2, "P10_12": 6.6}


def _solve3(rows, target, p):
    """Solve  x . rows = target  over F_p by Gauss-Jordan. Exact, no floats.

    `rows` are the basis quotient vectors as the ROWS of a 3x3 matrix, so the
    system is rows^T x = target.
    """
    n = len(rows)
    aug = [[int(rows[j][i]) % p for j in range(n)] + [int(target[i]) % p]
           for i in range(len(target))]
    piv_row = 0
    where = []
    for col in range(n):
        sel = next((r for r in range(piv_row, len(aug)) if aug[r][col] % p), None)
        if sel is None:
            where.append(None)
            continue
        aug[piv_row], aug[sel] = aug[sel], aug[piv_row]
        inv_p = pow(aug[piv_row][col], p - 2, p)
        aug[piv_row] = [(v * inv_p) % p for v in aug[piv_row]]
        for r in range(len(aug)):
            if r != piv_row and aug[r][col] % p:
                f = aug[r][col]
                aug[r] = [(a - f * b) % p for a, b in zip(aug[r], aug[piv_row])]
        where.append(piv_row)
        piv_row += 1
    for r in range(piv_row, len(aug)):
        if aug[r][n] % p:
            return None                      # inconsistent
    out = [0] * n
    for col in range(n):
        if where[col] is not None:
            out[col] = aug[where[col]][n] % p
    return out


def score(name, ambiguity_robust):
    """Deterministic simplicity score. LOWER is simpler.

    Weights are stated here rather than tuned: block count dominates, explicit
    bracket work is next, and runtime breaks ties. The one judgement encoded is
    that a candidate whose quotient image is INVARIANT under an unresolved
    source ambiguity is strictly preferable to one whose image moves -- a basis
    element that changes when the ambiguity is resolved is not a stable answer.
    """
    spec = PUBLISHED_DEGREE10[name]
    blocks = spec["blocks"]
    n_m = blocks.get("M", 0)
    n_1050 = blocks.get("N1050", 0)
    n_4125 = blocks.get("N4125", 0)
    brackets = spec.get("brackets", "")
    n_red = 1 if "RED" in brackets or "red" in brackets else 0
    n_black = brackets.count("BracketOp") + (
        1 if "explicit" in brackets else 0)
    # a RED symmetrisation over k slots expands to k! terms
    perm_terms = 6 if n_red else 1
    penalty = 0 if ambiguity_robust else 50
    return {
        "name": name,
        "m_factors": n_m, "n1050_factors": n_1050, "n4125_factors": n_4125,
        "black_ops": n_black, "red_ops": n_red,
        "expanded_permutation_terms": perm_terms,
        "runtime_seconds": RUNTIME.get(name),
        "ambiguity": spec.get("ambiguity"),
        "ambiguity_robust": ambiguity_robust,
        "score": (10 * (n_m + n_1050 + n_4125) + 5 * n_black + 8 * n_red
                  + perm_terms + (RUNTIME.get(name) or 0) + penalty),
    }


def main():
    payload = json.loads(MAP.read_text())
    per_prime = payload["per_prime"]

    # robustness: does the ambiguity variant land on the SAME quotient vector?
    robust = {}
    for prime, rec in per_prime.items():
        proj = rec["projections"]
        for key in list(proj):
            if "[" not in key:
                continue
            base = key.split("[")[0]
            if base not in proj or proj[base]["status"] != "solved":
                continue
            same = proj[base]["quotient_vector"] == proj[key]["quotient_vector"]
            robust[base] = robust.get(base, True) and same

    results = {}
    for prime, rec in sorted(per_prime.items()):
        p = int(prime)
        proj = rec["projections"]
        names = sorted(k for k in proj if "[" not in k
                       and proj[k]["status"] == "solved")
        vecs = {k: proj[k]["quotient_vector"] for k in names}
        nonzero = [k for k in names if any(v % p for v in vecs[k])]
        dim = rec["dim_Q10"]

        independent = []
        for triple in itertools.combinations(nonzero, 3):
            rows = np.asarray([vecs[k] for k in triple], dtype=np.int64) % p
            if rank_mod(rows, p) == 3:
                independent.append(list(triple))

        # removal minimality: which candidates appear in EVERY independent
        # triple? Those cannot be dropped from any basis at all.
        forced = [k for k in nonzero
                  if independent and all(k in t for t in independent)]

        results[prime] = {
            "dim_Q10": dim,
            "rank": rec["Q10_rank_from_published"],
            "nonzero_candidates": nonzero,
            "independent_triples": independent,
            "forced_members": forced,
        }
        print(f"prime {prime}: rank {rec['Q10_rank_from_published']}/{dim}, "
              f"nonzero {nonzero}, {len(independent)} independent triples, "
              f"forced {forced}")

    # consistency across primes before selecting anything
    triple_sets = [{tuple(t) for t in r["independent_triples"]}
                   for r in results.values()]
    consistent = all(s == triple_sets[0] for s in triple_sets)
    print(f"\nindependent-triple set agrees across primes: {consistent}")

    candidates = sorted(triple_sets[0]) if consistent and triple_sets else []
    scored = []
    for triple in candidates:
        parts = [score(n, robust.get(n, True)) for n in triple]
        scored.append({"triple": list(triple),
                       "total": sum(x["score"] for x in parts),
                       "members": parts})
    scored.sort(key=lambda x: (x["total"], x["triple"]))

    for entry in scored:
        flags = ",".join(
            f"{m['name']}{'' if m['ambiguity_robust'] else '*'}"
            for m in entry["members"])
        print(f"  {entry['total']:8.1f}  {flags}")

    preferred = scored[0] if scored else None

    # --- Level-A <-> Level-B change of basis --------------------------------
    # Express each Level-A quotient class Q10_A/B/C in the preferred Level-B
    # basis, by exact modular solve. Both sides are vectors in the SAME
    # 3-dimensional quotient coordinates, so this is a 3x3 solve per prime.
    control = json.loads(
        (ROOT / "results" / "intrinsic_candidates"
         / "degree10_positive_control.json").read_text())
    level_map = {}
    for prime, rec in sorted(per_prime.items()):
        p = int(prime)
        if preferred is None or prime not in control["per_prime"]:
            continue
        basis_rows = np.asarray(
            [per_prime[prime]["projections"][n]["quotient_vector"]
             for n in preferred["triple"]], dtype=object) % p
        entry = {}
        for label, target in zip(("Q10_A", "Q10_B", "Q10_C"),
                                 control["per_prime"][prime]["vectors"]):
            coeffs = _solve3(basis_rows, np.asarray(target, dtype=object) % p, p)
            entry[label] = None if coeffs is None else [int(c) for c in coeffs]
        level_map[prime] = entry
        print(f"  Level-A in Level-B at {prime}: "
              + ", ".join(f"{k}={v}" for k, v in sorted(entry.items())))
    basis = {
        "schema": 1,
        "claim": ("preferred under the simplicity criterion recorded in "
                  "scripts/build_Q10_levelB_basis.py, and minimal among the "
                  "twelve published equation-(4.24) candidates"),
        "not_claimed": ("canonical in any absolute sense; no larger class was "
                        "enumerated"),
        "preferred_basis": preferred["triple"] if preferred else None,
        "preferred_score": preferred["total"] if preferred else None,
        "ambiguity_robust": {k: bool(v) for k, v in sorted(robust.items())},
        "per_prime": results,
        "triples_agree_across_primes": consistent,
        "levelA_in_levelB": level_map,
    }
    OUT_BASIS.write_text(json.dumps(basis, indent=1, sort_keys=True) + "\n")
    OUT_SEARCH.write_text(json.dumps(
        {"schema": 1, "scored_triples": scored,
         "runtime_source": "measured at prime 32749, seed 5"},
        indent=1, sort_keys=True) + "\n")
    print(f"\npreferred basis: {basis['preferred_basis']}")
    print(f"wrote {OUT_BASIS.name}, {OUT_SEARCH.name}")


if __name__ == "__main__":
    main()
