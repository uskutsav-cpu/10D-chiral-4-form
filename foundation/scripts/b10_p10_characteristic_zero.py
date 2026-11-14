#!/usr/bin/env python3
"""Exact characteristic-zero dimension of `B10 + P10`, hence of `B10 cap P10`.

`dim(B10 cap P10) = dim B10 + dim P10 - dim(B10 + P10)`. The first two terms are
already pinned over `Q`: each space is spanned by exactly as many explicit
invariants as its modular rank, which forces equality. Only `dim(B10 + P10)` was
modular, so only

    dim_Q(B10 cap P10) <= 1

followed. This lifts the remaining term.

`P10` is spanned by two coordinate directions and is exact by construction. What
has to be lifted is `B10`: the atlas coordinates of the twelve published
candidates, which come from solving a linear system modulo `p`. The same
technique as `d10_characteristic_zero.py` applies -- CRT across primes, rational
reconstruction, validation at a held-out prime -- and it either produces a
certified rational coordinate matrix or it does not, in which case the bound
stands and says so.

    python scripts/b10_p10_characteristic_zero.py --solve   # recompute per prime
    python scripts/b10_p10_characteristic_zero.py           # lift and conclude
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from math import gcd
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

COORDS = ROOT / "results" / "degree10" / "B10_coordinates_per_prime.json"
OUT = ROOT / "results" / "degree10" / "B10_P10_intersection_exact.json"

from d10_characteristic_zero import crt, rational_reconstruct, rref  # noqa: E402


def solve_at_prime(prime: int) -> dict:
    """Atlas coordinates of every published candidate, modulo one prime."""
    from sdinv.forms import selfdual_projector, to_dense, random_form
    from sdinv.published_degree10_invariants import PUBLISHED_DEGREE10
    from stress_flow_closure import closure_span
    from test_M_only_quotients import (
        registry_items, evaluate_atlas_element, solve_exact,
    )

    # The degree-10 atlas ordering is prime-independent, so it is taken from
    # whatever is already recorded rather than from a per-prime closure
    # certificate. That decoupling is what allows extra primes to be added:
    # certificates exist for six primes, and six were not enough to lift the
    # three highest-height candidates.
    if COORDS.exists():
        store = json.loads(COORDS.read_text())
        if store.get("per_prime"):
            names = store["per_prime"][sorted(store["per_prime"])[0]]["basis"]
        else:
            names = None
    else:
        names = None
    if names is None:
        cert_path = (ROOT / "results" / "stress_flow" / "certificates" /
                     f"interacting_degree12_{prime}.json")
        if not cert_path.exists():
            return {}
        cert = json.loads(cert_path.read_text())
        seed = {4: ["I4_1"], 6: ["I6_2"], 8: ["I8_3", "I8_4", "I8_5", "I8_6"]}
        _, bmap, _ = closure_span(cert, seed, prime)
        names = bmap[10]

    items = registry_items()
    proj = selfdual_projector(10, 5, True, prime)
    forms = [to_dense((proj @ random_form(10, 5, np.random.default_rng(41000 + 11 * i),
                                          prime)) % prime, 10, 5, prime)
             for i in range(len(names) + 8)]
    cols = [[evaluate_atlas_element(items[nm], items, f, prime, {}) for f in forms]
            for nm in names]
    design = [[cols[j][i] for j in range(len(names))] for i in range(len(forms))]

    out = {"basis": names, "coordinates": {}, "unsolved": []}
    for cid, spec in PUBLISHED_DEGREE10.items():
        target = [spec["evaluator"](f, prime) for f in forms]
        vec, ok = solve_exact(design, target, prime)
        if ok:
            out["coordinates"][cid] = [int(v) % prime for v in vec]
        else:
            out["unsolved"].append(cid)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--solve", action="store_true",
                    help="recompute the per-prime coordinates (slow)")
    ap.add_argument("--primes", default="32693,32713,32717,32719,32749,32771")
    args = ap.parse_args()
    primes = [int(x) for x in args.primes.split(",") if x.strip()]

    if args.solve:
        store = json.loads(COORDS.read_text()) if COORDS.exists() else {"per_prime": {}}
        for p in primes:
            if str(p) in store["per_prime"]:
                print(f"  prime {p}: already recorded, skipped")
                continue
            rec = solve_at_prime(p)
            if not rec:
                print(f"  prime {p}: no closure certificate, skipped")
                continue
            store["per_prime"][str(p)] = rec
            COORDS.parent.mkdir(parents=True, exist_ok=True)
            COORDS.write_text(json.dumps(store, indent=1) + "\n")
            print(f"  prime {p}: solved {len(rec['coordinates'])}, "
                  f"unsolved {len(rec['unsolved'])}")
        return 0

    if not COORDS.exists():
        print(f"missing {COORDS.relative_to(ROOT)}; run with --solve first")
        return 1
    store = json.loads(COORDS.read_text())
    have = sorted(int(k) for k in store["per_prime"])
    if len(have) < 3:
        print(f"only {len(have)} prime(s) recorded; need at least 3 to lift")
        return 1

    fit, holdout = have[:-1], have[-1]
    basis = store["per_prime"][str(have[0])]["basis"]
    ids = sorted(store["per_prime"][str(have[0])]["coordinates"])

    lifted, failed, mismatches = {}, [], []
    for cid in ids:
        if not all(cid in store["per_prime"][str(p)]["coordinates"] for p in have):
            failed.append(f"{cid}: missing at some prime")
            continue
        cols = list(zip(*[store["per_prime"][str(p)]["coordinates"][cid] for p in fit]))
        vec, ok = [], True
        for c in cols:
            x, m = crt(list(c), fit)
            f = rational_reconstruct(x, m)
            if f is None:
                ok = False
                break
            vec.append(f)
        if not ok:
            failed.append(f"{cid}: no small rational lift")
            continue
        ref = store["per_prime"][str(holdout)]["coordinates"][cid]
        bad = False
        for j, f in enumerate(vec):
            den = f.denominator % holdout
            if den == 0:
                continue
            if (f.numerator * pow(den, -1, holdout)) % holdout != ref[j] % holdout:
                bad = True
                break
        if bad:
            mismatches.append(cid)
            continue
        lifted[cid] = vec

    n = len(basis)
    prod_idx = [i for i, nm in enumerate(basis) if nm.startswith("I4_1*")]
    p_rows = [[Fraction(1 if i == j else 0) for j in range(n)] for i in prod_idx]
    b_rows = [lifted[c] for c in sorted(lifted)]

    settled = not failed and not mismatches and len(lifted) == len(ids)
    rank_lifted = len(rref(b_rows)[1]) if b_rows else 0
    dim_p = len(rref(p_rows)[1])
    sum_lifted = len(rref(b_rows + p_rows)[1]) if b_rows else dim_p
    cap = rank_lifted + dim_p - sum_lifted

    print(f"primes: fitting {fit}, holdout {holdout}")
    print(f"lifted {len(lifted)}/{len(ids)} published candidates"
          + (f", failed {failed}" if failed else "")
          + (f", holdout mismatches {mismatches}" if mismatches else ""))
    if settled:
        print(f"over Q: dim B10 = {rank_lifted}, dim P10 = {dim_p}, "
              f"dim(B10+P10) = {sum_lifted}, dim(B10 cap P10) = {cap}")
        print("STATUS: exact")
    else:
        # These are ranks of the SUBSET that lifted. Calling them dim B10 would
        # be false: three unlifted rows can only raise the span.
        print(f"over Q, LIFTED SUBSET ONLY: rank {rank_lifted}, "
              f"rank with P10 {sum_lifted}, intersection of the subset {cap}")
        print("STATUS: NOT SETTLED — these are not dim B10 or dim(B10 cap P10); "
              "the modular bound stands")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "schema": 1,
        "generated_by": "scripts/b10_p10_characteristic_zero.py",
        "exact_domain": "Q, via CRT + rational reconstruction then Fraction arithmetic",
        "fitting_primes": fit,
        "holdout_prime": holdout,
        "n_published": len(ids),
        "n_lifted": len(lifted),
        "failed": failed,
        "holdout_mismatches": mismatches,
        "dim_B10_over_Q": rank_lifted if settled else None,
        "dim_P10_over_Q": dim_p,
        "dim_B10_plus_P10_over_Q": sum_lifted if settled else None,
        "dim_B10_cap_P10_over_Q": cap if settled else None,
        "lifted_subset_rank": rank_lifted,
        "lifted_subset_plus_P10_rank": sum_lifted,
        "lifted_subset_intersection": cap,
        "why_subset_values_are_not_dimensions": (
            "the unlifted rows can only enlarge the span, so the subset's rank "
            "is a lower bound on dim B10 and its intersection with P10 is not "
            "comparable to dim(B10 cap P10) in either direction"),
        "settled": settled,
        "permitted_wording": (
            f"dim_Q(B10 cap P10) = {cap}, established by lifting the published "
            f"coordinates to Q and validating at a held-out prime"
            if settled else
            "dim_Q(B10 cap P10) <= 1; the lift did not certify, so the modular "
            "value is not promoted to a characteristic-zero equality"),
    }, indent=1) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
