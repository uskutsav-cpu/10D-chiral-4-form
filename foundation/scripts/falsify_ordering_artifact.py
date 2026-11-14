#!/usr/bin/env python3
"""Falsification test 5: is the missing-direction set an ordering artifact?

The four degree-12 missing directions are the LAST FOUR labels in the atlas
ordering (I12_59..I12_62). That is exactly what an artifact of candidate
ordering would look like, so it must be excluded before any intrinsic search
targets those directions.

Two independent perturbations:

  1. shuffle the order in which certificate target rows are consumed by the
     closure fixed point;
  2. shuffle the order in which candidate directions are scanned.

The closure is a fixed point of "adjoin every activatable row", so its span
should not depend on either ordering. If it does, the result is an artifact
and the Phase 1 targets are wrong.
"""

import argparse
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from stress_flow_closure import closure  # noqa: E402

CERT = ROOT / "results" / "stress_flow" / "certificates"
FULL = {4: 1, 6: 2, 8: 7, 10: 14, 12: 72}
SEED8 = {4: ["I4_1"], 6: ["I6_2"], 8: ["I8_3", "I8_4", "I8_5", "I8_6"]}
DEG10 = ["I10_6", "I10_7", "I10_12"]
DEG12 = ["I12_59", "I12_60", "I12_61", "I12_62"]


def load(prime):
    with (CERT / f"interacting_degree12_{prime}.json").open() as stream:
        return json.load(stream)


def basis_at(certificate, degree):
    return next(t["basis"] for t in certificate["targets"]
                if t["field_degree"] == degree)


def scan(certificate, degree, seed, order):
    """Directions whose adjunction raises the closure at `degree`."""
    base, _ = closure(certificate, seed, order["prime"])
    found = []
    for name in order["candidates"]:
        trial = dict(seed)
        trial[degree] = [name]
        dims, _ = closure(certificate, trial, order["prime"])
        if dims[degree] > base[degree]:
            found.append(name)
    return base, found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime", type=int, default=32749)
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    certificate = load(args.prime)
    results = {"schema": 1, "prime": args.prime, "trials": args.trials,
               "degrees": {}}

    for degree, expected, seed in (
        (10, DEG10, dict(SEED8)),
        (12, DEG12, {**SEED8, 10: DEG10}),
    ):
        candidates = list(basis_at(certificate, degree))
        per_trial = []
        for trial in range(args.trials):
            rng = random.Random(1000 + trial)

            # perturbation 1: consume certificate rows in a shuffled order
            shuffled = dict(certificate)
            rows = list(certificate["targets"])
            rng.shuffle(rows)
            shuffled["targets"] = rows

            # perturbation 2: scan candidates in a shuffled order
            order_candidates = list(candidates)
            rng.shuffle(order_candidates)

            base, found = scan(
                shuffled, degree, seed,
                {"prime": args.prime, "candidates": order_candidates})
            per_trial.append({
                "trial": trial,
                "base_closure_at_degree": base[degree],
                "found_sorted": sorted(found),
                "matches_expected": sorted(found) == sorted(expected),
            })
            print(f"  degree {degree} trial {trial}: base={base[degree]} "
                  f"found={sorted(found)} "
                  f"{'OK' if sorted(found) == sorted(expected) else 'MISMATCH'}")

        stable = all(t["matches_expected"] for t in per_trial)
        results["degrees"][str(degree)] = {
            "expected": sorted(expected),
            "trials": per_trial,
            "stable_under_reordering": stable,
        }
        print(f"degree {degree}: stable under reordering = {stable}\n")

    verdict = all(v["stable_under_reordering"]
                  for v in results["degrees"].values())
    results["verdict"] = (
        "NOT an ordering artifact" if verdict else
        "ORDERING ARTIFACT DETECTED -- Phase 1 targets are invalid")
    print(f"VERDICT: {results['verdict']}")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(results, indent=1, sort_keys=True) + "\n")
        print(f"wrote {args.out}")
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
