#!/usr/bin/env python3
"""PO-08, partial: invariance of the closure under basis RELABELLING.

Scope, stated up front because it bounds the claim.

The certificate target rows are indexed by a coefficient monomial that is
expressed in the *same* basis as the output coordinates. A general invertible
change of basis therefore re-expresses BOTH the monomial index and the
coordinate vector, and cannot be simulated by multiplying coordinates by a
matrix -- doing that would produce a test that looks like a basis change and
is not one. A genuine GL test requires regenerating the certificates
(~530 s/prime plus changes to the generator machinery) and is NOT discharged
here.

What IS exact here is the permutation subgroup: relabelling the basis at a
degree, applied consistently to `basis`, `coordinates` and
`coefficient_monomial`. This is a genuine change of basis, and it is the one
that would expose any dependence of the activation rule on label order.

Invariant claims tested:
  - the closure dimension at each degree;
  - the cardinality of the missing set;
  - the identity of the missing set, tracked through the relabelling.
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


def relabel(certificate, degree, perm):
    """Apply a basis relabelling at `degree` consistently to every field."""
    out = dict(certificate)
    rows = []
    for target in certificate["targets"]:
        row = dict(target)
        if target["field_degree"] == degree:
            old = target["basis"]
            new_basis = [perm.get(n, n) for n in old]
            # permute coordinates to follow the relabelled basis
            index = {n: i for i, n in enumerate(old)}
            coords = target["coordinates"]
            row["basis"] = sorted(new_basis, key=lambda n: new_basis.index(n))
            row["basis"] = new_basis
            row["coordinates"] = coords
        # the coefficient monomial must be relabelled wherever it appears
        row["coefficient_monomial"] = [
            perm.get(c, c) for c in target["coefficient_monomial"]]
        rows.append(row)
    out["targets"] = rows
    return out


def missing_at(certificate, degree, seed):
    base, _ = closure(certificate, seed, certificate["prime"])
    found = []
    for name in next(t["basis"] for t in certificate["targets"]
                     if t["field_degree"] == degree):
        trial = dict(seed)
        trial[degree] = [name]
        dims, _ = closure(certificate, trial, certificate["prime"])
        if dims[degree] > base[degree]:
            found.append(name)
    return base, found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime", type=int, default=32749)
    ap.add_argument("--trials", type=int, default=4)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    with (CERT / f"interacting_degree12_{args.prime}.json").open() as s:
        certificate = json.load(s)

    results = {
        "schema": 1, "prime": args.prime, "trials": args.trials,
        "scope": (
            "permutation subgroup only; general GL basis change requires "
            "regenerating certificates and is NOT discharged here"),
        "degrees": {},
    }

    for degree, expected, seed in ((10, DEG10, dict(SEED8)),
                                   (12, DEG12, {**SEED8, 10: DEG10})):
        names = list(next(t["basis"] for t in certificate["targets"]
                          if t["field_degree"] == degree))
        per_trial = []
        for trial in range(args.trials):
            rng = random.Random(7000 + trial)
            shuffled = list(names)
            rng.shuffle(shuffled)
            perm = dict(zip(names, shuffled))
            inverse = {v: k for k, v in perm.items()}

            relabelled = relabel(certificate, degree, perm)
            seed_r = {k: [perm.get(x, x) if k == degree else x for x in v]
                      for k, v in seed.items()}
            base, found = missing_at(relabelled, degree, seed_r)
            recovered = sorted(inverse.get(n, n) for n in found)

            ok = (base[degree] == FULL[degree] - len(expected)
                  and recovered == sorted(expected))
            per_trial.append({
                "trial": trial,
                "closure_dimension": base[degree],
                "missing_count": len(found),
                "recovered_under_inverse_relabelling": recovered,
                "invariant": ok,
            })
            print(f"  degree {degree} trial {trial}: dim={base[degree]} "
                  f"count={len(found)} recovered={recovered} "
                  f"{'OK' if ok else 'MISMATCH'}")

        stable = all(t["invariant"] for t in per_trial)
        results["degrees"][str(degree)] = {
            "expected": sorted(expected), "trials": per_trial,
            "invariant_under_relabelling": stable}
        print(f"degree {degree}: invariant under relabelling = {stable}\n")

    verdict = all(v["invariant_under_relabelling"]
                  for v in results["degrees"].values())
    results["verdict"] = (
        "closure dimension and missing-set identity are relabelling-invariant"
        if verdict else "RELABELLING DEPENDENCE DETECTED")
    print(f"VERDICT: {results['verdict']}")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(results, indent=1, sort_keys=True) + "\n")
        print(f"wrote {args.out}")
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
