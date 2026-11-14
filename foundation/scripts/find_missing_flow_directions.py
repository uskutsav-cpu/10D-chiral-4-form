#!/usr/bin/env python3
"""Locate the directions a generalized stress flow cannot reach.

For a given seed, a basis direction is *missing* if adjoining it to the seed
strictly raises the iterative closure at its own degree. The missing set at a
degree has cardinality equal to the closure deficit there, and each element is
checked for non-redundancy by removal.

Verified on every available prime; a disagreement aborts rather than being
averaged away.
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from stress_flow_closure import closure  # noqa: E402

CERT = ROOT / "results" / "stress_flow" / "certificates"
PRIMES = (32749, 32719, 32717, 32693)
FULL = {4: 1, 6: 2, 8: 7, 10: 14, 12: 72}


def certificates():
    out = {}
    for prime in PRIMES:
        path = CERT / f"interacting_degree12_{prime}.json"
        if path.exists():
            with path.open() as stream:
                out[prime] = json.load(stream)
    if len(out) < 3:
        raise SystemExit("need at least three interacting certificates")
    return out


def basis_at(certificate, degree):
    return next(t["basis"] for t in certificate["targets"]
                if t["field_degree"] == degree)


def agree(certs, seed):
    """Closure dims; abort if primes disagree."""
    seen = None
    for prime, certificate in certs.items():
        dims, _ = closure(certificate, seed, prime)
        if seen is None:
            seen = dims
        elif dims != seen:
            raise SystemExit(
                f"PRIME DISAGREEMENT at {prime}: {dims} != {seen}")
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--degree", type=int, required=True, choices=(10, 12))
    ap.add_argument("--single-prime-scan", action="store_true",
                    help="scan candidates on one prime, then confirm the "
                         "selected set on all primes")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    certs = certificates()
    first = certs[PRIMES[0]]
    degree = args.degree

    seed = {4: ["I4_1"], 6: ["I6_2"],
            8: ["I8_3", "I8_4", "I8_5", "I8_6"]}
    if degree == 12:
        seed[10] = ["I10_6", "I10_7", "I10_12"]

    base = agree(certs, seed)
    deficit = FULL[degree] - base[degree]
    print(f"base closure: {base}")
    print(f"degree {degree} deficit: {deficit}")

    scan_certs = {PRIMES[0]: first} if args.single_prime_scan else certs
    raisers = []
    for name in basis_at(first, degree):
        trial = dict(seed)
        trial[degree] = [name]
        dims = agree(scan_certs, trial)
        if dims[degree] > base[degree]:
            raisers.append(name)
            print(f"  RAISES {name}")

    print(f"raisers: {raisers} (count {len(raisers)})")

    joint = dict(seed)
    joint[degree] = raisers
    closed = agree(certs, joint)
    print(f"with all raisers: {closed}")
    fully_closed = closed[degree] == FULL[degree]
    print(f"degree {degree} closed: {fully_closed}")

    removal = {}
    if fully_closed:
        for omitted in raisers:
            trial = dict(seed)
            trial[degree] = [x for x in raisers if x != omitted]
            dims = agree(certs, trial)
            removal[omitted] = dims[degree]
            print(f"  without {omitted}: {dims[degree]} "
                  f"({'REOPENS' if dims[degree] < FULL[degree] else 'NOT NEEDED'})")

    payload = {
        "schema": 1,
        "degree": degree,
        "primes_confirmed": list(certs),
        "seed": {str(k): v for k, v in seed.items()},
        "base_closure": {str(k): v for k, v in base.items()},
        "deficit": deficit,
        "missing_directions": raisers,
        "closure_with_all": {str(k): v for k, v in closed.items()},
        "degree_closed": fully_closed,
        "removal_minimality": removal,
        "caveat": (
            "These are SEED directions. Adjoining a flow GENERATOR is a "
            "different operation and is not established by this computation."),
    }
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
