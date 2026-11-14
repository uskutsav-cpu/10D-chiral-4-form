#!/usr/bin/env python3
"""Emit the machine-readable closure / minimality certificate.

Consolidates the reachability (closure) and minimality results into one
artifact with provenance, so they can be cited without re-running anything.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from stress_flow_closure import closure  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "results" / "stress_flow" / "certificates"
OUT = ROOT / "results" / "stress_flow" / "closure_and_minimality.json"
PRIMES = (32749, 32719, 32717, 32693)
DEGREES = (4, 6, 8, 10, 12)
FULL = {4: 1, 6: 2, 8: 7, 10: 14, 12: 72}
STATIC = {4: 1, 6: 1, 8: 2, 10: 2, 12: 4}
REQUIRED8 = ("I8_3", "I8_4", "I8_5", "I8_6")
INERT8 = ("I8_1", "I8_2", "I4_1^2")


def certificates():
    out = []
    for prime in PRIMES:
        path = CERT / f"interacting_degree12_{prime}.json"
        if path.exists():
            with path.open() as stream:
                out.append((prime, json.load(stream)))
    return out


def agree(seed, certs):
    """Closure dims, asserting every prime agrees."""
    seen = None
    for prime, cert in certs:
        dims, _ = closure(cert, seed, prime)
        if seen is None:
            seen = dims
        elif dims != seen:
            raise AssertionError(f"prime {prime} disagrees: {dims} != {seen}")
    return {str(k): v for k, v in seen.items()}


def main():
    certs = certificates()
    if len(certs) < 3:
        raise SystemExit("need at least three interacting certificates")

    seeds = {
        "free": {4: ["I4_1"]},
        "free_plus_K6": {4: ["I4_1"], 6: ["I6_2"]},
        "free_plus_J6": {4: ["I4_1"], 6: ["I6_1"]},
        "free_plus_K6_and_degree8_complement": {
            4: ["I4_1"], 6: ["I6_2"], 8: list(REQUIRED8)},
        "free_plus_degree8_complement_only": {
            4: ["I4_1"], 8: list(REQUIRED8)},
    }
    closures = {name: agree(seed, certs) for name, seed in seeds.items()}

    # minimality: drop one required degree-8 direction at a time
    removal = {}
    for omitted in REQUIRED8:
        kept = [d for d in REQUIRED8 if d != omitted]
        removal[f"without_{omitted}"] = agree(
            {4: ["I4_1"], 6: ["I6_2"], 8: kept}, certs)
    removal["without_K6"] = closures[
        "free_plus_degree8_complement_only"]

    inert = {}
    base = closures["free_plus_K6"]
    for name in INERT8:
        inert[name] = agree(
            {4: ["I4_1"], 6: ["I6_2"], 8: [name]}, certs) == base

    payload = {
        "schema": 1,
        "claim": (
            "Reachable (closure) dimensions of the pure stress flow by seed, "
            "and minimality of the generalized-flow completion through "
            "degree 8, exact over finite fields and identical on every prime."
        ),
        "caveat": (
            "Closure is reachability under seeding and INCLUDES Tr(tau) "
            "propagation. It is not the same as the new-forcing span in "
            "interacting_flow_equations.json, which excludes Tr(tau). The "
            "two answer different questions and must not be interchanged."
        ),
        "primes": [p for p, _ in certs],
        "samples_per_prime": 4,
        "degrees": list(DEGREES),
        "full_dimension": {str(k): v for k, v in FULL.items()},
        "static_stress_span": {str(k): v for k, v in STATIC.items()},
        "closures": closures,
        "minimality_by_removal": removal,
        "inert_directions_change_nothing": inert,
        "minimal_completion_through_degree8": ["K6", *REQUIRED8],
        "still_open": {"10": 3, "12": 4},
        "reproduction": (
            ".venv/bin/python scripts/emit_closure_certificate.py"
        ),
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as stream:
        json.dump(payload, stream, indent=1, sort_keys=True)
        stream.write("\n")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
