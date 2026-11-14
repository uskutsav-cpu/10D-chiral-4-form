"""Closure of the pure stress flow starting from a given seed.

The flow d V/d lambda = f(T[V], lambda) is expanded as a sum over trace
generators with arbitrary lambda-dependent coefficients. Each certificate
target is one (generator, field_degree, coefficient_monomial) row giving the
resulting coordinates in the degree basis.

A target is ACTIVE once every factor of its coefficient monomial lies in the
span already reached. Adding active rows can enlarge the span, which can
activate further targets, so we iterate to a fixed point. The fixed point is
the smallest stress-closed family containing the seed.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sdinv.exactmap import rank_mod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "results" / "stress_flow" / "certificates"
DEGREES = (4, 6, 8, 10, 12)


def monomial_factors(target):
    return [f for f in target["coefficient_monomial"] if f]


def closure(certificate, seed_ids, prime, verbose=False,
            return_span=False):
    """seed_ids: dict degree -> list of basis ids present in the seed."""
    basis = {}
    for target in certificate["targets"]:
        basis[target["field_degree"]] = target["basis"]

    # span[degree] = list of coordinate rows spanning what V can contain
    span = {d: [] for d in DEGREES}
    present = set()
    for degree, ids in seed_ids.items():
        for name in ids:
            row = [0] * len(basis[degree])
            row[basis[degree].index(name)] = 1
            span[degree].append(row)
            present.add(name)

    def in_span(degree, name):
        """Can the coefficient of basis direction `name` be nonzero in V?

        V ranges over the reached span, so its `name` coordinate can be made
        nonzero exactly when some spanning vector has a nonzero entry in that
        column. Requiring the unit vector e_name itself to lie in the span
        would be strictly stronger and would UNDERSTATE the closure.
        """
        if not span[degree]:
            return False
        idx = basis[degree].index(name)
        return any(row[idx] % prime for row in span[degree])

    for iteration in range(24):
        grew = False
        for target in certificate["targets"]:
            degree = target["field_degree"]
            factors = monomial_factors(target)
            # active iff every factor is a direction already reachable
            ok = True
            for f in factors:
                fdeg = next((d for d in DEGREES if f in basis.get(d, [])), None)
                if fdeg is None or not in_span(fdeg, f):
                    ok = False
                    break
            if not ok:
                continue
            row = [c % prime for c in target["coordinates"]]
            if not any(row):
                continue
            before = rank_mod(np.asarray(span[degree] or [[0] * len(basis[degree])],
                                         dtype=np.int64) % prime, prime)
            after = rank_mod(np.asarray(span[degree] + [row],
                                        dtype=np.int64) % prime, prime)
            if after > before:
                span[degree].append(row)
                grew = True
        if not grew:
            if verbose:
                print(f"  fixed point after {iteration} sweeps")
            break

    dims = {d: rank_mod(np.asarray(span[d] or [[0] * len(basis[d])],
                                   dtype=np.int64) % prime, prime)
            for d in DEGREES}
    if return_span:
        return dims, basis, span
    return dims, basis


def closure_span(certificate, seed_ids, prime):
    """Closure together with the spanning rows at each degree.

    Same fixed point as closure(); returns the actual subspace D_d rather than
    only its dimension, which is what the quotient construction Q_d = A_d / D_d
    needs.
    """
    dims, basis, span = closure(
        certificate, seed_ids, prime, return_span=True)
    return dims, basis, span


def main():
    primes = [32749, 32719, 32717, 32693]
    seeds = {
        "free seed  V = c*I4_1": {4: ["I4_1"]},
        "seeded with K6 (I6_2)": {4: ["I4_1"], 6: ["I6_2"]},
        "seeded with J6 (I6_1)": {4: ["I4_1"], 6: ["I6_1"]},
    }
    full = {4: 1, 6: 2, 8: 7, 10: 14, 12: 72}
    static_stress = {4: 1, 6: 1, 8: 2, 10: 2, 12: 4}

    for label, seed in seeds.items():
        print(f"=== {label} ===")
        results = {}
        for prime in primes:
            path = CERT / f"interacting_degree12_{prime}.json"
            if not path.exists():
                continue
            with path.open() as stream:
                cert = json.load(stream)
            dims, _ = closure(cert, seed, prime)
            results[prime] = dims
        ref = results[primes[0]]
        agree = all(r == ref for r in results.values())
        print(f"  {'degree':>7} {'full':>5} {'static':>7} {'closed':>7}")
        for d in DEGREES:
            print(f"  {d:>7} {full[d]:>5} {static_stress[d]:>7} {ref[d]:>7}")
        print(f"  identical across {len(results)} primes: {agree}")
        print()


if __name__ == "__main__":
    main()
