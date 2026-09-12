"""Fault-tolerant modular reconstruction utilities.

The Gaussian lattice candidate follows the standard lattice
< (M,0), (X,1) > construction. Candidates are treated as heuristic until
verified against modular images; verification, not lattice shortness, gates use.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd
from typing import Mapping, Sequence

from .exact import crt
from .modular_reconstruction import (
    asymmetric_profiles,
    profile_consensus_reconstruct,
    verify_fraction,
)


def _norm2(v: tuple[int, int]) -> int:
    return v[0] * v[0] + v[1] * v[1]


def _dot(a: tuple[int, int], b: tuple[int, int]) -> int:
    return a[0] * b[0] + a[1] * b[1]


def _nearest_integer_ratio(num: int, den: int) -> int:
    if den <= 0:
        raise ValueError("positive denominator required")
    # Exact nearest integer, ties away from zero is harmless for 2D Gauss.
    sign = -1 if num < 0 else 1
    x = abs(num)
    q, r = divmod(x, den)
    if 2 * r >= den:
        q += 1
    return sign * q


def gaussian_reduce_2d(u: tuple[int, int], v: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    """Gauss-reduce a rank-two integer lattice basis."""
    if u == (0, 0) or v == (0, 0):
        raise ValueError("nonzero lattice basis vectors required")
    u = (int(u[0]), int(u[1])); v = (int(v[0]), int(v[1]))
    while True:
        if _norm2(v) < _norm2(u):
            u, v = v, u
        q = _nearest_integer_ratio(_dot(u, v), _norm2(u))
        w = (v[0] - q * u[0], v[1] - q * u[1])
        if _norm2(w) >= _norm2(u):
            return u, v
        v = w


def gaussian_rational_candidate(residues_by_prime: Mapping[int, int]) -> Fraction | None:
    """Return a short-lattice rational candidate from all supplied images."""
    primes = tuple(sorted(int(p) for p in residues_by_prime))
    if not primes:
        raise ValueError("at least one modular image required")
    x, modulus = crt([int(residues_by_prime[p]) for p in primes], list(primes))
    shortest, _ = gaussian_reduce_2d((modulus, 0), (x, 1))
    a, b = shortest
    if b == 0:
        return None
    q = Fraction(a, b)
    # Sanity: the candidate vector must lie in the CRT lattice.
    if (q.numerator - x * q.denominator) % modulus:
        return None
    return q


def bad_prime_scores(value: Fraction, residues_by_prime: Mapping[int, int]) -> dict[int, int]:
    """0 for matching images, 1 for rejected/denominator-killing primes."""
    good, bad = verify_fraction(value, residues_by_prime)
    return {int(p): (0 if int(p) in good else 1) for p in residues_by_prime}


def leave_k_out_consensus(
    residues_by_prime: Mapping[int, int],
    *,
    max_bad_primes: int = 1,
    minimum_support: int | None = None,
) -> dict:
    """Search for a low-height rational supported by most modular images.

    Each leave-k-out subset is reconstructed under several asymmetric height
    windows.  All verified profile candidates vote.  The winner must clear the
    modular support threshold and have strictly more reconstruction votes than
    every competitor.  This avoids preferring a huge CRT-fitting fraction merely
    because it matches a corrupted residue as well.
    """
    primes = tuple(sorted(int(p) for p in residues_by_prime))
    if max_bad_primes < 0 or max_bad_primes >= len(primes):
        raise ValueError("invalid bad-prime budget")
    threshold = len(primes) - max_bad_primes if minimum_support is None else int(minimum_support)
    candidates: dict[Fraction, dict] = {}
    for k in range(max_bad_primes + 1):
        for excluded in combinations(primes, k):
            fit = tuple(p for p in primes if p not in excluded)
            result = profile_consensus_reconstruct(residues_by_prime, fit)
            for attempt in result.get("attempts", ()): 
                if attempt.get("status") != "verified" or attempt.get("value") is None:
                    continue
                q = Fraction(attempt["value"])
                good, bad = verify_fraction(q, residues_by_prime)
                if len(good) < threshold:
                    continue
                rec = candidates.setdefault(
                    q,
                    {
                        "votes": 0,
                        "good": good,
                        "bad": bad,
                        "height": max(abs(q.numerator), q.denominator),
                        "subsets": set(),
                    },
                )
                rec["votes"] += 1
                rec["subsets"].add(tuple(excluded))
                if len(good) > len(rec["good"]):
                    rec["good"], rec["bad"] = good, bad
    g = gaussian_rational_candidate(residues_by_prime)
    if g is not None:
        good, bad = verify_fraction(g, residues_by_prime)
        if len(good) >= threshold:
            rec = candidates.setdefault(
                g,
                {
                    "votes": 0,
                    "good": good,
                    "bad": bad,
                    "height": max(abs(g.numerator), g.denominator),
                    "subsets": set(),
                },
            )
            rec["votes"] += 1
    if not candidates:
        return {"status": "unresolved", "candidates": []}
    ordered = sorted(
        candidates.items(),
        key=lambda kv: (kv[1]["votes"], -kv[1]["height"], len(kv[1]["good"])),
        reverse=True,
    )
    best_q, best = ordered[0]
    second_votes = ordered[1][1]["votes"] if len(ordered) > 1 else -1
    if best["votes"] == second_votes:
        return {
            "status": "ambiguous",
            "candidates": [
                {
                    "value": str(q),
                    "votes": rec["votes"],
                    "support": len(rec["good"]),
                    "good": list(rec["good"]),
                    "bad": list(rec["bad"]),
                    "height": rec["height"],
                }
                for q, rec in ordered[:10]
            ],
        }
    return {
        "status": "verified_support_candidate",
        "value": str(best_q),
        "votes": best["votes"],
        "support": len(best["good"]),
        "good_primes": list(best["good"]),
        "bad_primes": list(best["bad"]),
        "scores": bad_prime_scores(best_q, residues_by_prime),
    }
