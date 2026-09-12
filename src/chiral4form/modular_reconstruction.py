"""Robust modular-to-rational reconstruction helpers.

This module deliberately separates *candidate recovery* from *certification*.
A rational candidate is never accepted unless its reduction is checked against
all requested verification primes.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from fractions import Fraction
from itertools import product
from math import gcd, isqrt
from typing import Iterable, Mapping, Sequence

from .exact import crt, rational_reconstruct


@dataclass(frozen=True)
class ReconstructionAttempt:
    status: str
    value: str | None
    numerator_bound: int
    denominator_bound: int
    fit_primes: tuple[int, ...]
    verified_primes: tuple[int, ...]
    rejected_primes: tuple[int, ...]
    modulus: int
    method: str = "bounded_asymmetric"

    def to_json(self) -> dict:
        return asdict(self)


def reduce_fraction(value: Fraction, prime: int) -> int:
    """Reduce a rational number modulo ``prime`` with denominator guard."""
    q = Fraction(value)
    if q.denominator % prime == 0:
        raise ValueError(f"prime {prime} divides reconstructed denominator")
    return q.numerator * pow(q.denominator, -1, prime) % prime


def verify_fraction(
    value: Fraction,
    residues_by_prime: Mapping[int, int],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return (verified, rejected) prime sets for a rational candidate."""
    good: list[int] = []
    bad: list[int] = []
    q = Fraction(value)
    for p, residue in sorted(residues_by_prime.items()):
        try:
            actual = reduce_fraction(q, p)
        except ValueError:
            bad.append(int(p))
            continue
        (good if actual == int(residue) % p else bad).append(int(p))
    return tuple(good), tuple(bad)


def symmetric_bound(modulus: int) -> int:
    if modulus <= 2:
        raise ValueError("modulus must exceed two")
    return isqrt((int(modulus) - 1) // 2)


def asymmetric_profiles(
    modulus: int,
    *,
    anchors: Sequence[int] = (1, 10**3, 10**6, 10**9, 10**12, 10**15, 10**18, 10**21, 10**24),
    include_symmetric: bool = True,
    safety_divisor: int = 4,
) -> list[tuple[int, int]]:
    """Generate deterministic (N,D) windows satisfying ``2*N*D < modulus``.

    The profiles explicitly test denominator-heavy and numerator-heavy
    possibilities instead of assuming equal heights.
    """
    m = int(modulus)
    if m <= 2 or safety_divisor < 3:
        raise ValueError("invalid reconstruction profile parameters")
    profiles: list[tuple[int, int]] = []
    if include_symmetric:
        b = symmetric_bound(m)
        # Stay below the uniqueness boundary.
        b = max(1, (b * 9) // 10)
        profiles.append((b, b))
    for n in anchors:
        n = int(n)
        if n < 1:
            continue
        d = (m - 1) // (safety_divisor * n)
        if d >= 1 and 2 * n * d < m:
            profiles.append((n, d))
    for d in anchors:
        d = int(d)
        if d < 1:
            continue
        n = (m - 1) // (safety_divisor * d)
        if n >= 1 and 2 * n * d < m:
            profiles.append((n, d))
    # Preserve order, remove duplicate windows.
    return list(dict.fromkeys(profiles))


def reconstruct_asymmetric(
    residues_by_prime: Mapping[int, int],
    fit_primes: Sequence[int],
    numerator_bound: int,
    denominator_bound: int,
    *,
    verification_primes: Sequence[int] = (),
) -> ReconstructionAttempt:
    """Reconstruct one fraction and verify it on disjoint primes."""
    fit = tuple(int(p) for p in fit_primes)
    if not fit or len(set(fit)) != len(fit):
        raise ValueError("fit primes must be nonempty and distinct")
    verify = tuple(int(p) for p in verification_primes)
    if set(fit) & set(verify):
        raise ValueError("verification primes must be disjoint from fit primes")
    if any(p not in residues_by_prime for p in fit + verify):
        raise ValueError("missing modular residue")
    a, modulus = crt([int(residues_by_prime[p]) for p in fit], list(fit))
    try:
        q = rational_reconstruct(
            a,
            modulus,
            int(numerator_bound),
            int(denominator_bound),
        )
    except ValueError:
        return ReconstructionAttempt(
            status="no_candidate",
            value=None,
            numerator_bound=int(numerator_bound),
            denominator_bound=int(denominator_bound),
            fit_primes=fit,
            verified_primes=(),
            rejected_primes=verify,
            modulus=modulus,
        )
    good, bad = verify_fraction(q, {p: residues_by_prime[p] for p in verify})
    status = "verified" if not bad else "holdout_rejected"
    return ReconstructionAttempt(
        status=status,
        value=str(q),
        numerator_bound=int(numerator_bound),
        denominator_bound=int(denominator_bound),
        fit_primes=fit,
        verified_primes=good,
        rejected_primes=bad,
        modulus=modulus,
    )


def profile_consensus_reconstruct(
    residues_by_prime: Mapping[int, int],
    fit_primes: Sequence[int],
    *,
    verification_primes: Sequence[int] = (),
    profiles: Sequence[tuple[int, int]] | None = None,
) -> dict:
    """Try many asymmetric height windows and require candidate consensus."""
    fit = tuple(int(p) for p in fit_primes)
    _, modulus = crt([int(residues_by_prime[p]) for p in fit], list(fit))
    windows = list(profiles) if profiles is not None else asymmetric_profiles(modulus)
    attempts: list[ReconstructionAttempt] = []
    candidates: dict[Fraction, list[ReconstructionAttempt]] = {}
    for n, d in windows:
        attempt = reconstruct_asymmetric(
            residues_by_prime,
            fit,
            n,
            d,
            verification_primes=verification_primes,
        )
        attempts.append(attempt)
        if attempt.status == "verified" and attempt.value is not None:
            q = Fraction(attempt.value)
            candidates.setdefault(q, []).append(attempt)
    if len(candidates) == 1:
        value = next(iter(candidates))
        return {
            "status": "verified_unique_consensus",
            "value": str(value),
            "supporting_profiles": len(candidates[value]),
            "attempts": [a.to_json() for a in attempts],
        }
    return {
        "status": "ambiguous" if candidates else "unresolved",
        "values": sorted(str(q) for q in candidates),
        "attempts": [a.to_json() for a in attempts],
    }


def lcm(a: int, b: int) -> int:
    return abs(a // gcd(a, b) * b) if a and b else 0


def common_denominator(values: Iterable[Fraction]) -> int:
    d = 1
    for q in values:
        d = lcm(d, Fraction(q).denominator)
    return d


def sage_style_echelon_height_bound(
    candidate_rows: Sequence[Sequence[Fraction]],
    *,
    source_height: int,
    modulus_product: int,
    ncols: int,
) -> dict:
    """Evaluate Sage's sufficient multimodular echelon proof inequality.

    After clearing a common denominator d from the candidate echelon form E,
    check ``H(dE) * ncols(A) * H(A) < product(primes)``.  The caller must
    provide a valid characteristic-zero source-matrix height bound H(A).
    """
    if source_height < 1 or modulus_product < 2 or ncols < 1:
        raise ValueError("positive source height, modulus, and column count required")
    flat = [Fraction(x) for row in candidate_rows for x in row]
    d = common_denominator(flat)
    height = max((abs(int(q * d)) for q in flat), default=0)
    lhs = height * int(ncols) * int(source_height)
    return {
        "common_denominator": d,
        "cleared_candidate_height": height,
        "source_height": int(source_height),
        "ncols": int(ncols),
        "modulus_product": int(modulus_product),
        "lhs": lhs,
        "proved": lhs < int(modulus_product),
    }
