"""Lift stabilized modular row spaces with verification and fail-closed semantics."""
from __future__ import annotations

from fractions import Fraction
from math import gcd
from typing import Mapping, Sequence

from .exact import rref_q
from .finite_field import rref
from .fault_tolerant import leave_k_out_consensus
from .modular_reconstruction import (
    common_denominator,
    profile_consensus_reconstruct,
    reduce_fraction,
)


def canonical_row_spaces(
    rows_by_prime: Mapping[int, Sequence[Sequence[int]]],
    *,
    ncols: int,
    expected_rank: int | None = None,
) -> tuple[dict[int, list[list[int]]], tuple[int, ...]]:
    spaces: dict[int, list[list[int]]] = {}
    pivots_ref: tuple[int, ...] | None = None
    for p, rows in sorted(rows_by_prime.items()):
        rr, pivots = rref(rows, int(p), ncols=ncols)
        basis = rr[:len(pivots)]
        if expected_rank is not None and len(pivots) != expected_rank:
            raise ValueError(f"rank mismatch at prime {p}: {len(pivots)} != {expected_rank}")
        piv = tuple(pivots)
        if pivots_ref is None:
            pivots_ref = piv
        elif piv != pivots_ref:
            raise ValueError(f"pivot pattern differs at prime {p}")
        spaces[int(p)] = basis
    if pivots_ref is None:
        raise ValueError("nonempty modular row spaces required")
    return spaces, pivots_ref


def _residues_for_entry(spaces, row, col) -> dict[int, int]:
    return {int(p): int(matrix[row][col]) for p, matrix in spaces.items()}


def lift_stable_row_space(
    spaces: Mapping[int, Sequence[Sequence[int]]],
    pivots: Sequence[int],
    *,
    ncols: int,
    holdout_primes: Sequence[int] = (),
    max_bad_primes: int = 0,
) -> dict:
    """Lift a canonical modular RREF, verifying every recovered entry.

    Ordinary asymmetric reconstruction is attempted first. If that fails and
    ``max_bad_primes`` is positive, a leave-k-out/Gaussian fault-tolerant
    candidate may be used, but theorem status remains gated by modular support.
    """
    primes = tuple(sorted(int(p) for p in spaces))
    holdout = tuple(sorted(int(p) for p in holdout_primes))
    if any(p not in spaces for p in holdout):
        raise ValueError("holdout prime absent from modular spaces")
    fit = tuple(p for p in primes if p not in holdout)
    if not fit:
        raise ValueError("at least one fit prime required")
    rank = len(pivots)
    nonpivots = tuple(c for c in range(ncols) if c not in pivots)
    rational = [[Fraction(0) for _ in range(ncols)] for _ in range(rank)]
    for r, c in enumerate(pivots):
        rational[r][c] = Fraction(1)
    unresolved = []
    bad_prime_votes: dict[int, int] = {p: 0 for p in primes}
    methods: dict[str, int] = {}
    entry_supports: list[tuple[int, ...]] = []
    for r in range(rank):
        for c in nonpivots:
            residues = _residues_for_entry(spaces, r, c)
            ordinary = profile_consensus_reconstruct(
                residues,
                fit,
                verification_primes=holdout,
            )
            q: Fraction | None = None
            method = None
            support = primes
            if ordinary.get("status") == "verified_unique_consensus":
                q = Fraction(ordinary["value"]); method = "asymmetric_profiles"
            elif max_bad_primes > 0:
                robust = leave_k_out_consensus(residues, max_bad_primes=max_bad_primes)
                if robust.get("status") == "verified_support_candidate":
                    q = Fraction(robust["value"]); method = "fault_tolerant_support"
                    support = tuple(robust["good_primes"])
                    for p in robust["bad_primes"]:
                        bad_prime_votes[int(p)] += 1
            if q is None:
                unresolved.append({"row": r, "column": c, "ordinary_status": ordinary.get("status")})
                continue
            # Exact reduction check on every claimed supporting prime.
            for p in support:
                if reduce_fraction(q, p) != residues[p] % p:
                    raise AssertionError("internal reconstruction verification failure")
            rational[r][c] = q
            entry_supports.append(tuple(sorted(int(p) for p in support)))
            methods[method] = methods.get(method, 0) + 1
    if unresolved:
        return {
            "status": "unresolved",
            "rank": rank,
            "pivots": list(pivots),
            "nonpivot_count": len(nonpivots),
            "unresolved": unresolved,
            "method_counts": methods,
            "bad_prime_votes": bad_prime_votes,
            "qq_candidate": False,
            "qq_verified": False,
            "all_primes_verified": False,
        }
    # Matrix-level normalization: clear denominators per row and record heights.
    row_denominators = [common_denominator(row) for row in rational]
    integer_rows = [[int(q * d) for q in row] for row, d in zip(rational, row_denominators)]
    rr_q, piv_q = rref_q(rational, ncols=ncols)
    if tuple(piv_q) != tuple(pivots):
        raise ValueError("lifted QQ row space has a different pivot pattern")
    all_primes = tuple(primes)
    all_primes_verified = all(support == all_primes for support in entry_supports)
    return {
        "status": "lifted_and_verified" if all_primes_verified else "fault_tolerant_candidate",
        "rank": rank,
        "pivots": list(pivots),
        "row_denominators": row_denominators,
        "integer_rows": integer_rows,
        "rational_rows": [[str(q) for q in row] for row in rational],
        "method_counts": methods,
        "bad_prime_votes": bad_prime_votes,
        "qq_candidate": True,
        "qq_verified": all_primes_verified,
        "all_primes_verified": all_primes_verified,
    }
