"""Exact finite-field change-of-basis and obstruction certificates."""

from fractions import Fraction
from math import gcd, isqrt

import numpy as np

from .forms import random_form, selfdual_projector, to_dense
from .modp import inv
from .stress import (
    five_form_moment,
    matrix_trace_power,
    symmetric_inner,
    stress_correction_r,
    stress_mixed,
)


DEFAULT_PRIMES = (32749, 32719, 32693, 32717, 32771)
DEFAULT_SAMPLE_SEEDS = tuple(range(20260901, 20260919))


def rank_mod(matrix, mod):
    """Exact rank of a small dense matrix over ``F_mod``."""
    work = np.asarray(matrix, dtype=np.int64).copy() % mod
    if work.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    row = 0
    for column in range(work.shape[1]):
        pivots = np.nonzero(work[row:, column])[0]
        if not len(pivots):
            continue
        pivot = row + int(pivots[0])
        if pivot != row:
            work[[row, pivot]] = work[[pivot, row]]
        work[row] = work[row] * inv(work[row, column], mod) % mod
        for other in range(work.shape[0]):
            if other != row and work[other, column]:
                work[other] = (
                    work[other] - int(work[other, column]) * work[row]
                ) % mod
        row += 1
        if row == work.shape[0]:
            break
    return int(row)


def solve_full_column_rank(matrix, target, mod):
    """Solve ``matrix @ x = target`` over ``F_mod``.

    The coefficient matrix may be tall but must have full column rank, which
    makes the solution unique.  Inconsistency is reported explicitly.
    """
    matrix = np.asarray(matrix, dtype=np.int64) % mod
    target = np.asarray(target, dtype=np.int64).reshape(-1) % mod
    if matrix.ndim != 2 or matrix.shape[0] != target.shape[0]:
        raise ValueError("incompatible matrix and target dimensions")
    rows, columns = matrix.shape
    augmented = np.column_stack((matrix, target)) % mod
    pivot_rows = {}
    row = 0
    for column in range(columns):
        candidates = np.nonzero(augmented[row:, column])[0]
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        if pivot != row:
            augmented[[row, pivot]] = augmented[[pivot, row]]
        augmented[row] = (
            augmented[row] * inv(augmented[row, column], mod)
        ) % mod
        for other in range(rows):
            if other != row and augmented[other, column]:
                augmented[other] = (
                    augmented[other]
                    - int(augmented[other, column]) * augmented[row]
                ) % mod
        pivot_rows[column] = row
        row += 1
        if row == rows:
            break
    if len(pivot_rows) != columns:
        raise ValueError(
            f"basis matrix has rank {len(pivot_rows)}, expected {columns}")
    for residual in augmented:
        if not np.any(residual[:columns]) and residual[columns]:
            raise ValueError("target is outside the supplied basis span")
    solution = np.zeros(columns, dtype=np.int64)
    for column, pivot_row in pivot_rows.items():
        solution[column] = augmented[pivot_row, columns]
    if not np.array_equal(matrix @ solution % mod, target):
        raise AssertionError("internal modular solve verification failed")
    return solution


def independent_row_indices(matrix, mod):
    """Deterministically select a maximal independent set of matrix rows."""
    matrix = np.asarray(matrix, dtype=np.int64) % mod
    selected = []
    current_rank = 0
    for index in range(matrix.shape[0]):
        candidate = selected + [index]
        next_rank = rank_mod(matrix[candidate], mod)
        if next_rank > current_rank:
            selected.append(index)
            current_rank = next_rank
        if current_rank == matrix.shape[1]:
            break
    return selected


def _crt_pair(a, modulus_a, b, modulus_b):
    """Chinese-remainder combination for coprime moduli."""
    if gcd(modulus_a, modulus_b) != 1:
        raise ValueError("CRT moduli must be coprime")
    step = ((b - a) * pow(modulus_a, -1, modulus_b)) % modulus_b
    modulus = modulus_a * modulus_b
    return int((a + modulus_a * step) % modulus), int(modulus)


def crt(residues, moduli):
    residues = [int(x) for x in residues]
    moduli = [int(x) for x in moduli]
    if not residues or len(residues) != len(moduli):
        raise ValueError("residue/modulus count mismatch")
    value, modulus = residues[0] % moduli[0], moduli[0]
    for residue, next_modulus in zip(residues[1:], moduli[1:]):
        value, modulus = _crt_pair(
            value, modulus, residue % next_modulus, next_modulus)
    return value, modulus


def rational_reconstruct(residue, modulus):
    """Recover the unique small rational represented modulo ``modulus``.

    Numerator and denominator are bounded by ``floor(sqrt(modulus/2))``.
    This is sufficient for the low-degree integer-normalized contractions;
    every reconstructed result is rechecked under each prime.
    """
    residue = int(residue) % int(modulus)
    modulus = int(modulus)
    if residue == 0:
        return Fraction(0, 1)
    bound = isqrt(modulus // 2)
    old_r, r = modulus, residue
    old_t, t = 0, 1
    while r and r > bound:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_t, t = t, old_t - quotient * t
    if r == 0 or t == 0:
        raise ValueError("rational reconstruction failed")
    numerator, denominator = int(r), int(t)
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    common = gcd(abs(numerator), denominator)
    numerator //= common
    denominator //= common
    if abs(numerator) > bound or denominator > bound:
        raise ValueError("reconstructed fraction exceeds uniqueness bound")
    if gcd(denominator, modulus) != 1:
        raise ValueError("reconstructed denominator is not invertible")
    if (numerator * pow(denominator, -1, modulus)) % modulus != residue:
        raise ValueError("reconstructed fraction has the wrong residue")
    return Fraction(numerator, denominator)


def reconstruct_vector(vectors, primes):
    """CRT and rationally reconstruct matching modular vectors."""
    if len(vectors) != len(primes):
        raise ValueError("one coefficient vector is required per prime")
    lengths = {len(vector) for vector in vectors}
    if len(lengths) != 1:
        raise ValueError("coefficient vectors have different lengths")
    result = []
    for column in range(len(vectors[0])):
        residue, modulus = crt(
            [int(vector[column]) for vector in vectors], primes)
        try:
            fraction = rational_reconstruct(residue, modulus)
        except ValueError as exc:
            raise ValueError(
                f"rational reconstruction failed at column {column}: {exc}"
            ) from exc
        for vector, prime in zip(vectors, primes):
            expected = (
                fraction.numerator
                * inv(fraction.denominator, prime)
            ) % prime
            if expected != int(vector[column]) % prime:
                raise AssertionError(
                    "rational reconstruction failed prime verification")
        result.append(fraction)
    return result


def fraction_record(value):
    value = Fraction(value)
    return {
        "numerator": int(value.numerator),
        "denominator": int(value.denominator),
        "text": str(value),
    }


def sample_selfdual_five_form(seed, prime):
    projector = selfdual_projector(10, 5, True, prime)
    rng = np.random.default_rng(int(seed))
    compact = (
        projector @ random_form(10, 5, rng, prime)
    ) % prime
    return to_dense(compact, 10, 5, prime)


def trace_targets(five_form, mod, extra_targets=()):
    """Evaluate the unnormalized stress-building contractions.

    The free stress tensor is ``M/48``.  Keeping the M-normalized targets
    gives small rational change-of-basis coefficients; the physical powers
    of 48 are attached separately to the artifact.
    """
    _, mixed = five_form_moment(five_form, mod)
    traces = {
        power: matrix_trace_power(mixed, power, mod)
        for power in range(2, 11)
    }
    result = {
        "tr_M2": traces[2],
        "tr_M3": traces[3],
        "tr_M4": traces[4],
        "tr_M5": traces[5],
        "tr_M2^2": traces[2] * traces[2] % mod,
        "tr_M2*tr_M3": traces[2] * traces[3] % mod,
    }
    extra_targets = set(extra_targets)
    if {"tr_M2R", "paper_I8"} & extra_targets:
        r_tensor = stress_correction_r(five_form, mod)
        r_mixed = stress_mixed(r_tensor, mod)
        if "paper_I8" in extra_targets:
            m_lower, _ = five_form_moment(five_form, mod)
            result["paper_I8"] = symmetric_inner(
                m_lower, r_tensor, mod)
        if "tr_M2R" in extra_targets:
            m2 = mixed @ mixed % mod
            result["tr_M2R"] = int(np.trace(m2 @ r_mixed % mod) % mod)
    return result


STRESS_TARGETS = {
    4: ("tr_M2",),
    6: ("tr_M3",),
    8: ("tr_M4", "tr_M2^2"),
    10: ("tr_M5", "tr_M2*tr_M3"),
}


def compute_degree_map(registry, degree, primes=DEFAULT_PRIMES,
                       sample_seeds=DEFAULT_SAMPLE_SEEDS,
                       extra_targets=()):
    """Compute and reconstruct a degree-specific exact basis map."""
    degree = int(degree)
    targets = tuple(STRESS_TARGETS[degree]) + tuple(extra_targets)
    items = registry.basis(degree)
    if len(sample_seeds) <= len(items):
        raise ValueError(
            f"degree {degree} requires more than {len(items)} samples "
            "to reserve at least one held-out validation row"
        )
    per_prime = {}
    modular_solutions = {target: [] for target in targets}

    for prime in primes:
        basis_rows = []
        target_columns = {target: [] for target in targets}
        for seed in sample_seeds:
            five_form = sample_selfdual_five_form(seed, prime)
            basis_rows.append(
                registry.evaluate_degree(degree, five_form, prime))
            values = trace_targets(
                five_form,
                prime,
                extra_targets=targets,
            )
            for target in targets:
                target_columns[target].append(values[target])
        basis_matrix = np.asarray(basis_rows, dtype=np.int64) % prime
        basis_rank = rank_mod(basis_matrix, prime)
        if basis_rank != len(items):
            raise RuntimeError(
                f"degree {degree}, prime {prime}: sampled basis rank "
                f"{basis_rank}, expected {len(items)}")
        fit_indices = independent_row_indices(basis_matrix, prime)
        if len(fit_indices) != len(items):
            raise RuntimeError("failed to select a full-rank fit sample set")
        validation_indices = [
            index for index in range(len(sample_seeds))
            if index not in set(fit_indices)
        ]
        solutions = {}
        for target in targets:
            solution = solve_full_column_rank(
                basis_matrix[fit_indices],
                np.asarray(target_columns[target], dtype=np.int64)[fit_indices],
                prime,
            )
            all_values = (
                basis_matrix @ solution % prime
            )
            if not np.array_equal(
                    all_values,
                    np.asarray(target_columns[target], dtype=np.int64) % prime):
                raise AssertionError(
                    f"{target} failed held-out sample validation")
            solutions[target] = [int(x) for x in solution]
            modular_solutions[target].append(solution)
        per_prime[str(prime)] = {
            "sample_count": len(sample_seeds),
            "sample_seeds": [int(seed) for seed in sample_seeds],
            "fit_sample_seeds": [
                int(sample_seeds[index]) for index in fit_indices],
            "held_out_sample_seeds": [
                int(sample_seeds[index]) for index in validation_indices],
            "held_out_validation_passed": True,
            "basis_rank": basis_rank,
            "sample_values": {
                "basis": [
                    [int(value) for value in row]
                    for row in basis_matrix
                ],
                "targets": {
                    target: [int(value) for value in target_columns[target]]
                    for target in targets
                },
            },
            "solutions": solutions,
        }

    reconstructed = {
        target: reconstruct_vector(modular_solutions[target], primes)
        for target in targets
    }
    return {
        "degree": degree,
        "basis": [item.id for item in items],
        "basis_dimension": len(items),
        "targets": {
            target: [fraction_record(value) for value in coefficients]
            for target, coefficients in reconstructed.items()
        },
        "primes": [int(prime) for prime in primes],
        "per_prime": per_prime,
    }


def select_standard_complement(stress_rows):
    """Select standard-coordinate rows which complement ``stress_rows``."""
    rational_rows = [
        [Fraction(value) for value in row] for row in stress_rows]
    if not rational_rows:
        raise ValueError("at least one stress row is required")
    width = len(rational_rows[0])
    if any(len(row) != width for row in rational_rows):
        raise ValueError("stress rows have inconsistent widths")

    # Rational Gaussian rank is tiny here; Fraction keeps this proof exact.
    def rational_rank(rows):
        work = [list(row) for row in rows]
        pivot_row = 0
        for column in range(width):
            pivot = next(
                (row for row in range(pivot_row, len(work))
                 if work[row][column]),
                None,
            )
            if pivot is None:
                continue
            work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
            scale = work[pivot_row][column]
            work[pivot_row] = [entry / scale for entry in work[pivot_row]]
            for row in range(len(work)):
                if row != pivot_row and work[row][column]:
                    scale = work[row][column]
                    work[row] = [
                        left - scale * right
                        for left, right in zip(work[row], work[pivot_row])
                    ]
            pivot_row += 1
            if pivot_row == len(work):
                break
        return pivot_row

    selected = []
    rows = list(rational_rows)
    current_rank = rational_rank(rows)
    for column in range(width):
        candidate = [
            Fraction(int(index == column), 1) for index in range(width)]
        next_rank = rational_rank(rows + [candidate])
        if next_rank > current_rank:
            rows.append(candidate)
            selected.append(column)
            current_rank = next_rank
        if current_rank == width:
            break
    if current_rank != width:
        raise AssertionError("failed to construct a complete complement")
    return selected


def rational_target_rows(degree_map):
    """Decode reconstructed target rows from a map artifact."""
    rows = []
    for target in STRESS_TARGETS[int(degree_map["degree"])]:
        rows.append([
            Fraction(
                entry["numerator"],
                entry["denominator"],
            )
            for entry in degree_map["targets"][target]
        ])
    return rows
