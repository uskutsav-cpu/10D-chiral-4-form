"""Interface for a future independently supplied spinor implementation.

No mentor or third-party code lives here. An adapter implements the protocol
below and returns one value matrix per degree on the same compact five-form
samples used by the trace backend. Column spaces are compared exactly over a
finite field, so different bases and normalizations are harmless.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np

from .modp import P, RankSieve


class SpinorInvariantBackend(Protocol):
    """Minimal interface expected from an external spinor implementation."""

    name: str
    attribution: str

    def evaluate_degree(self, five_form_components, degree, prime):
        """Return shape (samples, invariants_at_degree) values modulo prime."""
        ...


def _matrix(values, prime):
    matrix = np.asarray(values, dtype=np.int64) % prime
    if matrix.ndim != 2:
        raise ValueError("invariant values must be a 2-D matrix")
    return matrix


def exact_column_rank(values, prime=P):
    """Exact finite-field rank of the columns of a sample-value matrix."""
    matrix = _matrix(values, prime)
    sieve = RankSieve(matrix.shape[0], prime)
    for column in matrix.T:
        sieve.add(column)
    return sieve.rank


def compare_column_spaces(trace_values, spinor_values, prime=P):
    """Compare two invariant bases on corresponding samples over F_p.

    Equal column spaces need not have the same number or order of columns.
    Equality holds exactly when both ranks equal the rank of their union.
    """
    trace = _matrix(trace_values, prime)
    spinor = _matrix(spinor_values, prime)
    if trace.shape[0] != spinor.shape[0]:
        raise ValueError("trace and spinor matrices use different samples")
    trace_rank = exact_column_rank(trace, prime)
    spinor_rank = exact_column_rank(spinor, prime)
    union_rank = exact_column_rank(
        np.concatenate((trace, spinor), axis=1), prime)
    return {
        "prime": int(prime),
        "samples": int(trace.shape[0]),
        "trace_columns": int(trace.shape[1]),
        "spinor_columns": int(spinor.shape[1]),
        "trace_rank": trace_rank,
        "spinor_rank": spinor_rank,
        "union_rank": union_rank,
        "equal_column_spaces": (
            trace_rank == spinor_rank == union_rank),
    }
