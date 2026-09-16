from chiral4form.degree12_exact_targets import (
    B13_SAMPLE_BOUND,
    FULL_PRIMES,
    GRAPH12_SAMPLE_BOUND,
    TRACE_PAIR_BOUND,
)
from chiral4form.degree12_exact_source_rows import VAR_INDEX


def test_exact_interpolation_bounds_are_conservative_and_named_by_scope():
    assert GRAPH12_SAMPLE_BOUND == 10**30
    assert B13_SAMPLE_BOUND == 10**30
    assert TRACE_PAIR_BOUND > B13_SAMPLE_BOUND


def test_full_checkpoint_prime_set_is_distinct():
    assert len(FULL_PRIMES) == 8
    assert len(set(FULL_PRIMES)) == 8


def test_reduced_variable_order():
    assert VAR_INDEX == {"A": 0, "B": 1, "C": 2, "D": 3, "U": 4, "V": 5}
