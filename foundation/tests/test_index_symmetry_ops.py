"""The nested-bracket engine, including the order-matters regression.

The source states that red-bracket operations act UPON the black-bracket
antisymmetrisations. If that ordering could be ignored without consequence the
convention would be vacuous, so one test constructs an explicit example where
the two orders disagree.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.index_symmetry_ops import (  # noqa: E402
    BLACK, RED, BracketOp, BracketProgram, apply_bracket, permutation_sign)

MOD = 32749


def _rand(shape, seed):
    rng = np.random.default_rng(seed)
    return rng.integers(0, MOD, size=shape, dtype=np.int64)


def test_permutation_sign_basic():
    assert permutation_sign([0, 1, 2]) == 1
    assert permutation_sign([1, 0, 2]) == -1
    assert permutation_sign([1, 2, 0]) == 1
    assert permutation_sign([2, 1, 0]) == -1


def test_antisymmetriser_is_idempotent_and_kills_symmetric_part():
    t = _rand((6, 6), 1)
    sym = (t + t.T) % MOD
    op = BracketOp("antisym", (0, 1))
    assert np.array_equal(apply_bracket(sym, op, MOD) % MOD,
                          np.zeros_like(sym))
    a = apply_bracket(t, op, MOD)
    assert np.array_equal(apply_bracket(a, op, MOD) % MOD, a % MOD)


def test_symmetriser_is_idempotent_and_kills_antisymmetric_part():
    t = _rand((6, 6), 2)
    anti = (t - t.T) % MOD
    op = BracketOp("sym", (0, 1))
    assert np.array_equal(apply_bracket(anti, op, MOD) % MOD,
                          np.zeros_like(anti))
    s = apply_bracket(t, op, MOD)
    assert np.array_equal(apply_bracket(s, op, MOD) % MOD, s % MOD)


def test_three_slot_antisymmetriser_matches_explicit_sum():
    t = _rand((4, 4, 4), 3)
    got = apply_bracket(t, BracketOp("antisym", (0, 1, 2)), MOD)
    want = np.zeros_like(t)
    for perm, sign in ((( 0,1,2),1), ((1,2,0),1), ((2,0,1),1),
                       ((1,0,2),-1), ((0,2,1),-1), ((2,1,0),-1)):
        want = (want + sign * np.transpose(t, perm)) % MOD
    want = (want * pow(6, MOD - 2, MOD)) % MOD
    assert np.array_equal(got % MOD, want % MOD)


def test_antisymmetrising_a_repeated_slot_is_rejected():
    try:
        BracketOp("antisym", (0, 0))
    except ValueError:
        return
    raise AssertionError("repeated slot must be rejected")


def test_unknown_kind_and_stage_are_rejected():
    for bad in (lambda: BracketOp("skew", (0, 1)),
                lambda: BracketOp("sym", (0, 1), stage="green")):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("invalid bracket accepted")


def test_normalisation_convention_is_explicit():
    t = _rand((5, 5), 7)
    norm = apply_bracket(t, BracketOp("antisym", (0, 1), normalized=True), MOD)
    raw = apply_bracket(t, BracketOp("antisym", (0, 1), normalized=False), MOD)
    assert np.array_equal((2 * norm) % MOD, raw % MOD)


def test_stage_order_black_before_red():
    prog = BracketProgram([
        BracketOp("sym", (0, 1), stage=RED, source="red"),
        BracketOp("antisym", (1, 2), stage=BLACK, source="black"),
    ])
    black, red = prog.stages()
    assert [o.source for o in black] == ["black"]
    assert [o.source for o in red] == ["red"]


def test_red_before_black_gives_a_different_answer():
    """The regression that makes the source convention non-vacuous.

    If this ever passes with equality, the ordering statement in eq (4.24)
    would carry no content and every formula built on it would be suspect.
    """
    t = _rand((4, 4, 4), 11)
    prog = BracketProgram([
        BracketOp("antisym", (0, 1), stage=BLACK),
        BracketOp("sym", (1, 2), stage=RED),
    ])
    correct = prog.apply(t, MOD, reverse_stages=False)
    swapped = prog.apply(t, MOD, reverse_stages=True)
    assert not np.array_equal(correct % MOD, swapped % MOD), (
        "black-then-red and red-then-black agreed, which would make the "
        "source's ordering statement vacuous")


def test_serialization_is_deterministic_and_records_stage():
    prog = BracketProgram([
        BracketOp("antisym", (0, 1, 2), stage=BLACK),
        BracketOp("sym", (3, 4), stage=RED, normalized=False),
    ])
    s = prog.serialize()
    assert s == prog.serialize()
    assert "black:antisym[0, 1, 2]/norm" in s
    assert "red:sym[3, 4]/raw" in s


def test_no_overflow_for_high_rank_tensors():
    """Every intermediate stays below the modulus, so int64 cannot wrap."""
    t = _rand((6,) * 5, 13)
    out = apply_bracket(t, BracketOp("antisym", (0, 1, 2, 3)), MOD)
    assert out.dtype == np.int64
    assert int(out.max()) < MOD and int(out.min()) >= 0
