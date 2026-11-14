"""G-10: the free stress-tensor trace vanishes at quadratic order.

The closure's leading-degree bookkeeping rests on one physics statement:

    Tr(tau) begins at field degree 4, not 2, because the free stress tensor of a
    self-dual five-form is traceless.

Everything downstream depends on it. If Tr(tau) began at degree 2, every
generated target would land in the wrong graded piece and `dim D10 = 11` --
hence `dim Q10 = 3` -- would be wrong with it.

This file does not import the flow code, the closure, or the leading-degree
table. It re-derives the statement from the five-form and the Hodge star alone.

The argument
------------
The free stress tensor of a p-form is quadratic in F, so its trace is a degree-2
scalar. Whatever the improvement convention, every candidate trace is a multiple
of the single scalar

    <F, F> = F_{m1..m5} F^{m1..m5},

because the only quadratic Lorentz scalar available is that full contraction:
for the standard p-form stress tensor

    T_{mn} ~ F_{m...} F_n^{...} - c eta_{mn} <F,F>,

the trace is (1 - c d) <F,F>, a multiple of <F,F> for ANY c. So the vanishing is
independent of the improvement term and of the overall normalisation.

For a middle-degree form in even dimension, F ^ F is a top form built from two
copies of an odd-degree form, hence identically zero. Since <F, *F> is
proportional to F ^ F, it vanishes; and on either eigenspace of the star,
F = +-*F gives <F,F> = +-<F,*F> = 0.

The tests below check the computational half -- that <F,F> vanishes identically
on both eigenspaces and does NOT vanish generically -- which is what makes the
argument's conclusion true of this implementation and not only on paper.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sdinv.forms import (                       # noqa: E402
    hodge_matrix, random_form, to_dense, check_star_squared,
)
from sdinv.modp import inv                      # noqa: E402

PRIME = 32749
DIM, DEG = 10, 5
ETA = np.diag([-1] + [1] * (DIM - 1)).astype(np.int64)


def _projectors(p):
    H = hodge_matrix(DIM, DEG, True, p)
    n = H.shape[0]
    I = np.eye(n, dtype=np.int64)
    half = inv(2, p)
    return H, ((I + H) * half) % p, ((I - H) * half) % p


def _full_contraction(F, p):
    """<F,F> with all five indices raised by the Lorentzian metric."""
    G = F
    for ax in range(DEG):
        G = np.moveaxis(np.tensordot(ETA, G, axes=([1], [ax])), 0, ax)
    return int(np.tensordot(F, G, axes=DEG) % p)


def _sample(projector, p, seed):
    raw = random_form(DIM, DEG, np.random.default_rng(seed), p)
    v = raw if projector is None else (projector @ raw) % p
    return v, to_dense(v, DIM, DEG, p)


def test_star_squared_is_plus_one():
    """The whole self-duality setup needs *^2 = +1 on middle forms."""
    assert check_star_squared(DIM, DEG, True, PRIME) == 1


@pytest.mark.parametrize("seed", range(8))
def test_quadratic_trace_vanishes_on_self_dual(seed):
    H, P_sd, _ = _projectors(PRIME)
    v, F = _sample(P_sd, PRIME, 5000 + seed)
    assert np.array_equal((H @ v) % PRIME, v % PRIME), "sample is not self-dual"
    assert _full_contraction(F, PRIME) % PRIME == 0, (
        "the degree-2 part of the stress trace does not vanish on a self-dual "
        "five-form; the leading-degree bookkeeping and dim Q10 = 3 depend on it")


@pytest.mark.parametrize("seed", range(8))
def test_quadratic_trace_vanishes_on_anti_self_dual(seed):
    """The argument predicts it on BOTH eigenspaces; check the other one too."""
    H, _, P_asd = _projectors(PRIME)
    v, F = _sample(P_asd, PRIME, 5000 + seed)
    assert np.array_equal((H @ v) % PRIME, (-v) % PRIME)
    assert _full_contraction(F, PRIME) % PRIME == 0


@pytest.mark.parametrize("seed", range(8))
def test_it_does_not_vanish_generically(seed):
    """The control. Without it the two tests above prove nothing."""
    _, F = _sample(None, PRIME, 5000 + seed)
    assert _full_contraction(F, PRIME) % PRIME != 0, (
        "<F,F> vanished on an unprojected five-form, so the vanishing above is "
        "not evidence of self-duality doing any work")


@pytest.mark.parametrize("p", [32749, 32719, 32713, 32707])
def test_conclusion_holds_at_several_primes(p):
    H, P_sd, _ = _projectors(p)
    for seed in range(3):
        _, F = _sample(P_sd, p, 7000 + seed)
        assert _full_contraction(F, p) % p == 0


def test_leading_degree_table_agrees_with_the_derivation():
    """The table says 4. The derivation says 'not 2'. They must not diverge.

    This is the only place the production table is read, and it is read to be
    checked against an independent argument -- not to supply one.
    """
    from sdinv.formal_flow import TRACE_GENERATORS
    tr_tau = next(g for g in TRACE_GENERATORS if g.id == "tr_tau")
    assert tr_tau.leading_field_degree != 2, (
        "the table assigns Tr(tau) leading degree 2, contradicting the vanishing "
        "of the quadratic trace on self-dual five-forms")
    assert tr_tau.leading_field_degree == 4, (
        f"the table assigns Tr(tau) leading degree "
        f"{tr_tau.leading_field_degree}; the derivation gives 4")
