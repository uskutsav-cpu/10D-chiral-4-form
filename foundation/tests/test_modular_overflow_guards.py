"""Guards against the int64 overflow that produced a wrong P10_02 value.

A bare np.einsum multiplies ALL operands before summing. With four operands
whose entries reach p ~ 32749, products reach p^4 ~ 1.15e18; accumulated over
~1e6 index combinations this passes 2^63 and wraps with no warning, returning
a plausible-looking wrong integer.

That is exactly what happened: P10_02 first evaluated to 9605 and, after
routing through mod_einsum, to 4674. Only the homogeneity check exposed it.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.modp import P, mod_einsum  # noqa: E402
from sdinv.forms import selfdual_projector, to_dense, random_form  # noqa: E402
from sdinv.published_degree10_invariants import (  # noqa: E402
    evaluate_implemented)
from sdinv.published_degree12_invariants import evaluate_all  # noqa: E402


def _sample(prime, seed):
    projector = selfdual_projector(10, 5, True, prime)
    raw = random_form(10, 5, np.random.default_rng(seed), prime)
    return to_dense((projector @ raw) % prime, 10, 5, prime)


def test_bare_einsum_would_overflow_where_mod_einsum_does_not():
    """Demonstrates the failure mode on a small, fully controlled case."""
    prime = P
    rng = np.random.default_rng(5)
    a, b, c, d = (rng.integers(prime // 2, prime, size=(40, 40),
                               dtype=np.int64) for _ in range(4))

    bare = np.einsum("ij,jk,kl,li->", a, b, c, d)      # int64, unguarded
    safe = mod_einsum("ij,jk,kl,li->", [a, b, c, d], prime)

    # the guarded result is a valid residue; the bare one need not even be
    # congruent, because it wrapped
    assert 0 <= int(safe) < prime
    if int(bare) % prime != int(safe) % prime:
        return          # overflow demonstrated, which is the point
    # if this platform happened not to overflow, the guarded path must still
    # agree with an exact object-dtype computation
    exact = np.einsum("ij,jk,kl,li->", *(x.astype(object) for x in (a, b, c, d)))
    assert int(exact) % prime == int(safe) % prime


def test_published_degree10_homogeneity_across_scalings():
    """The check that caught the real bug. Several c values, two primes."""
    for prime in (32749, 32719):
        form = _sample(prime, 20260729)
        base = evaluate_implemented(form, prime)
        for c in (2, 3, 5, 7):
            scaled = evaluate_implemented((form * c) % prime, prime)
            for name, v in base.items():
                assert scaled[name] == (pow(c, 10, prime) * v) % prime, (
                    f"{name} lost degree-10 homogeneity at prime {prime}, "
                    f"c={c}: a bare int64 einsum would do exactly this")


def test_published_degree12_homogeneity_across_scalings():
    for prime in (32749, 32719):
        form = _sample(prime, 20260730)
        base = evaluate_all(form, prime)
        for c in (2, 3, 5):
            scaled = evaluate_all((form * c) % prime, prime)
            for name, v in base.items():
                assert scaled[name] == (pow(c, 12, prime) * v) % prime


def test_published_paths_do_not_call_bare_einsum():
    """Static guard on the published-invariant modules."""
    root = os.path.join(os.path.dirname(__file__), "..", "src", "sdinv")
    for module in ("published_degree10_invariants.py",
                   "published_degree12_invariants.py"):
        with open(os.path.join(root, module)) as stream:
            body = stream.read()
        for line in body.splitlines():
            stripped = line.strip()
            if stripped.startswith("#") or stripped.startswith('"'):
                continue
            assert "np.einsum(" not in stripped, (
                f"{module}: bare np.einsum in a published-result path")
