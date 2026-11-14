# Modular contraction safety

## The failure that motivated this

`P10_02` was first implemented with a bare `np.einsum` over four int64
operands. Entries reach `p ~ 32749`, so products reach `p^4 ~ 1.15e18`;
accumulated over ~1e6 index combinations this passes `2^63` and **wraps with
no warning**.

    before (bare np.einsum):  P10_02 = 9605
    after  (mod_einsum):      P10_02 = 4674

The overflowed value was a perfectly plausible residue. Nothing about it looked
wrong. It was caught only by the homogeneity check `F -> cF ~ c^10`, which the
wrapped value failed.

This is the same trap the README documents as having produced three wrong
answers historically. It is easy to walk into and invisible without a
structural check.

## Audit of every bare `np.einsum` in `src/sdinv`

| location | operands | verdict |
|---|---|---|
| `modp.py:176,184,192` | 2 | **safe by construction** — these are inside `mod_einsum` itself, which contracts strictly pairwise with an explicit `(p-1)^2 * nterms < 2^63` guard and reduces after every step |
| `stress.py:111,132,174,250` | 2 | **safe and deliberate** — `backend="reference"` oracles. The code documents the bound: "At most 10^4 terms of size < mod^2 enter each output, safely below int64." Kept as independent oracles and must not be optimised away |
| `published_degree10_invariants.py` | — | **none**; routed through `mod_einsum`, enforced by a static test |
| `published_degree12_invariants.py` | — | **none** |
| `index_symmetry_ops.py` | — | **none**; permutation sums reduce mod `p` after every term, so intermediates never exceed `p` regardless of tensor rank |

No unsafe call remains in a scientific result path.

## Rules adopted

1. Any contraction with **three or more** operands of unbounded modular
   magnitude goes through `mod_einsum`.
2. Two-operand reference oracles may use `np.einsum` **only** with a written
   magnitude bound, as `stress.py` has.
3. Every published invariant carries a homogeneity test at several scalings —
   this is the check that actually catches wrapping.
4. Reference and optimised backends must agree where both exist.
5. A static test forbids `np.einsum` in the published-invariant modules.

## Why homogeneity is the right detector

An overflow produces a valid-looking residue, so no range check finds it. But
wrapping is not equivariant under `F -> cF`: the wrapped value cannot scale as
`c^degree`. Testing several `c` at two primes catches it immediately, which is
precisely what happened here.

Tests: `tests/test_modular_overflow_guards.py`.
