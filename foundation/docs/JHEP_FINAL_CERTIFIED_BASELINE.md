# Certified baseline — frozen

Machine-readable: `results/jhep/final_certified_baseline.json`, which re-derives
every value below from the certificates and records 18 checks, **0 failed**.

## Rank 81

| requirement | value |
|---|---|
| candidates scheduled / evaluated | 83 / 83 |
| evaluation errors | 0 |
| interrupted | 0 |
| structurally rejected (silently skipped) | 0 |
| Jacobian zero rows | 0 |
| Euler homogeneity | 83 / 83, none failed |
| exact modular rank | 81 |
| certified sample points | 2 |
| explicit minor | 81 x 81, determinant nonzero |
| independent determinant routines | agree |

**The two bounds, kept apart.**

*Lower.* The coordinate basis is integral, so the Jacobian is the reduction of an
integer matrix and `rank_{F_p} <= rank_Q`. The observed modular rank is therefore
an **unconditional** characteristic-zero lower bound: `rank_Q >= 81`. The
explicit non-vanishing minor witnesses this independently of the rank routine,
and a single prime suffices for it — an integer determinant that vanished would
vanish modulo every prime.

*Upper.* `126 - 45 = 81`, from the trivial generic stabiliser. This is
**analytic and from the literature**. It is not ours and is cited.

Together they pin the generic characteristic-zero rank at 81. The manuscript
states both explicitly, which is the condition under which it may say so.

**Four things this does not say**, and the manuscript does not:

- the rank at one modular sample is not the characteristic-zero rank;
- the characteristic-zero lower bound is not the generic rank;
- the generic rank is not established by any finite set of points;
- 83 candidates of rank 81 means the selection is functionally **dependent** —
  they are not algebraically independent and are never described as such.

## Bridge

Exact over `F_p` at both primes, no floating point anywhere:

- Clifford relations exact;
- the oscillator frame's real signature is `(5,5)` — **split**, not Euclidean —
  computed from all one hundred anticommutators of the archive's own operators;
- `*^2 = +1` on five-forms, so real self-dual five-forms exist there;
- gamma-traceless dimension 126; forward-map rank 126;
- kernel = the anti-self-dual 126 and image = the gamma-traceless 126, each by
  two-way span equality rather than dimension counting;
- left inverse composes to the self-dual projector exactly;
- equivariant under `GL(5)` with character `det(A)`, and under products of
  Clifford reflections, which generate the full orthogonal group.

Over `R`, `(5,5)` and `(1,9)` are inequivalent real forms; the transition is
exact over `F_p` and over `C`, and every component-level comparison is stated at
that level and no stronger.

## Degree 8

At four primes, on a hashed common registry of self-dual five-forms:

    dim A_8^tensor = dim A_8^spinor = dim(A_8^tensor + A_8^spinor) = 7

with two-way containment and a change-of-basis map validated on holdout samples.
The union rank is the load-bearing number: two *different* seven-dimensional
subspaces would also give equal ranks, and a union of eight.

Family ablation:

| family removed | remaining rank |
|---|---:|
| port graphs | 7 |
| products | 7 |
| structured | 7 |
| **tensor words** | **6** |

**Permitted wording, and used:** tensor words are indispensable *within the
tested candidate-family decomposition*. No universal uniqueness is claimed — a
family not tested here could in principle supply the same direction.

## What is frozen

These three results are not revisited by later phases. Any change to them would
require a demonstrated defect, and none is known.
