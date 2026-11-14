# Stress-flow classification through degree 12

Status: exact, multi-prime, reproducible. Every number below is computed over
finite fields with four deterministic samples per prime and agrees on every
prime tested. Nothing here relies on floating point.

Reproduce:

    .venv/bin/python scripts/stress_flow_closure.py
    .venv/bin/python -m pytest tests/test_k6_flow_role.py \
                              tests/test_stress_flow_closure.py -q

## 1. The two objects that must not be confused

**Static free-stress span.** The span of the scalars built from the free
traceless INZ stress tensor, degree by degree. Recorded in
`results/stress_flow/dimension_table.json`.

**Stress-flow closure.** The smallest subspace closed under

    dV/dlambda = f(T[V], lambda)

containing a given seed, where `f` runs over analytic functions of the
characteristic traces of the *fully interacting* `T[V]` with arbitrary
lambda-dependent coefficients.

These are different, and the difference is large:

| degree | full dimension | static span | free-seed closure |
|---:|---:|---:|---:|
| 4 | 1 | 1 | 1 |
| 6 | 2 | 1 | **1** |
| 8 | 7 | 2 | 3 |
| 10 | 14 | 2 | 11 |
| 12 | 72 | 4 | 67 |

At degrees 8, 10 and 12 the closure is strictly larger than the static span
(3 > 2, 11 > 2, 67 > 4). **Any no-go argument built on the static span alone
is invalid.** This is pinned by
`test_closure_strictly_exceeds_the_static_span`.

At the same time the closure is a *proper* subspace at every degree above 4,
so a genuine obstruction does exist — it is simply much smaller than the
static numbers suggest. Pinned by
`test_closure_is_a_proper_subspace_at_every_degree_above_four`.

## 2. The role of K6

Of the options posed in the objective, the supported statement is **C**:
*K6 can occur only if it is present in the seed.* Neither A nor B holds.

### The exact mechanism

At field degree six only three trace generators can contribute, namely those
with `leading_field_degree <= 6`: `tr_tau`, `tr_tau2`, `tr_tau3`. The
certificates enumerate every resulting row, so the degree-six analysis is
exhaustive rather than a sample. The four rows, in the ordered repository
basis `(I6_1, I6_2)`, identical over all four primes:

| target | coordinates | note |
|---|---|---|
| `tr_tau3 \| d=6 \| c=` | `[10927, 0]` | pure generator on the free seed |
| `tr_tau2 \| d=6 \| c=I4_1` | `[0, 0]` | vanishes identically |
| `tr_tau \| d=6 \| c=I6_1` | `[40, 0]` | Euler, reproduces J6 |
| `tr_tau \| d=6 \| c=I6_2` | `[0, 40]` | Euler, reproduces K6 |

Exactly one row carries a nonzero `I6_2` coordinate, and it is the row whose
coefficient monomial *is* `I6_2`. Its value 40 is `10*(6-2)`, matching the
recorded identity `Tr(tau)[V_d] = 10*(d-2)*V_d` — that is, `Tr(tau)` is the
homogeneity (Euler) operator and does nothing but rescale what is already
present.

### The consequence, in intrinsic coordinates

Write `q6` for the intrinsic quotient coordinate normalised by `q6([K6]) = 1`;
in repository coordinates `q6(c1 I6_1 + c2 I6_2) = (125/3) c2` (see
`docs/intrinsic_sextic_basis.md`). The degree-six component of the flow gives

    d q6 / d lambda = 10*(6-2) * a(lambda) * q6 = 40 * a(lambda) * q6

where `a(lambda)` is the coefficient of the `Tr(tau)` generator. This is a
linear **homogeneous** ODE, so

    q6(lambda) = q6(0) * exp( 40 * integral of a )

Because the exponential never vanishes:

- `q6(0) = 0` implies `q6(lambda) = 0` for all lambda — K6 is never created;
- `q6(0) != 0` implies `q6(lambda) != 0` for all lambda — K6 is never destroyed.

**K6 is transported, never created.** It is a conserved *flag*, not an
obstruction that some clever choice of `f` could switch on.

This is why option A is wrong: K6 need not vanish in every pure stress flow.
It vanishes in every pure stress flow *from a seed without it*, which is a
strictly weaker and quite different statement.

### Confirmation at the closure level

The first-order argument above is confirmed by the full fixed-point closure,
which converges after 2 sweeps:

- free seed `V = c*I4_1`: degree-six closure is **1** — K6 absent;
- seed with K6 added: degree-six closure is **2**, and the degree-12 closure
  rises by exactly one, 67 -> 68;
- seed with J6 added: closure unchanged at (1,1,3,11,67), because J6 already
  lies in the free-seed closure.

So seeding K6 is not inert — it enlarges the reachable set at degree 12 as
well — but nothing inside a K6-free family ever reaches it.

## 2b. The same result from the assembled equations (independent confirmation)

`results/stress_flow/interacting_flow_equations.json` is the **rationally
reconstructed** coefficient system for `dV/dlambda = f(T, lambda)`, fit on
primes 32749, 32719, 32693, 32771, 32713 and validated on the independent
holdout prime 32717, with `all_modular_and_rational_holdouts_passed: true`.

It reports the **new forcing** span — the directions the flow can genuinely
*create* — computed with `Tr(tau)` excluded, for exactly the reason
established above: that generator only propagates what is already present.

| degree | full | new forcing | quotient (the obstruction) |
|---:|---:|---:|---:|
| 4 | 1 | 1 | **0** |
| 6 | 2 | 1 | **1** |
| 8 | 7 | 3 | 4 |
| 10 | 14 | 5 | 9 |
| 12 | 72 | 21 | 51 |

Degree 4 closes completely. Every degree above it leaves a genuine
obstruction, and at degree 6 that obstruction is exactly one-dimensional —
the class `[K6]`.

This was produced by `assemble_interacting_stress_adapted.py` over rational
arithmetic; the closure table in §1 was produced by `stress_flow_closure.py`
over modular arithmetic with a different algorithm and a different notion
(reachability rather than forcing). They agree exactly at degrees 4, 6 and 8
(1, 1, 3). The degree-10 and degree-12 entries differ (11 vs 5, 67 vs 21)
precisely because the closure *includes* `Tr(tau)` propagation while the
forcing span excludes it — the two numbers answer different questions and
should not be quoted interchangeably.

## 3. Closed families

**The free-seed stress-closed family** has dimensions `(1, 1, 3, 11, 67)` at
degrees `(4, 6, 8, 10, 12)`. It is the smallest pure-stress-closed family
containing the free Lagrangian, and it is stable: adding J6 changes nothing.

**Conformal / homogeneous flows.** `Tr(tau)` acts as `10*(d-2)`, so it
annihilates degree 2 and is the exact obstruction to conformality at every
higher degree. A flow whose generator omits `Tr(tau)` cannot transport K6 at
all: with `a = 0` the ODE above gives `dq6/dlambda = 0` identically, so `q6`
is not merely preserved in sign but frozen.

**What is not yet settled.** Whether the ModMax-like model lies in the
free-seed family, and the complete solution of the degree-8/10/12 coefficient
equations, require the rational coefficient system rather than the modular
ranks used here. See `docs/assumptions_limitations_and_open_questions.md`.

## 4. Falsification gates in place

- `test_pure_trM3_generator_carries_no_k6` asserts `Tr(tau^3)` has zero K6
  coordinate **and** nonzero J6 coordinate, so the no-creation statement is
  substantive rather than vacuous.
- `test_degree6_generator_set_is_exhaustive` fails if any eligible generator
  produces no row, which would mean the degree-six enumeration was incomplete.
- `test_seeding_j6_changes_nothing` fails if the closure is accidentally
  sensitive to a redundant seed direction.
