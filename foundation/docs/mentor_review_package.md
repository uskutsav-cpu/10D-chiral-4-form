# Mentor review package

One-page orientation for review. Every claim links to the artifact and the
test that guards it. Read
`docs/assumptions_limitations_and_open_questions.md` alongside this — several
results are deliberately narrower than they might first appear.

## 1. Headline results

### 1.1 The static free-stress subalgebra through degree 12

| degree | full | stress span | quotient |
|---:|---:|---:|---:|
| 4 | 1 | 1 | 0 |
| 6 | 2 | 1 | 1 |
| 8 | 7 | 2 | 5 |
| 10 | 14 | 2 | 12 |
| 12 | 72 | 4 | 68 |

Exact over 15 primes, 4 samples each, same standard complement at every
prime. Artifact `results/stress_flow/dimension_table.json`.

### 1.2 The dynamical closure is much larger than the static span

| degree | static span | free-seed closure | full |
|---:|---:|---:|---:|
| 4 | 1 | 1 | 1 |
| 6 | 1 | 1 | 2 |
| 8 | 2 | **3** | 7 |
| 10 | 2 | **11** | 14 |
| 12 | 4 | **67** | 72 |

**This is the single most important correction in this phase.** A no-go
argument built on the static span would be badly wrong at degrees 10 and 12
(2 vs 11, 4 vs 67). A genuine obstruction still exists — the closure is a
proper subspace everywhere above degree 4 — but it is far smaller than the
static numbers suggest.

### 1.3 The role of K6: transported, never created

Writing `q6` for the intrinsic quotient coordinate,

    d q6 / d lambda = 40 * a(lambda) * q6
    =>  q6(lambda) = q6(0) * exp( 40 * integral of a )

A **linear homogeneous** ODE. Therefore:

- `q6(0) = 0` ⟹ `q6 ≡ 0`: K6 is never generated;
- `q6(0) ≠ 0` ⟹ `q6 ≠ 0` for all lambda: K6 is never destroyed.

Of the options posed, **C** is correct: *K6 can occur only if it is present
in the seed.* Option A — *"K6 must vanish in every pure stress flow"* — is
**false**, and this distinction is the part most likely to be misquoted.

The mechanism is exhaustive, not sampled: only three trace generators have
`leading_field_degree <= 6`, all three produce rows, and exactly one row
carries a nonzero K6 coordinate — the one whose coefficient monomial is K6
itself, with coefficient `40 = 10*(6-2)`, i.e. pure homogeneity.

### 1.4 Minimal generalized flow through degree 8

Five directions, each provably non-redundant: **K6, I8_3, I8_4, I8_5, I8_6**.
Minimality shown by removal — dropping any one reopens a gap. `I8_1`, `I8_2`,
`I4_1^2` are inert (already inside the closure). Degrees 10 and 12 remain
open by 3 and 4.

### 1.5 The exact flow coefficient equations

`results/stress_flow/interacting_flow_equations.json` gives the **rationally
reconstructed** system for `dV/dlambda = f(T, lambda)` through degree 12, fit
on five primes and validated on an independent holdout
(`all_modular_and_rational_holdouts_passed: true`). Its **new forcing** span
— what the flow can create, with the homogeneity generator `Tr(tau)` excluded
— is:

| degree | full | new forcing | obstruction |
|---:|---:|---:|---:|
| 4 | 1 | 1 | 0 |
| 6 | 2 | 1 | **1** |
| 8 | 7 | 3 | 4 |
| 10 | 14 | 5 | 9 |
| 12 | 72 | 21 | 51 |

This is an **independent confirmation** of §1.3: a different script, a
different algorithm, rational rather than modular arithmetic, and the same
one-dimensional sextic obstruction.

## 2. The negative result that matters

**`Tr(M^6)` has no certified rational lift.** With 15 primes (CRT modulus
≈ 5.2e67), 29 of 72 columns exceed the uniqueness bound. Coefficients were
**not** guessed; `Tr(M^6)` was adopted as a stress-adapted basis element
instead.

The contrast with §1.5 is the informative part: the *flow* coefficients
lifted cleanly at five primes, while the degree-12 stress row is still
unlifted at fifteen. That points to a genuine height problem specific to
`Tr(M^6)`, not to insufficient effort — which is why the recommended fix is
an analytic identity rather than more primes.

Please do not let a write-up quote rational `Tr(M^6)` coordinates. The
4-dimensional *rank* is certified; the coordinates are not.

## 3. Integrity of the published foundation

The degree-12 atlas is untouched by this phase:

- all 8 original `src/sdinv/*.py` files byte-identical to public `de696ca`;
- `results/10d_order12.json`, `results/degree12_benchmarks.json` unchanged;
- diff from `de696ca` is **1702 insertions, 0 deletions**, entirely in new
  modules;
- semantic fingerprint `26b61c44…5c25` reproduces.

## 4. Test inventory

| area | file |
|---|---|
| atlas order and certificate | `tests/test_degree12.py` |
| intrinsic sextic J6/K6 | `tests/test_sextic.py` |
| interacting stress, published formulas | `tests/test_stress_flow.py` |
| formal flow algebra | `tests/test_formal_flow.py` |
| static degree-12 certificates | `tests/test_static_stress_degree12.py` |
| **K6 dynamical role** | `tests/test_k6_flow_role.py` |
| **flow closure** | `tests/test_stress_flow_closure.py` |
| **minimality by removal** | `tests/test_minimal_generalized_flow.py` |

Falsification gates are deliberate: `test_pure_trM3_generator_carries_no_k6`
requires `Tr(tau^3)` to be nonzero at degree 6, so the no-creation claim is
substantive rather than vacuous; `test_degree6_generator_set_is_exhaustive`
fails if any eligible generator produced no row.

## 5. Questions for the mentor

1. Is there a known analytic identity for `Tr(M^6)` in the atlas generators?
   That single identity removes the largest limitation in this phase.
2. Equation (2.36) of arXiv:2509.14351v2 prints a trace sign inconsistent
   with (2.33) and with the reproduced `V(I4)` formula. Which is intended?
3. Is the degree-8 complement `{I8_3..I8_6}` known intrinsically, as `K6` is
   at degree 6?
4. Does the published change of basis from `(Sigma_1, Sigma_2)` to
   `(Tr(M^3), K6)` exist anywhere? Ours is inferred, not sourced.
