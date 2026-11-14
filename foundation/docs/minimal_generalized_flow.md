# Minimal generalized flow

Question: if the pure stress flow does not close, what is the smallest set of
additional intrinsic five-form directions needed?

    dV/dlambda = f(T, S, lambda)

Computed exactly over finite fields, four samples per prime, agreeing on
every prime tested. Reproduce with `scripts/stress_flow_closure.py`; pinned
by `tests/test_minimal_generalized_flow.py`.

## 1. The gap to be closed

Free-seed closure against the full invariant space:

| degree | full | free-seed closure | gap |
|---:|---:|---:|---:|
| 4 | 1 | 1 | 0 |
| 6 | 2 | 1 | **1** |
| 8 | 7 | 3 | **4** |
| 10 | 14 | 11 | **3** |
| 12 | 72 | 67 | **5** |

Degree 4 already closes. Everything above it does not.

## 2. Adding directions one at a time

| added to the seed | closure | remaining gap |
|---|---|---|
| — (free) | (1, 1, 3, 11, 67) | (0, 1, 4, 3, 5) |
| K6 | (1, 2, 3, 11, 68) | (0, 0, 4, 3, 4) |
| K6, I8_3 | (1, 2, 4, 11, 68) | (0, 0, 3, 3, 4) |
| K6, I8_3, I8_4 | (1, 2, 5, 11, 68) | (0, 0, 2, 3, 4) |
| K6, I8_3, I8_4, I8_5 | (1, 2, 6, 11, 68) | (0, 0, 1, 3, 4) |
| K6, I8_3…I8_6 | (1, 2, **7**, 11, 68) | (0, 0, **0**, 3, 4) |

## 3. Minimality through degree 8

Five directions are needed to close degrees 6 and 8: `K6` and the four
degree-8 complement directions `I8_3, I8_4, I8_5, I8_6`.

Minimality is established by removal, not by assertion:

- **Each contributes exactly one dimension.** The closure at degree 8 climbs
  3 -> 4 -> 5 -> 6 -> 7 as the four are added one at a time, so no two of
  them are redundant with each other.
- **`I8_1`, `I8_2` and `I4_1^2` are inert.** Adding any of them changes
  nothing, because they already lie inside the 3-dimensional free-seed
  closure at degree 8. They are not part of a minimal set.
- **K6 is not substitutable.** Seeding `I8_3…I8_6` *without* K6 gives
  (1, **1**, 7, 11, **67**): degree 6 stays open and degree 12 stays at 67.
  So K6 contributes independently at two different degrees, and no amount of
  degree-8 data replaces it.

## 4. What is still not closed

Even with all five directions, degrees 10 and 12 retain gaps of **3** and
**4**. Closing those requires additional directions at degree 10 and above,
which this seeding analysis can identify but which have not been reduced to
an intrinsic tensorial form (the degree-12 complement is 68-dimensional and
carries only graph labels at present).

So the honest statement through degree 12 is:

> The pure stress flow from the free seed does not close. Through degree 8
> the minimal completion requires exactly five additional directions —
> K6 together with `I8_3, I8_4, I8_5, I8_6` — and each is provably
> non-redundant. Degrees 10 and 12 remain open by 3 and 4 dimensions
> respectively.

## 5. Caveat on what "adding S" means here

This analysis computes closure under **enlarging the seed**: allowing `V` to
contain a direction and asking what the flow then reaches. That is exactly
computable from the existing certificates.

It is *not* the same as adding `S` as an independent **generator** of `f`,
which would let the flow produce `S` times an arbitrary function directly.
That variant needs certificate rows that do not exist yet — the current
certificates enumerate trace generators only. The seeding result bounds the
generator result from one side (a generator can do at least what a seed can),
but the two are not proven equal here. This distinction is deliberately left
open rather than papered over.
