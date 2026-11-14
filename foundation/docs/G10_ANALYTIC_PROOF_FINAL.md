# G-10 — the quadratic stress trace vanishes

Certificate: `results/stress_flow/G10_publication_certificate.json`.
Verifier: `tests/test_G10_trace_activation.py` (30 tests, imports no flow code).
Counterfactual: `results/stress_flow/G10_counterfactual.json`.

## The statement

> `Tr(tau)` begins at field degree 4, not 2, because the free stress tensor of a
> self-dual five-form is traceless.

The closure's leading-degree bookkeeping rests on it, and through that so does
`dim_Q D10 = 11` and `dim_Q Q10 = 3`.

## The argument

**1. The trace is a degree-2 scalar.** The free stress tensor of a `p`-form is
quadratic in `F`.

**2. Every candidate trace is a multiple of one scalar.** For

    T_{mn} ~ F_{m...} F_n^{...} - c eta_{mn} <F,F>

the trace is `(1 - c d) <F,F>`. This holds for **any** improvement coefficient
`c`, so the conclusion does not depend on the improvement term or on the overall
normalisation — only on whether `<F,F>` vanishes. That is what makes the
statement robust across formulations rather than tied to one action.

**3. `F ^ F = 0` identically.** It is a top form built from two copies of an
odd-degree form in even dimension, so it is its own negative.

**4. Hence `<F,F> = 0` on either eigenspace.** `<F, *F>` is proportional to
`F ^ F` and so vanishes; with `F = ±*F`, `<F,F> = ±<F,*F> = 0`.

Therefore the degree-2 part of `Tr(tau)` vanishes identically and the first
possible contribution is degree 4.

## Verification, with a control

`tests/test_G10_trace_activation.py` re-derives this from the five-form and the
Hodge star alone. It does not import the flow code, the closure, or the
leading-degree table — it reads the table once, at the end, to check it against
the derivation rather than to take a value from it.

| case | `<F,F>` | duality verified |
|---|---|---|
| self-dual | **0** at every sample | `*F = F` |
| anti-self-dual | **0** at every sample | `*F = -F` |
| **generic (control)** | **non-zero** | — |

`*^2 = +1` is checked, and the conclusion is re-tested at primes
`32749, 32719, 32713, 32707`.

The control is the part that matters. Without it, "the contraction is zero"
would be equally consistent with a contraction routine that returns zero for
everything.

An earlier version of this test used the wrong control: it called
`selfdual_projector(10, 5, False, p)`, where the third argument selects
*Lorentzian versus Euclidean*, not self-dual versus anti-self-dual. The apparent
non-vanishing on "anti-self-dual" forms was an artefact of a Euclidean
projector, where `*^2 = -1` and `(1+*)/2` is not a projector at all. Corrected,
both eigenspaces vanish, exactly as the argument predicts.

## The counterfactual — the statement is load-bearing

Re-running the exact rational closure with `tr_tau` targets forced
unconditionally active:

| | degree 4 | 6 | 8 | 10 | `Q10` |
|---|---:|---:|---:|---:|---:|
| as derived | 1 | 2 | 7 | **11** | **3** |
| counterfactual | 1 | 2 | 7 | **14** | **0** |

The paper's central number moves from 3 to 0. G-10 is not a technicality.

The counterfactual also reproduces **exactly** the historical error in which the
unconditional span gave rank 14 and quotient 0. It is retained as a negative
regression fixture, so that error cannot return unnoticed.

## What remains for the mentor

The derivation above settles the mathematics given the conventions. What a
coauthor is still needed for:

- confirmation that the free theory and stress-tensor convention used here are
  the intended ones;
- confirmation that no formulation in use — PST, auxiliary-field, clone,
  Hamiltonian — changes the quadratic trace.

Step 2 narrows this considerably: since the vanishing is independent of the
improvement coefficient, a differing formulation would have to change the
*quadratic scalar available*, not merely the trace convention.

## Status

`CERTIFIED` — analytically derived, computationally verified with a control, and
shown load-bearing by counterfactual. Mentor confirmation of conventions remains
outstanding and is item **G-10** on the decision form.
