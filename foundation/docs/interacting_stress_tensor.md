# The fully interacting stress tensor

Implementation: `src/sdinv/stress.py`, `src/sdinv/interaction.py`,
`src/sdinv/formal_flow.py`
Gates: `tests/test_stress_flow.py`, `tests/test_formal_flow.py`

## 1. Why the interacting tensor, not the free one

The objective is explicit that the classification must use the fully
interacting `T[V]`, not the free or `V(I4)`-specialised tensor. The
difference is not cosmetic: the free tensor is traceless, whereas the
interacting tensor has

    Tr(T)|_{V_d} = (5/24)(d - 2) V_d

which is nonzero for every `d >= 4`. That trace is exactly what transports
the sextic quotient class in `docs/stress_flow_classification.md`. A
free-tensor-only analysis would miss the entire mechanism and would wrongly
conclude that K6 is unreachable *and* unpreservable.

## 2. The anti-self-dual projection — a correctness fix, not a refinement

An ambient graph derivative is **not** automatically the derivative on the
constrained self-dual field space. The gradient must be projected onto the
appropriate anti-self-dual subspace before it means anything.

This error was hidden for a while by an accident: the `I4` derivative happens
already to lie in the correct subspace, so every `V(I4)` check passed while
the general case was wrong. Anything beyond `V(I4)` — which is the entire
point of this phase — needs the projection.

Guarded by `test_registry_gradients_are_anti_selfdual_and_obey_euler`, which
checks both that gradients are anti-self-dual and that they satisfy the Euler
relation. Euler contractions are unchanged by the projection, which is why
the bug was silent.

## 3. Normalisation

    tau = 48 * T

Arbitrary generator coefficient functions absorb the corresponding constant
powers of 48, so this rescaling is free. It is recorded explicitly in every
certificate under `normalization`, and each trace generator additionally
carries its own `physical_T_rescaling` string, e.g.

    product(Tr(tau^k)) = 48^3 * product(Tr(T^k))     for tr_tau*tr_tau2

so physical coefficients can be recovered without re-deriving the powers.
The static artifact likewise carries `physical_free_stress_coordinates`
alongside raw coordinates, related by `48^6` at degree 12.

**Any comparison with published formulas must fix this normalisation first.**

## 4. Trace generators

Eighteen generators are enumerated, being all products of
`Tr(tau^k)` with total leading field degree at most 12:

| leading degree | generators |
|---:|---|
| 4 | `tr_tau`, `tr_tau2` |
| 6 | `tr_tau3` |
| 8 | `tr_tau4`, `tr_tau^2`, `tr_tau*tr_tau2`, `tr_tau2^2` |
| 10 | `tr_tau5`, `tr_tau*tr_tau3`, `tr_tau2*tr_tau3` |
| 12 | `tr_tau6`, `tr_tau*tr_tau4`, `tr_tau2*tr_tau4`, `tr_tau3^2`, `tr_tau^3`, `tr_tau^2*tr_tau2`, `tr_tau*tr_tau2^2`, `tr_tau2^3` |

`leading_field_degree` is what makes the degree-wise analysis exhaustive
rather than a sample: at field degree 6 only the first three can contribute,
and all three produce rows, so the degree-6 enumeration is complete. This is
asserted by `test_degree6_generator_set_is_exhaustive`.

Note that `tr_tau` has leading field degree **4**, not 2 — the free stress
tensor is traceless, so the trace's leading contribution only appears once
the interaction is switched on.

## 5. Certificate structure

Each `results/stress_flow/certificates/interacting_degree12_<prime>.json`
records 192 target rows indexed by
`(generator, field_degree, coefficient_monomial)`, distributed as

| degree | 4 | 6 | 8 | 10 | 12 |
|---|---:|---:|---:|---:|---:|
| targets | 2 | 4 | 15 | 37 | 134 |

with a fit/holdout column split (72 fit, 8 holdout at degree 12) and
`all_holdouts_passed: true`. Roughly 530 s per prime.

## 6. Reproduce

    .venv/bin/python scripts/interacting_flow_degree12.py compute \
        --prime PRIME \
        --static-checkpoint work/static-degree12/PRIME.checkpoint.json \
        --checkpoint work/interacting-flow/PRIME.checkpoint.json \
        --out results/stress_flow/certificates/interacting_degree12_PRIME.json
