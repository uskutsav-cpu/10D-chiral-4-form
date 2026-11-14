# Reproduction of the published D=10 formulas

Every equation below is reproduced in exact finite-field arithmetic and
guarded by a test. Sources are pinned to specific arXiv versions; see
`docs/stress_flow_definition_audit.md` for the full convention audit.

Primary source: Hutomo, Lechner, Sorokin, *On non-linear chiral 4-form
theories in D=10*, arXiv:2509.14351v2, equations (2.18)–(2.36), (3.1)–(3.16).

## 1. What is reproduced, and by which test

| published object | test |
|---|---|
| M and N conventions, the 4125 contraction | `test_paper_m_n_and_relevant_4125_contraction_conventions` |
| full 4125 projector: symmetries, traces, contraction | `test_full_4125_projector_symmetries_traces_and_relevant_contraction` |
| eq (3.3) ≡ eq (3.4) | `test_equations_3_3_and_3_4_are_exactly_equivalent` |
| general interacting T reduces to the V(I4) formula | `test_general_interacting_stress_reduces_to_v_i4_formula` |
| gradient normalisation | `test_registry_gradient_has_paper_normalization` |
| published I8 / I12 structures, ModMax stress square | `test_published_i8_i12_and_modmax_stress_square` |
| eq (3.16) normalisation identity | `test_equation_3_16_normalization_identity` |
| Lorentz-boost covariance of M, R and the ModMax T | `test_lorentz_boost_covariance_of_m_r_and_modmax_t` |
| Cayley–Hamilton interface for stress traces | `test_stress_trace_cayley_hamilton_interface` |
| V = c·I4 trace expansion through degree 10 | `test_v_equals_c_i4_trace_expansion_through_degree10` |
| reference vs optimized backend agreement | `test_simple_reference_m_and_n_match_optimized_backend` |

Both a reference and an optimized backend are implemented and required to
agree, so an optimisation bug cannot hide behind a single implementation.

## 2. The equation (2.36) sign discrepancy — preserved, not smoothed

This is a genuine internal inconsistency in the pinned source and is recorded
rather than silently corrected.

Equation (2.33), with `eta = (-,+,...,+)`, has a traceless first bracket and a
`-g_{mu nu}/4!` term. Taking the trace algebraically gives

    Tr(T) = -(5/12) * ( V - (Lambda . dV)/2 )

**Equation (2.36) of arXiv:2509.14351v2 displays the opposite overall sign.**

Evidence for the minus sign being the one that follows from the derivation:

1. it is what (2.33) yields when the trace is taken directly;
2. the independently reproduced `V(I4)` formula carries the same minus sign.

The repository implements the (2.33)-derived convention and says so in
`src/sdinv/stress.py::interacting_trace_formula`, which names the discrepancy
in its docstring instead of quietly flipping a sign to force agreement.

**The zero-trace (conformal) condition is unaffected**, since it reads

    Lambda . dV = 2 V

either way.

## 3. Homogeneity, and an independent cross-check

For a homogeneous component `V_d` of field degree `d`, Euler's theorem plus
the (2.33) convention gives

    Tr(T)|_{V_d} = (5/24) * (d - 2) * V_d

The interacting certificates are computed in the rescaled normalisation
`tau = 48*T` and independently record

    Tr(tau)[V_d] = 10 * (d - 2) * V_d

These agree exactly: `48 * 5/24 = 10`. Two separate parts of the codebase —
the analytic trace formula and the numerically reduced certificate rows —
land on the same coefficient. At degree 6 it gives 40, which is precisely the
value observed in the degree-six flow rows underpinning the K6 result in
`docs/stress_flow_classification.md`.

## 4. Consequence for conformality

Since `Tr(T)|_{V_d}` vanishes only at `d = 2`, **every nonzero analytic
polynomial interaction of degree `d >= 4` is nonconformal.** Conformal
examples such as `sqrt(I4)`, `(I6)^(1/3)` or invariant ratios are homogeneous
of field degree 2 but are non-polynomial or singular in polynomial
coordinates. The analytic and conformal classifications must therefore be
kept separate, and this repository classifies the analytic case.
