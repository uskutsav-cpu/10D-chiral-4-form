# Assumptions, limitations and open questions

This document exists so that no claim in this repository is read as stronger
than the evidence behind it. Anything listed here is **not** settled.

## 1. Assumptions behind the classification

1. **Analytic polynomial interactions.** `V` is expanded as a polynomial in
   the invariant basis through field degree 12. Non-polynomial interactions
   (`sqrt(I4)`, `(I6)^(1/3)`, invariant ratios) are outside the analysis.
   This matters because those are exactly where the conformal examples live.
2. **Scalar stress generators.** The flow generator `f` is a function of the
   characteristic traces `Tr(tau^k)` and their products. Non-scalar or
   derivative-dependent generators are not considered.
3. **Truncation at degree 12.** Everything is exact *through* degree 12.
   Nothing is claimed at degree 14 or beyond, and the pattern through 12 is
   not evidence about what happens above it.
4. **Fixed conventions.** `eta = (-,+,...,+)`, `tau = 48*T`, and the
   equation-(2.33) sign convention. Results are convention-dependent in
   normalisation, though the *quotient class* statements are not (§3 of
   `docs/intrinsic_sextic_basis.md`).

## 2. Known limitations

### 2.1 Tr(M^6) has no certified rational lift

With 15 primes (CRT modulus ≈ 5.2e67), **29 of 72 columns** fall outside the
uniqueness bound. The degree-12 stress row is therefore expressed in a
stress-adapted basis with `Tr(M^6)` as a basis element, not in the atlas
basis with rational coordinates. The *rank* statement (4-dimensional span) is
certified; the *coordinates* are not. See
`docs/stress_subalgebra_through_degree12.md` §3.

### 2.2 The interacting flow equations DO lift — with enough primes

This limitation was real at three primes and has since been removed, which is
worth recording because it calibrates how much the `Tr(M^6)` failure means.

With three fit primes, `assemble_interacting_stress_adapted.py` failed at
column 1 (`reconstructed fraction exceeds uniqueness bound`). With **five**
fit primes (32749, 32719, 32693, 32771, 32713) it succeeded, and the result
passed the independent holdout prime 32717 —
`all_modular_and_rational_holdouts_passed: true`. The artifact is
`results/stress_flow/interacting_flow_equations.json`.

So the flow coefficient system **is** available in exact rational form. The
contrast with `Tr(M^6)` (§2.1), which is still unlifted at fifteen primes, is
informative: the flow coefficients have modest height while the degree-12
stress row does not. That is evidence the `Tr(M^6)` obstruction is a genuine
height problem rather than a shortage of effort, and it strengthens the case
for finding an analytic identity instead of adding primes (open question 1).

### 2.3 Seeding is not the same as adding a generator

The minimal-flow analysis computes closure under enlarging the **seed**. The
objective's `f(T, S, lambda)` — adding `S` as an independent **generator** —
is a strictly different computation requiring certificate rows that do not
exist. The seeding result bounds it from one side only. See
`docs/minimal_generalized_flow.md` §5.

### 2.4 Degrees 10 and 12 are not closed

Even with the five-direction completion, gaps of 3 and 4 remain at degrees 10
and 12. The missing directions are identified only as coordinates in a graph
basis; they have not been reduced to intrinsic tensorial form.

### 2.5 The K6 ↔ spinor identification is ours

arXiv:2509.14350v2 proves `(Sigma_1, Sigma_2)` is a sextic basis but does not
publish the change of basis to `(Tr(M^3), K6)`. Our identification is not
independently confirmed against the source.

### 2.6 Source inconsistency at equation (2.36)

The pinned source prints a sign inconsistent with equation (2.33) and with
the reproduced `V(I4)` formula. We implement the (2.33)-derived convention
and document both. This should be confirmed with the authors rather than
treated as settled.

## 3. Open questions

1. **An analytic identity for `Tr(M^6)`** in terms of the atlas generators
   would sidestep the reconstruction failure entirely. This is the highest
   value open item — it converts an unbounded prime hunt into algebra.
2. **How many primes would actually suffice** for the `Tr(M^6)` lift? The
   observed heights suggest a substantially larger budget, but no bound has
   been computed. Computing that bound is cheap and would tell us whether
   brute force is even viable.
3. **Does the ModMax-like model lie in the free-seed closed family?** This
   needs the rational coefficient system, not just modular ranks.
4. **Is the degree-8 complement `{I8_3..I8_6}` intrinsically meaningful?**
   Four directions is suspiciously structured; there may be a tensorial
   characterisation as there was at degree 6.
5. **Behaviour above degree 12.** Whether the closure gap grows, stabilises,
   or closes is entirely unknown.
6. **Non-polynomial and conformal flows**, excluded by assumption 1, are
   where the physically interesting conformal models sit.

## 4. What a referee should check first

- That the degree-12 atlas fingerprint still reproduces
  (`26b61c44…`) — everything downstream depends on it.
- That `Tr(M^6)` coordinates are **not** quoted as rationals anywhere.
- That the static span (1,1,2,2,4) is never presented as the dynamical
  reachable set (1,1,3,11,67).
- That the K6 statement is quoted as *"transported, never created"* and not
  as *"K6 must vanish"*, which is false for a seeded flow.
