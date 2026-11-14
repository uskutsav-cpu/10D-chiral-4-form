# Stress-flow definition and convention audit

## Scope and sources

This audit fixes the meaning of every stress-flow statement in this
repository.  It distinguishes the public baseline commit `04175461f...` from
the stronger implementation on `stress-flow/classification-through-degree12`.
The primary sources are:

- Hutomo, Lechner and Sorokin, [*On non-linear chiral 4-form theories in
  D=10*](https://arxiv.org/abs/2509.14351v2), especially equations
  (2.18)--(2.36) and (3.1)--(3.16);
- Ferko, Kuzenko, Lechner, Sorokin and Tartaglino-Mazzucchelli,
  [*Interacting Chiral Form Field Theories and TTbar-like Flows in Six and
  Higher Dimensions*](https://arxiv.org/abs/2402.06947), especially the
  distinctions among INZ, PST, clone, Hamiltonian and Lagrangian flows;
- Ferko, Kuzenko, Smith and Tartaglino-Mazzucchelli,
  [*Duality-Invariant Non-linear Electrodynamics and Stress Tensor
  Flows*](https://arxiv.org/abs/2309.04253), for the lower-dimensional
  all-seed result that motivates, but does not prove, its D=10 analogue;
- Kuzenko, [*Interacting p-form gauge theories: New
  developments*](https://arxiv.org/abs/2504.01421), as a review of the
  formulation dictionary.

The local pinned v2 source and PDF checksums are recorded in
`scripts/stress_flow_pipeline.py`.  Formula claims below are checked against
that pinned version, not against an unversioned transcription.

## Fields, signature and degree

The code uses

\[
\eta_{\mu\nu}=\operatorname{diag}(-1,+1,\ldots,+1),\qquad
\varepsilon^{01\ldots9}=-1,\qquad *\Lambda=\Lambda .
\]

`Lambda` is the auxiliary self-dual five-form of the INZ formulation.  It is
not an arbitrary physical five-form.  On the physical gauge-field equation
of motion, the paper identifies it with

\[
\Lambda_5=-F_5^+ .
\]

All degree labels in the atlas mean homogeneous polynomial degree in the
components of `Lambda` (equivalently in `F5+` after the on-shell
identification).  Thus `Lambda` has degree 1, `M` and the free stress tensor
have degree 2, `I4` has degree 4, and `I6_1` and `I6_2` have degree 6.
Coupling order is a separate filtration and must not be confused with field
degree.

## Which energy-momentum tensor is meant?

The tensor in `src/sdinv/stress.py` follows the Hilbert tensor obtained from
the INZ action.  Its fully interacting auxiliary-field form is

\[
T_{\mu\nu}[\mathcal V]
=\frac1{2\cdot4!}\left(
\Lambda_\mu{}^{\rho(4)}\Lambda_{\nu\rho(4)}
-25\,\mathcal V_{\Lambda\,\mu}{}^{\rho(4)}
       \mathcal V_{\Lambda\,\nu\rho(4)}
\right)
-\frac1{4!}g_{\mu\nu}
\left(\mathcal V-\frac12\Lambda\cdot\mathcal V_\Lambda\right).
\]

Here
\(\mathcal V_{\Lambda\,\mu(5)}
=\partial\mathcal V/\partial\Lambda^{\mu(5)}\) is anti-self-dual.
The paper reaches this formula by using the non-dynamical `Lambda` equation
to eliminate `B` in favor of `Lambda`.  It is therefore:

- fully interacting in the arbitrary function `V(Lambda)`;
- on the auxiliary-field equation used in equation (2.32);
- manifestly Lorentz covariant and independent of the PST vector after that
  substitution;
- not, at this stage, a statement that the physical `A4` equation has been
  imposed;
- equal to an `F5+` expression after imposing `Lambda=-F5+`;
- expressible in terms of the full physical `F5` only after also using the
  nonlinear self-duality equation.

At public commit `04175461f...`, the code implemented the following slices:

1. `free_stress`: \(T_{\mu\nu}=M_{\mu\nu}/48\);
2. `stress_v_i4` and `stress_v_i4_raw`: the interacting tensor only for
   \(\mathcal V=\mathcal V(I_4)\), in the decomposed equations (3.4) and the
   unprojected equation (3.3);
3. `modmax_stress`: the conformal specialization
   \(\mathcal V=b\sqrt{I_4}\).

It did **not** accept the derivative of a general multi-invariant
interaction.  The classification branch adds `interacting_stress`,
`interaction_gradient_i4`, and exact graph/product derivatives.  Equation
(2.33) is calibrated by reproducing equation (3.3) over three prime fields.

## Which equations are imposed computationally?

Atlas and stress samples obey algebraic Hodge self-duality exactly:

\[
*\Lambda=\Lambda .
\]

They do not contain spacetime coordinates or derivatives.  Consequently the
sample calculations do not impose:

- the differential Bianchi identity \(dF=0\);
- the gauge-potential equation of motion;
- the nonlinear physical self-duality equation for a chosen `V`;
- stress-tensor conservation \(\partial^\mu T_{\mu\nu}=0\).

Conservation follows from diffeomorphism invariance on the appropriate full
equations of motion, but it cannot be tested on an isolated algebraic
five-form sample.  Repository tests therefore check covariance, symmetry,
trace, Hodge type and formula equivalence; they must not be cited as a direct
numerical conservation test.

## Stress scalars in the public low-degree map

The public static map used the **free**, traceless tensor

\[
T^{(0)}_{\mu\nu}=M_{\mu\nu}/48.
\]

It generated powers `Tr(T^k)` for \(k=2,\ldots,10\), subject to the
ten-dimensional Cayley--Hamilton limit.  `Tr(T)` was omitted because it
vanishes identically for this free self-dual tensor.  This omission is
harmless for that static free map but is not harmless for an interacting
classification: for a nonconformal interaction, `Tr(T[V])` contains direct
information about `V`.

The independent homogeneous free-stress monomials actually used through
degree 10 were:

| field degree | free-stress scalar rows |
|---:|---|
| 4 | \(\operatorname{Tr}T^2\) |
| 6 | \(\operatorname{Tr}T^3\) |
| 8 | \(\operatorname{Tr}T^4,\;(\operatorname{Tr}T^2)^2\) |
| 10 | \(\operatorname{Tr}T^5,\;\operatorname{Tr}T^2\operatorname{Tr}T^3\) |

Equivalently, before the factors of \(48^{-k}\), these are
\(\operatorname{Tr}M^2\), \(\operatorname{Tr}M^3\),
\(\operatorname{Tr}M^4\), \((\operatorname{Tr}M^2)^2\),
\(\operatorname{Tr}M^5\), and
\(\operatorname{Tr}M^2\operatorname{Tr}M^3\).

At degree 12 the corresponding free-stress candidates are

\[
\operatorname{Tr}T^6,\quad
\operatorname{Tr}T^4\operatorname{Tr}T^2,\quad
(\operatorname{Tr}T^3)^2,\quad
(\operatorname{Tr}T^2)^3 .
\]

Their exact rank and coordinates in the 72-dimensional degree-12 atlas are
not inherited from the degree-10 artifact and must be computed separately.

The classification branch changes `stress_traces` to expose
\(\operatorname{Tr}T,\ldots,\operatorname{Tr}T^{10}\).  This is a deliberate
interface extension; it does not retroactively change the public
free-stress map.

## Products and full invariant spaces

The homogeneous five-form bases include lower products, not only connected
graph representatives:

- degree 8: `I4_1^2`;
- degree 10: `I4_1*I6_1`, `I4_1*I6_2`;
- degree 12: `I4_1^3`, `I6_1^2`, `I6_1*I6_2`, `I6_2^2`, and
  `I4_1*I8_i` for \(i=1,\ldots,6\).

The verified dimensions \(1,2,7,14,72\) at degrees \(4,6,8,10,12\)
refer to these complete homogeneous value spaces.  The number 81 instead
refers to cumulative generic Jacobian rank through degree 12.  These two
notions are not interchangeable.

## Trace and conformal homogeneity

The first bracket in equation (2.33) is traceless for self-dual `Lambda` and
anti-self-dual `V_Lambda`.  Directly tracing the displayed equation (2.33)
with the stated metric gives

\[
\operatorname{Tr}T
=-\frac5{12}
\left(\mathcal V-\frac12\Lambda\cdot\mathcal V_\Lambda\right).
\]

The pinned v2 source prints a **plus** sign in equation (2.36).  This is
internally inconsistent with the minus sign multiplying \(g_{\mu\nu}\) in
equation (2.33), and with the independently reproduced \(V(I_4)\) formula,
which also carries the minus sign.  The zero-trace condition is unaffected:

\[
\Lambda\cdot\mathcal V_\Lambda=2\mathcal V .
\]

For a homogeneous polynomial component \(\mathcal V_d\) of field degree
\(d\), Euler's theorem gives

\[
\left.\operatorname{Tr}T\right|_{\mathcal V_d}
=-\frac5{24}(2-d)\mathcal V_d
=\frac5{24}(d-2)\mathcal V_d
\]

under the equation-(2.33) sign convention.  Therefore every nonzero analytic
polynomial interaction of degree \(d\ge4\) is nonconformal.  Conformal
examples such as \(\sqrt{I_4}\), \((I_6)^{1/3}\), or invariant ratios are
homogeneous of field degree 2 but are non-polynomial or singular in the
polynomial coordinates.  Analytic and conformal classifications must
therefore be kept separate.

## What the old obstruction statement did and did not prove

The exact public result

\[
\operatorname{Tr}M^3=\frac{32}{3}I_{6,1}
\]

and the rank certificate for `I6_2` prove that the degree-6 free-stress
static subspace is one-dimensional inside the two-dimensional sextic
invariant space.  They prove linearized non-membership of the second quotient
class for a polynomial flow generator evaluated on the free seed.

They do **not** prove that the class:

- is absent from `Tr(T[V])` for a seed already containing it;
- cannot be generated indirectly by higher interacting traces;
- is an obstruction for every admissible seed;
- survives nonanalytic conformal generators or invariant ratios;
- is formulation-independent without checking the permitted on-shell field
  redefinitions.

In particular, once `Tr(T)` is admitted, a nonconformal homogeneous seed
component is visible directly through the trace/homogeneity relation.  The
role of the sextic quotient must therefore be classified separately for a
free seed, a general seed, conformal flows, nonconformal analytic flows, and
broader nonanalytic generators.

## Audited status

| claim | evidence | status |
|---|---|---|
| public low-degree dimensions \(1,1,2,2\) | exact artifact and held-out rank certificates | verified for the free static map through degree 10 |
| `I6_2` outside the free sextic stress row | five-prime rank-increment certificates | verified, basis label still to be replaced by an intrinsic quotient formula |
| equations (3.3) and (3.4) agree | independent tensor paths, three primes | verified |
| equation (2.33) reduces to equation (3.3) | explicit anti-self-dual derivative, three primes | verified on this branch |
| ModMax \(T\), \(I_8\), \(I_{12}\), and \(T^2\) | compact/expanded/reference paths, three primes | verified algebraically |
| conservation | no spacetime-dependent implementation | not computationally tested |
| fully interacting closure through degree 12 | no reduced coefficient system yet | not established |
| universal `I6_2` obstruction | contradicted as an inference from the static map alone | not established |

