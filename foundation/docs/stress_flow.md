# Exact low-degree stress-flow map

## Scope and result

This branch starts the physics stage without changing the verified degree-10
branch.  It constructs the exact map from scalar invariants of the **free**
chiral four-form stress tensor to the committed five-form basis through
degree 10, reproduces the D=10 ModMax calculation of
Hutomo--Lechner--Sorokin, and performs a precisely scoped perturbative closure
pilot.

The homogeneous dimensions are:

| five-form degree | full verified value space | free-stress-generated space | complement |
|---:|---:|---:|---:|
| 4 | 1 | 1 | 0 |
| 6 | 2 | 1 | 1 |
| 8 | 7 | 2 | 5 |
| 10 | 14 | 2 | 12 |

The first nonzero complement occurs at degree 6.  In the committed basis it
is represented by `I6_2`; the stress direction is
`Tr(M^3)=(32/3) I6_1`.

These are statements about the verified homogeneous value spaces and the
free-seed polynomial stress algebra.  They are not a claim that the 81
functionally independent invariants generate the full polynomial ring, nor a
proof of nonlinear all-coupling flow closure.

## Paper conventions and normalization

The implementation follows arXiv:2509.14351v2:

- \(\eta_{\mu\nu}=\operatorname{diag}(-1,+1,\ldots,+1)\);
- \(\varepsilon^{01\ldots9}=-1\);
- \(\Lambda_5={*}\Lambda_5\);
- \(M_\mu{}^\nu=\Lambda_{\mu\rho_1\ldots\rho_4}
  \Lambda^{\nu\rho_1\ldots\rho_4}\);
- \(N_{\mu_1\mu_2\mu_3,\nu_1\nu_2\nu_3}
  =\Lambda_{\mu_1\mu_2\mu_3\rho\sigma}
  \Lambda_{\nu_1\nu_2\nu_3}{}^{\rho\sigma}\);
- \(T^{\rm free}_{\mu\nu}=M_{\mu\nu}/(2\cdot4!)=M_{\mu\nu}/48\).

The full six-index 1,050 and 4,125 components are implemented directly from
equations (2.15) and (2.17).  Tests check the 4,125 block symmetries, vanishing
traces, four-index antisymmetrization identity, and its agreement with the
paper's relevant symmetric two-index contraction
\((N^{(4125)}MM)_{\mu\nu}\).  The faster production stress path obtains that
contraction through the equivalent equations (2.17), (3.3), and (3.4):

\[
R_{\mu\nu}
=2(M^3)_{\mu\nu}-\frac12\eta_{\mu\nu}\operatorname{Tr}(M^3)
-12(NMM)_{\mu\nu}-\frac9{14}I_4M_{\mu\nu},
\]

\[
(N^{(4125)}MM)_{\mu\nu}
=\frac1{12}\left[
\frac57\left((M^3)_{\mu\nu}
-\frac1{10}\eta_{\mu\nu}\operatorname{Tr}(M^3)\right)-R_{\mu\nu}
\right].
\]

Both tensors are checked to be symmetric and traceless over three independent
prime fields.  The optimized implementation agrees with a direct
two-operand reference implementation.

For \(\mathcal V=\mathcal V(I_4)\), the implemented tensor is

\[
T_{\mu\nu}
=\frac1{48}M_{\mu\nu}
\left(1-\frac{96}{7}\mathcal V_I^2 I_4\right)
-\frac1{24}\eta_{\mu\nu}(\mathcal V-2I_4\mathcal V_I)
+\frac{\mathcal V_I^2}{3}R_{\mu\nu}.
\]

This is tested exactly against the unprojected-\(N\) expression in equation
(3.3).

Equation (3.16) uses the \(4!\)-rescaled INZ density: the action in equation
(2.7) carries an overall \(1/4!\).  This accounts for the apparent factor of
24 if the action density and the rescaled flow density are compared without
adjustment.

## ModMax reproduction

For \(\mathcal V=b\sqrt{I_4}\), the square root cancels from the stress tensor:

\[
T_{\mu\nu}
=\frac1{48}\left(1-\frac{24b^2}{7}\right)M_{\mu\nu}
+\frac{b^2}{12I_4}R_{\mu\nu}.
\]

The displayed paper invariants are independently reproduced in both compact
and expanded form:

\[
I_8=M^{\mu\nu}R_{\mu\nu},\qquad
I_{12}=R^{\mu\nu}R_{\mu\nu}.
\]

Thus

\[
T_{\mu\nu}T^{\mu\nu}
=\frac{I_4}{4(4!)^2}\left(1-\frac{24b^2}{7}\right)^2
+\frac{b^2}{288I_4}\left(1-\frac{24b^2}{7}\right)I_8
+\frac{b^4}{144I_4^2}I_{12}.
\]

The compact and expanded equations (3.14)--(3.15), the direct square of
\(T\), and the right side above agree at three deterministic samples under
each of three primes.  With
\(b=-\frac12\sqrt{7/6}\tanh(\gamma/2)\), subtracting the \(I_8\) and \(I_{12}\)
terms gives exactly the modified root-flow equation (3.16).

In the committed degree-8 basis, the paper's structure is

\[
I_8^{\rm paper}
=2I_{8,1}-12I_{8,2}-\frac{18}{7}I_{4,1}^2.
\]

The \(I_{8,2}\) coefficient explicitly shows why the correction is not in the
two-dimensional free-stress subspace
\(\operatorname{span}\{\operatorname{Tr}(M^4),(\operatorname{Tr}M^2)^2\}\).

## Exact change of basis

The unnormalized maps are:

\[
\operatorname{Tr}(M^2)=2I_{4,1},
\qquad
\operatorname{Tr}(M^3)=\frac{32}{3}I_{6,1},
\]

\[
\operatorname{Tr}(M^4)=I_{8,1},
\qquad
(\operatorname{Tr}M^2)^2=4I_{4,1}^2,
\]

\[
\operatorname{Tr}(M^2)\operatorname{Tr}(M^3)
=\frac{64}{3}I_{4,1}I_{6,1}.
\]

For degree 10, in the ordered basis

```text
I10_1, ..., I10_12, I4_1*I6_1, I4_1*I6_2
```

the row for \(\operatorname{Tr}(M^5)\) is

```text
1/27097 * [
  289944576, -15485184, -65968128, -11038464,
  55185408, 0, 0, -23929344, -42799104,
  336752640, 21378816, 0, 1124000, -433728
]
```

where \(48107520/3871=336752640/27097\).  The machine artifact stores every
fraction in reduced form and includes its residue under every prime.

Physical free-stress rows follow by dividing a degree-\(2k\) trace by
\(48^k\).  The implementation exposes
\(\operatorname{Tr}(T^2),\ldots,\operatorname{Tr}(T^{10})\), the
Cayley--Hamilton limit for a \(10\times10\) matrix.  Only powers two through
five have leading five-form degree at most 10; powers six through ten begin
at degrees 12 through 20.

## Exact complements and obstruction certificates

The chosen exact complement bases are:

- degree 6: `I6_2`;
- degree 8: `I8_2`, `I8_3`, `I8_4`, `I8_5`, `I8_6`;
- degree 10: `I10_1`, ..., `I10_12`.

This choice is a direct-sum complement, not an assertion that a stress row has
zero coefficients on those graph representatives.  At degree 10,
\(\operatorname{Tr}(M^5)\) necessarily mixes primitive and product
coordinates.

For every listed complement element, the artifact saves the exact value
matrix at 18 deterministic self-dual samples.  Each change-of-basis row is
solved on a deterministic full-rank subset and checked on the remaining
held-out samples.  Under each of five primes the artifact also records

```text
rank(stress values) = r
rank(stress values | obstruction value) = r + 1
```

No non-membership conclusion is based on a failed fit, SVD, or floating-point
tolerance.

## Flow-closure pilot

For a general linearized deformation

\[
\partial_\lambda\mathcal V
=\sum_{d=4,6,8,10}\sum_a c_{d,a}I_{d,a},
\]

free-seed polynomial stress closure is equivalent, degree by degree, to the
coefficient vector lying in the saved stress-row span.  This gives a
nontrivial family with dimensions \(1,1,2,2\) at degrees \(4,6,8,10\).
Closure removes all quotient directions but does not uniquely fix the
remaining stress parameters.  Its first obstruction is `I6_2`.

For the nonlinear one-coupling ansatz \(\mathcal V=c_4I_4\), the exact
homogeneous tensor is

\[
T=\frac{M}{48}+\frac{c_4I_4}{24}\mathbf1
+c_4^2\left(-\frac27I_4M+\frac13R\right).
\]

Direct polynomial matrix multiplication verifies every contribution through
five-form degree 10.  For example,

\[
\left.\operatorname{Tr}(T^2)\right|_{d=8}
=c_4^2\left(\frac{I_8}{72}+\frac{11I_4^2}{2016}\right),
\]

\[
\left.\operatorname{Tr}(T^3)\right|_{d=10}
=\frac{c_4^2}{48^2}
\left(\operatorname{Tr}(M^2R)-\frac67I_4\operatorname{Tr}(M^3)\right).
\]

The exact basis coordinates for these terms are in the machine artifact.
Mixed derivative-square contributions for a completely general
multi-invariant nonlinear ansatz are not classified here; that is explicitly
left unresolved rather than inferred from the linearized calculation.

## Degree-12 import interface

The ten degree-12 product directions are already registered:

- six `I4_1*I8_i` products;
- `I6_1^2`, `I6_1*I6_2`, `I6_2^2`;
- `I4_1^3`.

The registry requires exactly 62 concrete primitive directions and uses
unambiguous upper-triangle graph records, so order-12 labels do not depend on
the compact one-digit vertex format.  Importing them will expose a
72-dimensional homogeneous degree-12 basis without changing the evaluator.
The paper's rational \(I_{12}=R^2\) structure is already implemented and
tested, but it is not mislabeled as one of those polynomial primitive slots.

## Reproduction

```bash
.venv/bin/python -m pytest tests -q

.venv/bin/python scripts/stress_flow_pipeline.py \
  --out results/stress_flow_exact_low_degree.json
```

The generated artifact is
`results/stress_flow_exact_low_degree.json`.  It contains all five primes,
all 18 sample seeds, modular value matrices, modular solutions, reconstructed
rational matrices, complement certificates, the closure pilot, and the
degree-12 import schema.  Its arXiv v2 source and PDF checksums are pinned
constants, so regeneration is byte-for-byte independent of ignored local
paper downloads.  The CLI rejects a CRT modulus too small for the verified
degree-10 coefficient bound and requires at least one held-out sample beyond
the 14 degree-10 fit rows.  Custom moduli are restricted to the backend's
small-prime exact-arithmetic range and may not divide fixed map denominators.
