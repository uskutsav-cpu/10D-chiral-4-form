# Bounded integer-residual certificate for the degree-eight map

## Why this is stronger than a held-out-prime fit

Ordinary interpolation at a finite set of tensors does not prove that an
unknown polynomial identity holds. CRT reconstruction also does not fix that
without a height/completeness argument. This certificate instead uses two
separate finite, checkable steps.

### Step 1: rational residuals vanish at fixed bounded tensors

Nine inputs are fixed by 126 electric components in `{-1,0,1}`. The magnetic
components are determined by real Lorentzian self-duality, and all 252
canonical components still have absolute value at most one. The **same
integer tensors** are evaluated modulo seven distinct primes:

```
30011, 30013, 30029, 30047, 30059, 30071, 30089.
```

A degree-n graph invariant has `5n/2` contracted index labels and hence value
bounded by `10^(5n/2)` on those tensors. A derivative at a specified dense
slot has at most `n * 10^(5(n-2)/2)` terms. Antisymmetrization and Hodge
projection are averages, so do not increase that component bound.

A constrained derivative has denominator dividing `2*5! = 240`. In every
stress-generator coefficient needed through degree eight, derivatives occur
at most quadratically, and only the quartic derivative is needed outside the
exact trace identity. Thus `D=240^2=57600` clears all target denominators.
It also clears the denominators of the proposed rational basis-map rows.

For example, write `G4=10^10`, `G8=10^20`, `Mmax=10^4`, and
`D4max=4*10^5`. The bilinear derivative bound is
`B44max=10^4*D4max^2`. The largest conservative target bound is

```
2*10^2*Mmax*25*B44max + 10*(2*G4)^2 = 8.4 * 10^22.
```

Every proposed map row has coefficient absolute sum at most 400, so its
value is bounded by `400*G8=4*10^22`. The absolute integer residual is
therefore bounded by

```
B = 57600 * (8.4*10^22 + 4*10^22).
```

The prime product exceeds `2B`. If every modular residual is zero, its
integer value must be zero, hence the rational residual is exactly zero at
that tensor. The code records the complete bound arithmetic, the fixed
inputs, every target/candidate value, all residues, and the modulus product.
No assumption about unknown reconstructed coefficient heights is used.

### Step 2: these tensor evaluations are injective

External invariant-theory input gives homogeneous dimensions at most
`1,2,7` at field degrees `4,6,8`, respectively
(arXiv:2509.14350v2, also consistent with the HLS low-order discussion).
The supplied graph/product evaluation matrices at the nine fixed tensors
have modular ranks `1,2,7`. Nonzero modular minors imply those same rational
matrices have full column rank. The matching upper bounds imply that these
basis functions span the entire corresponding invariant spaces.

Each target residual is itself a homogeneous Lorentz invariant in that
space. Because the evaluation map on that space is injective and the
residual vanishes at every chosen tensor, it is the zero invariant. This
establishes the declared degree-eight map over QQ **conditional on the
external invariant-theory bounds and the physics/convention inputs**.

## Verification levels

`verify-map8` rechecks the explicit arithmetic certificate, expected monomial
envelope, fixed-input bounds, residue equations, injective evaluation ranks,
and content hash. It does not trust a stored `passed` flag.

`prove-map8` performs/replays the tensor evaluations with content-addressed
checkpoints. `verify-map8 --reevaluate-tensors` recomputes the actual tensor
values from the recorded electric components without trusting stored
sample values. Both modes are available; the latter is intentionally not
part of every fast CI job.

As with any computer-assisted calculation, checking the evaluator and its
normalizations matters. The tensor implementation has Lorentz-boost,
self-duality, constrained-gradient/Euler, integer-vs-BLAS, and separate
multi-prime regression tests. The exact-BLAS path is allowed only when all
integer dot products fit strictly below `2^52`; a strict int64 path is also
available. Memory figures are conservative array estimates, not OS quotas.

## Limits

The published homogeneous dimension bounds are external inputs, not proved
by this repository. The HLS action/stress sign convention and the relation
between auxiliary self-dual fields and physical fields need expert review.
The certificate is not a statement about degree ten/twelve, a global
presentation of the invariant ring, supersymmetry, causality, or novelty.
