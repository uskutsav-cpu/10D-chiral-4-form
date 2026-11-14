# The intrinsic sextic basis (J6, K6)

Artifact: `results/stress_flow/change_of_basis/sextic_intrinsic.json`
Gates: `tests/test_sextic.py` (4 tests)
Longer derivation: `docs/I6_obstruction_proof.md`

The repository's degree-six labels `I6_1`, `I6_2` are graph names, not
tensorial objects. Any statement about "the obstruction being `I6_2`" is
therefore basis-dependent and not publishable as stated. This document
records the intrinsic replacement.

## 1. Definitions

    J6 = Tr(M^3) = M_mu^nu M_nu^rho M_rho^mu

    N1050 = Lambda^mn_[abc Lambda_de]fmn

    K6 = N1050_[abc,de]f  N1050^[abc,]_[ghi]  N1050^[def,gi]h

with normalised antisymmetrisation throughout.

## 2. Exact change of basis

Read column-wise against the ordered repository basis `(I6_1, I6_2)`:

    ( J6  K6 ) = ( I6_1  I6_2 ) * [ 32/3   -1/1125 ]
                                  [ 0       3/125  ]

so

    J6 = (32/3) I6_1
    K6 = -(1/1125) I6_1 + (3/125) I6_2

Inverse:

    [ 3/32   1/288  ]
    [ 0      125/3  ]

Determinant `32/125`, nonzero, which is simultaneously an independence proof
and an exact inverse. Checked: `(32/3)(3/125) - (-1/1125)(0) = 32/125`, and
`(3/32)(125/3) = 125/32`.

Note the consequences for naming: `I6_1` is **not** the paper-normalised
`Tr(M^3)`; it equals `3 J6 / 32`. And `I6_2` is **not** `K6`; it differs by
normalisation *and* a `J6` admixture.

## 3. Quotient coordinate

Normalising `q6([K6]) = 1`,

    q6(a J6 + b K6) = b      i.e.      q6(c1 I6_1 + c2 I6_2) = (125/3) c2

Under an invertible change of sextic basis the coordinate column and the row
covector transform contragrediently, so `q6(V6)`, its vanishing, and the
quotient class are unchanged. Under a field rescaling `Lambda -> s*Lambda`
both `J6` and `K6` scale by `s^6`, so membership in `<J6>` and non-vanishing
of the quotient class are preserved.

This is what makes the dynamical statement in
`docs/stress_flow_classification.md` basis-independent rather than an
artifact of the graph labelling.

## 4. Validation

Exact finite-field evaluation over primes **32749, 32719, 32693**, four
deterministic samples each, comparing three independent computations:

- the eight-term contraction-graph expansion of `K6`;
- an independently implemented dense `N^(1050)` cubic contraction;
- the verified registry basis.

On every sample `direct_1050_value == K_1050`, exactly. The graph expansion
carries explicit orientation signs; isomorphic multigraphs are **not**
identified without checking orientation, because contractions of odd forms
can pick up a sign under graph isomorphism.

## 5. Caveat carried from the source

Recorded in the artifact and repeated here because it limits what may be
claimed: arXiv:2509.14350v2 proves that `(Sigma_1, Sigma_2)` is a sextic
basis in the spinor formalism, but **does not publish the change of basis**
to `(Tr(M^3), K6)`. The identification of `K6` with the second spinor
invariant is ours, not the paper's, and is not independently confirmed.
