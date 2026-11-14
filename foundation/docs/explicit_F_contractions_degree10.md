# Level-A representatives at degree 10

**Status: Level A COMPLETE for all three degree-10 quotient classes.**
Level B (M / N^(1050) / N^(4125)) and Level C (canonical compact) are **not
derived**.

Artifact: `results/intrinsic_candidates/explicit_F_contractions.json`
Gate: `tests/test_graph_to_tensor.py`

## 1. Why a graph is already an intrinsic representative

A verified contraction graph defines a Lorentz-scalar polynomial in F. What
was missing was not the mathematics but the *writing down*: an explicit
Einstein-index formula, independent of any numerical coordinate basis.

An earlier report in this program described the intrinsic rank as "0 of 3".
That was wrong. It conflated *no compact named formula* — true — with *no
intrinsic representative* — false. The translation below is deterministic and
was always available.

## 2. The translation rule

For a graph on n vertices with every vertex of valence 5:

    I  =  [ prod_v  F_{ s(v,1) ... s(v,5) } ]  x  [ prod_pairs  eta^{ s_a s_b } ]

- one self-dual five-form per vertex;
- one contracted index pair per unit of edge multiplicity;
- one inverse metric per pair;
- dummy names assigned canonically by edge order.

At degree 10 this gives **25 metric factors** and 50 index slots.

Implementation note: the evaluator uses **one einsum letter per contracted
pair**, applying the metric by raising one slot, rather than carrying eta as an
explicit operand. Writing eta explicitly needs two letters per pair — 60 at
degree 12 — which exceeds the 52-letter budget and adds 30 operands. The
equivalent one-letter form runs 40x faster and in 3.4x less memory
(0.1 s / 807 MB versus 4.0 s / 2.7 GB at degree 10).

## 3. The three classes

| intrinsic ID | graph-basis label | metric factors |
|---|---|---|
| `Q10_A` | I10_6 | 25 |
| `Q10_B` | I10_7 | 25 |
| `Q10_C` | I10_12 | 25 |

The intrinsic IDs are deliberately distinct from the graph labels, which are
basis-dependent coordinates. The correspondence is recorded, not identified.

Structure of `Q10_A`, abbreviated:

    F^(0)_[i0 i1 i2 i3 i4]  F^(1)_[i5 i6 i7 i8 i9]
    F^(2)_[i10 i11 i12 i13 i14]  F^(3)_[i15 i16 i17 i18 i19]  ...

with 25 metric contractions `eta^{i_k j_k}`. The full specification — every
slot, every pair, and LaTeX — is in the artifact.

## 4. Validation

Each class: 3 primes (32749, 32719, 32693) x 6 samples (4 fitting + 2 fresh)
= **18 dense-vs-graph checks**, all agreeing exactly, with **no sign
discrepancy**.

The dense evaluator is an independent code path: it builds einsum subscripts
from the dummy-index assignment and applies the metric by raising one slot per
pair, rather than using the repository's slot-planner. Agreement is therefore
a genuine check on the translation, not a tautology.

Homogeneity verified: `F -> cF` scales the scalar by `c^10` on every class.

## 5. What Level A is not

These are 10-fold contractions with 25 metric factors. They are explicit,
unambiguous and coordinate-independent — and that is all. They are **not**
compact, **not** canonical, and **not** conceptually explanatory. Any claim
that the degree-10 quotient has been *understood* requires Level B at minimum,
which is not derived.
