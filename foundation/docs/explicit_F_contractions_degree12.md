# Level-A representatives at degree 12

**Status: Level A COMPLETE for all four degree-12 quotient classes.**
Level B (M / N^(1050) / N^(4125)) and Level C (canonical compact) are **not
derived**.

Artifact: `results/intrinsic_candidates/explicit_F_contractions.json`
Gate: `tests/test_graph_to_tensor.py`

## 1. The translation rule

For a graph on 12 vertices, every vertex of valence 5:

    I  =  [ prod_v  F_{ s(v,1) ... s(v,5) } ]  x  [ prod_pairs  eta^{ s_a s_b } ]

giving **30 metric factors** and 60 index slots. Dummy names are assigned
canonically by edge order, so the output is reproducible.

The label dialect differs from degree 10: vertex indices reach two digits, so
edges are written `0-11^4` rather than concatenated. The parser handles both.

## 2. The four classes

| intrinsic ID | graph-basis label | metric factors |
|---|---|---|
| `Q12_A` | I12_59 | 30 |
| `Q12_B` | I12_60 | 30 |
| `Q12_C` | I12_61 | 30 |
| `Q12_D` | I12_62 | 30 |

The intrinsic IDs are deliberately distinct from the graph labels, which are
basis-dependent coordinates selected by an echelon form. The correspondence is
*recorded*, not *identified*.

**Note on `Q12_C` and `Q12_D`.** These correspond to I12_61 and I12_62, the two
candidates the atlas records under `discovery.functional_dependencies` —
polynomially independent members of the rank-72 basis, but not new generic
functional directions. That both are also dynamically unreachable is conjecture
**CJ-01**, with no mechanism proposed, and nothing here explains it.

## 3. Validation

Each class: 3 primes (32749, 32719, 32693) x 6 samples (4 fitting + 2 fresh)
= **18 dense-vs-graph checks**, all agreeing exactly, no sign discrepancy.
Across all seven classes at both degrees: **126 checks, 126 agreements**.

The dense evaluator is an independent code path — einsum subscripts built from
the dummy-index assignment, metric applied by raising one slot per pair —
rather than the repository's slot-planner.

Homogeneity verified: `F -> cF` scales by `c^12` on every class.

## 4. Implementation note that made degree 12 feasible

Carrying `eta` as an explicit einsum operand needs two letters per contracted
pair: 60 at degree 12, which exceeds the 52-letter budget and fails outright.
It also adds 30 operands. Using one letter per pair, with the metric applied by
raising one slot, is the same contraction and reduced degree 10 from
4.0 s / 2.7 GB to 0.1 s / 807 MB, with degree 12 at 0.4 s.

That 3.4x memory reduction is what brings the Stage-2 quadratic-block search
within reach of an 8 GiB machine at all; the earlier formulation exceeded the
1.5 GB working ceiling on a *single* evaluation.

## 5. What Level A is not

These are 12-fold contractions with 30 metric factors. Explicit, unambiguous,
coordinate-independent — and nothing more. They are **not** compact, **not**
canonical, and **not** explanatory. Claiming the degree-12 quotient is
*understood* requires Level B at minimum, which is not derived.
