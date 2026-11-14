# Degree-10 dimension dictionary

Every number below dimensions a *different* object. Several are numerically
equal, and one apparent agreement turned out to be a coincidence when it was
finally computed.

## 1. The numbers

| n | dimensions what | grading | arithmetic | source |
|---:|---|---|---|---|
| **126** | components of a self-dual 5-form in D=10 | — | exact | analytic |
| **45** | dim SO(1,9) | — | exact | analytic |
| **81** | generic **functional** dimension = 126 − 45 | cumulative, all degrees | analytic (upper bound); exact modular (lower bound) | literature for the bound; this work for the matching computation |
| **14** | degree-10 **atlas** `A10`: polynomial invariants of F-degree 10 | graded, degree 10 | exact modular | trace repo |
| **12** | **published-candidate span** `B10` inside that 14 | graded, degree 10 | exact modular | trace repo |
| **12** | **graph-generator span** `G10` = span(I10_1..I10_12) | graded, degree 10 | exact modular | trace repo |
| **2** | **product subspace** `P10` = span(I4_1·I6_1, I4_1·I6_2) | graded, degree 10 | exact modular | trace repo |
| **11** | **reachable closure** `D10` (stress-flow generated) | graded, degree 10 | exact modular | trace repo |
| **3** | **Q10** = `A10 / D10` | graded, degree 10 | exact modular | trace repo |

## 2. Three splittings of the same 14 — and which one is real

    14 = 12 + 2      published span B10   +  complement       (trace)
    14 = 12 + 2      graph generators G10 +  products P10     (trace)
    14 = 11 + 3      reachable closure D10 + quotient Q10     (trace)

The third line is **unrelated** to the first two: `Q10` is a quotient by the
*stress-flow reachable* subspace, not by products. Nothing about `Q10` follows
from a primitive/product decomposition, and "independent" must never be used for
both.

The first two lines both read `12 + 2`, which invited the identification "the
published candidates span the primitive part". **That identification is false**,
and this is the correction that matters most in this file.

## 3. The resolved product/primitive question

Computed exactly at two primes
(`results/intrinsic_candidates/degree10_published_product_intersection.json`):

    dim B10 = 12,  dim P10 = 2
    dim(B10 ∩ P10) = 1        dim(B10 + P10) = 13   ← not 14

So the published span **contains one product direction** and misses one atlas
direction. Its primitive content is **11, not 12**. The `12 = 12` match with the
spinor-side primitive count is a dimensional coincidence, not a structural one.

The decomposition that *is* structural:

    P10 ∩ G10 = 0,   P10 + G10 = 14     ⟹     A10 = G10 ⊕ P10   (12 ⊕ 2)

with `G10` the span of the twelve graph generators. So the spinor-side
`12 primitive + 2 product` splitting corresponds to `G10 ⊕ P10`, **not** to
`B10 ⊕ complement`. `B10` is a third, distinct 12-dimensional subspace meeting
`G10` in only 10 dimensions.

Also verified rather than assumed: `P10 ⊆ D10`. Products contribute nothing to
`Q10`.

## 4. Complete pairwise incidence at degree 10

From `results/intrinsic_candidates/degree10_space_incidence.json`, identical at
both primes:

| pair | dim A | dim B | A ∩ B | A + B | containment |
|---|---:|---:|---:|---:|---|
| A10, B10 | 14 | 12 | 12 | 14 | B10 ⊂ A10 |
| A10, D10 | 14 | 11 | 11 | 14 | D10 ⊂ A10 |
| A10, G10 | 14 | 12 | 12 | 14 | G10 ⊂ A10 |
| A10, P10 | 14 | 2 | 2 | 14 | P10 ⊂ A10 |
| B10, G10 | 12 | 12 | 10 | 14 | neither |
| D10, B10 | 11 | 12 | 9 | 14 | neither |
| D10, G10 | 11 | 12 | 9 | 14 | neither |
| P10, B10 | 2 | 12 | **1** | 13 | neither |
| P10, D10 | 2 | 11 | 2 | 11 | **P10 ⊂ D10** |
| P10, G10 | 2 | 12 | **0** | 14 | complementary |

## 5. Spinor-side numbers, and how they now compare

| n | dimensions what | arithmetic |
|---:|---|---|
| 126 | gamma-traceless symmetric chiral square `Sym²(16) ⊖ 10` | exact modular (bridge) |
| 1, 2, 6, 14 | spinor **evaluation rank** at F-degree 4, 6, 8, 10 | exact modular (bridge) |
| 1, 3, 9, 21, 81 | **cumulative Jacobian rank** at F-degree 4, 6, 8, 10, 12 | exact modular (`results/rank81/certificate.json`) |

These are different quantities and must not be read off one another. Two
corrections to an earlier revision of this table, both of which were exactly the
error this document exists to prevent:

1. The cumulative Jacobian ranks were given as `1, 3, 8, …`. The certificate
   records `1, 3, 9, 21, 81`. The per-degree block ranks are `1, 2, 6, 12, 62`
   against candidate counts `1, 2, 6, 12, 62` — every block is full rank.
2. The text then compared `1 + 2 + 7 = 10` against that rank and concluded two
   functional relations exist by degree 8. Both halves are wrong. The `7` is the
   **trace-side** graded dimension; the Jacobian is taken over the **spinor**
   candidate family, which has 6 at degree 8. Mixing the two sides in a single
   subtraction is not a meaningful quantity.

Where the functional relations actually appear: nowhere before degree 12. The
full selection has 83 candidates and Jacobian rank 81, so exactly **two**
functional dependencies exist among the selected functions, and both are
degree-12 phenomena. The candidates are not an independent set and are never
described as one.

## 6. Words that are banned without a qualifier

| word | must be qualified as |
|---|---|
| rank | atlas rank / published rank / quotient rank / evaluation rank / Jacobian rank |
| independent | linearly independent / functionally (algebraically) independent |
| primitive | non-product generator — and **never** "published candidate" |
| dimension | of which space, at which degree, cumulative or graded |
| product | the lower-degree product subspace `P10`, spanned by the two explicit entries |

## 7. Change log

- **2026-08-01.** Section 5 corrected. The cumulative Jacobian ranks were wrong
  (`1, 3, 8, …` for `1, 3, 9, 21, 81`) and the accompanying text subtracted a
  spinor-side rank from a trace-side dimension to conclude that two functional
  relations exist by degree 8. Neither the number nor the comparison was sound.
  The two functional dependencies are real but they are degree-12, and they come
  from 83 candidates having rank 81. Recorded rather than quietly overwritten,
  because this document's whole purpose is to stop that class of mistake and it
  had made one.
- **2026-07-31.** Section 3 previously said the product/primitive
  relation was uncomputed and must not be asserted. It has since been computed:
  the tempting identification is **false**, and the true structural
  correspondence is `A10 = G10 ⊕ P10`. Sections 4 and 5 added. The
  "UNVERIFIED — do not assert" wording is removed because the question is now
  settled, not because the caution was wrong.
