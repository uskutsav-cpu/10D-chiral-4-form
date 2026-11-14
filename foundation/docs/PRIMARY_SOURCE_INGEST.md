# Primary source ingest

**Correction of an earlier claim.** A previous session reported equation (4.24)
as an external blocker on the grounds that the paper was unavailable. That was
wrong — both versions are publicly downloadable. The formulas below are
transcribed from the rendered source, not reconstructed from memory.

PDFs are **not** committed. Only hashes, citations, page references and
independently written implementation notes are.

## 1. Sources

| role | version | pages | sha256 (first 32) |
|---|---|---:|---|
| **source of record** | J. Phys. A 59 (2026) 065203, doi:10.1088/1751-8121/ae3bb8 | 35 | `0437ffc2a5fce3cfcaf45601a972a966` |
| cross-check | arXiv:2509.14350v2 | 52 | `ec51d9e22c4d75651e2024d09b17562b` |

Cederwall, Hutomo, Kuzenko, Lechner, Sorokin, *Some remarks on invariants*.
Downloaded 2026-07-30. Full manifest:
`results/intrinsic_candidates/primary_source_manifest.json`.

Note the `file(1)` header on the arXiv PDF misreports the page count as 10;
PyMuPDF reads 52 correctly. Equation locations were found by searching the
extracted text, not by trusting the reported pagination.

| equation | arXiv v2 PDF page | journal PDF page |
|---|---:|---:|
| (4.24) — twelve degree-10 candidates | 25 | 17 |
| (4.25) — three degree-12 structures | 26 | 18 |

## 2. The bracket convention, quoted exactly

> "where the symmetrization and/or anti-symmetrization of the indices within
> the red brackets is made upon the anti-symmetrization within the black
> brackets"

This ordering is essential: red-bracket operations act **after** the
black-bracket antisymmetrizations. Any implementation of (4.24) that applies
them in the other order is wrong, and the colour distinction is exactly what
plain text extraction loses.

## 3. Equation (4.25) — transcribed and IMPLEMENTED

    I^(1)_12 = tr M^6

    I^(2)_12 = (M^3)^{mu nu} N^(4125)_{mu a1 a2, nu b1 b2} M^{a1 b1} M^{a2 b2}
             = (M^3)^{mu nu} (N^(4125) M M)_{mu nu}

    I^(3)_12 = (N^(4125) M M)^{mu nu} (N^(4125) M M)_{mu nu}

The paper introduces these as invariants that "have appeared in certain models
of the non-linear self-dual 5-form theory".

Implemented in `src/sdinv/published_degree12_invariants.py` as `P12_01`,
`P12_02`, `P12_03`, built on the repository's already-tested `n4125_mm`
(which *is* `(N^(4125) M M)_{mu nu}`), `symmetric_inner` and
`matrix_trace_power`. Nothing was re-derived; no convention was re-chosen.

Values at seed 20260729, prime 32749: `19148`, `20011`, `10961`.
Homogeneity `F -> cF ~ c^12` verified for all three.

## 4. Equation (4.24) — transcribed, NOT yet implemented

Twelve candidates `I^(1..12)_10`, built by "choosing and contracting different
irreps in the decompositions of symmetric products of M^(54), N^(1050) and
N^(4125)". The first is simply `tr M^5`; the remainder involve nested
black/red bracket structures on `N^(1050)` and `N^(4125)`.

Full extracted text is in the manifest. Implementation requires a second
independent parsing pass against **rendered page images**, because the
red/black bracket distinction does not survive text extraction — the extracted
stream shows bracket characters without colour, so the operation ordering is
not recoverable from text alone.

**Status: transcription obtained, implementation deferred pending image-based
verification of bracket ordering.** Implementing from the text stream alone
would risk exactly the silent convention error this program forbids.

## 5. A discrepancy in the primary source — recorded, not resolved

Both versions state:

> "At the 10th order there are 12 linearly independent invariants and at order
> 12 there are **64**. ... together with the lower-order invariants their number
> is **83**"

These are inconsistent:

    1 + 2 + 6 + 12 + 64 = 85   (not 83)
    1 + 2 + 6 + 12 + 62 = 83   (consistent)

The repository's atlas has **62** primitives at degree 12, which reproduces the
paper's own stated total of 83. The "64" appears to be a typo present in both
the arXiv and journal versions.

**No convention or count was altered to fit.** This is logged for the authors.

## 6. Independent corroboration of an existing repository result

The paper states that 83 candidates against a ring dimension of 81 implies
"at least 2 non-linear relations". The repository independently found
cumulative functional rank 81 with exactly two functional dependencies,
`I12_61` and `I12_62` (C-ATLAS-04).

**These are related but not identical statements.** The paper conjectures
*non-linear relations*; the repository established *Jacobian dependence*. A
Jacobian dependence is not a polynomial identity — that gap is PO-03 and
remains open. The agreement in the number 2 is striking and worth recording,
but it is not a proof that the two statements are the same.
