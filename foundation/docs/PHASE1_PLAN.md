# Phase 1 implementation plan

**Objective**: identify *intrinsically* and prove necessary the seven remaining
generalized-flow directions — three at degree 10, four at degree 12.

**Status entering this plan**: the deficits are located and removal-minimal
(commit `226c4d8`). What remains is the part the Phase 1 gate actually
requires and that no amount of scanning supplies.

## 1. The precise problem

`I10_6`, `I10_7`, `I10_12`, `I12_59`, `I12_60`, `I12_61`, `I12_62` are **graph
labels in one basis**. A label is not a tensor. This is exactly the position
`I6_2` occupied before it was replaced by

    K6 = N1050_[abc,de]f N1050^[abc,]_[ghi] N1050^[def,gi]h

Two facts make the task tractable rather than open-ended:

- the closure is **coordinate-aligned** (C-MIN-04): the reachable set is
  spanned by 68 of the 72 basis vectors, so the quotient has clean coordinate
  representatives;
- the quotient is small — 3-dimensional at degree 10, 4-dimensional at
  degree 12.

But the *representatives* are basis-dependent. The intrinsic object is the
quotient class, and what must be exhibited is a tensor expression whose class
spans it.

## 2. Strategy — span-solving, not guessing

The successful K6 route was: propose a tensor, evaluate it exactly, solve for
its coordinates. Generalise that into a search that either succeeds or
returns a precise negative.

**Step A — build an intrinsic candidate library** at degrees 10 and 12 from
structures that are manifestly intrinsic:

| family | source |
|---|---|
| `Tr(M^k)` and products with total degree 2k | `src/sdinv/stress.py` |
| contractions of `N^(1050)` and its irreducible components | `src/sdinv/sextic.py` (`composite_n1050`) |
| products of lower-degree intrinsics (J6, K6, I4) | already available |
| stress-adapted characteristic invariants | `dimension_table.json` stress bases |
| mixed M–N contractions | to be enumerated |

**Step B — evaluate the library exactly** at the four deterministic samples
over ≥5 fit primes, in the same normalisation (`tau = 48T`) as the certificates.

**Step C — solve.** Let `V` be the span of the library in the degree-d value
space and `Q` the quotient representatives. For each missing direction, solve
for its coordinates in `V`. Three outcomes, all informative:

1. **exact solution** → the direction is intrinsic, expression recorded,
   verified on holdout primes and fresh samples;
2. **solution exists but coefficients exceed the CRT bound** → same situation
   as `Tr(M⁶)`; keep the direction as a primitive intrinsic element and prove
   span statements without coordinates (Phase 2 route B);
3. **no solution** → the library is too small. Report the **searched space and
   bounds** (standard 15) and state that a new tensor structure is required.
   This is a legitimate negative result, not a failure to be hidden.

Outcome 3 is the case where mathematical insight is genuinely required, and
the plan must not pretend otherwise.

## 3. Whether one generator produces several directions

Testable directly and cheaply with existing machinery: adjoin each candidate
singly and record which deficits collapse. Concretely, does any degree-10
direction, once adjoined, reduce the degree-12 deficit below 4? The degree-10
scan already showed it does **not** (deg12 stayed 68 for every degree-10
addition), which is itself a recorded result — the deficits at the two degrees
are dynamically independent.

The full **subset closure lattice** over the 7 directions is 2⁷ = 128 subsets
× 96 s ≈ 3.5 h on four primes, or ~50 min on one prime with a four-prime
confirmation of the extremes. Shardable per subset. Worth doing because the
lattice, not the cardinality, is the real minimality statement.

## 4. Minimality, properly

Present status is removal-minimality in one fixed basis. The gate and PO-08
require more:

| sense | method | cost |
|---|---|---|
| removal | done | — |
| basis change | 20 random invertible grading-preserving transforms, re-scan | ~35 min |
| alternative representatives | pick different coset representatives, confirm same quotient dimension | ~20 min |
| ≥5 fit primes | extend from 4 to 6 certificates | 2 × 530 s |
| ≥2 holdout primes | already discharged for the flow equations (PO-06); repeat for closure | ~20 min |
| fresh samples | new deterministic seeds, full recompute | 530 s/prime |

All estimable, all fit in RAM, all shardable.

## 5. Falsification tests owed by Phase 1

From the objective's list, these bear directly on Phase 1 and are **not yet
run**:

- **test 4** — random basis transformation (also PO-08);
- **test 5** — reorder candidates and shards, confirming the missing sets are
  not an artifact of atlas ordering. This matters more than it looks: the four
  degree-12 directions are the *last four labels*, which is precisely what an
  ordering artifact would look like;
- **test 7** — special non-generic field configurations (also C-SCOPE-02);
- **test 18** — every claimed minimal set fails when an element is removed
  (done for removal, owed for the lattice).

Test 5 should be run **before** any intrinsic search, because if the sets are
ordering artifacts the search targets are wrong.

## 6. Ordering of work

1. falsification test 5 (reordering) — cheap, and it validates the target set;
2. basis-change minimality (PO-08) — cheap, discharges an obligation;
3. extend to 6 fit primes + 2 holdouts for the closure;
4. build the intrinsic candidate library;
5. span-solve for all seven directions;
6. subset closure lattice;
7. write `minimality_proof_through_degree12.md` and certificates.

Steps 1–3 and 6 are mechanical and estimable. Steps 4–5 are where the phase
either completes or produces a precise negative.

## 7. Gate

Phase 1 is **not** complete until all seven directions have intrinsic
expressions *or* a stated-bounds negative result, and minimality is certified
under basis change as well as removal. Locating the deficits — already done —
is not the gate.
