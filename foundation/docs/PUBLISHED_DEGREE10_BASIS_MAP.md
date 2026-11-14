# The twelve published degree-10 candidates in the atlas and in Q10

## 1. Result

**The twelve equation-(4.24) candidates span Q10 completely: rank 3 of 3.**

Artifact `results/intrinsic_candidates/published_degree10_map.json`, schema 2.

| prime | role | Q10 rank | atlas rank | runtime |
|---|---|---:|---:|---:|
| 32749 | fit | **3 / 3** | **12 / 14** | 483 s |
| 32717 | holdout | **3 / 3** | **12 / 14** | 429 s |

Quotient vectors, per candidate:

| candidate | 32749 | 32717 | reaches Q10 |
|---|---|---|---|
| P10_01 … P10_08 | `[0,0,0]` | `[0,0,0]` | no |
| **P10_09** | `[28088, 575, 0]` | `[30992, 15284, 0]` | **yes** |
| **P10_10** | `[26346, 9877, 9547]` | `[15788, 13584, 14050]` | **yes** |
| **P10_11** | `[15468, 21541, 0]` | `[1539, 18934, 0]` | **yes** |
| **P10_12** | `[27310, 1310, 0]` | `[9234, 28791, 0]` | **yes** |

Every one of the thirty-two entries has `status = "solved"`: each candidate
lies in the atlas span and its coordinates were obtained exactly. A candidate
that merely failed to solve would also contribute rank 0, for an entirely
uninformative reason.

## 2. A correction worth keeping

An earlier state of this document reported that every published candidate had
zero image in Q10, and drew a cross-degree pattern from it together with the
degree-12 result.

**That was wrong**, and the way it was wrong is instructive. Only the five
simplest candidates were implemented at the time — P10_01, 02, 03, 06, 07 —
and all five do project to zero. So do P10_04, P10_05 and P10_08. The four
that reach the quotient are precisely the four that were hardest to transcribe
and were implemented last: P10_09, carrying a red bracket, and P10_10, P10_11,
P10_12, carrying five `N^(1050)` blocks each.

The subset that had been implemented was selected for ease of implementation,
which is correlated with structural simplicity, which is exactly what
determines whether a candidate can reach the quotient. **A null result over a
subset chosen for convenience is not evidence about the whole set.**
`test_only_the_four_hardest_candidates_reach_the_quotient` now pins which
candidates carry the quotient, so a transcription error that moved a nonzero
image between candidates would be caught even though the rank would not change.

The degree-12 statement (`C-P12-01`) is unaffected: it concerns three specific
structures from equation (4.25) and stands as recorded. What is withdrawn is
the cross-degree *pattern* claim, which rested on the incomplete degree-10 set.

## 3. Corroboration of the paper's own count

The published atlas rank is **12 of 14** at both primes. The paper states that
"at the 10th order there are 12 linearly independent invariants". The twelve
candidates span a 12-dimensional subspace of the 14-dimensional degree-10
atlas — the remaining two atlas directions being products of lower-degree
invariants rather than order-10 primitives.

This is an independent numerical corroboration of a stated claim, obtained
without using that claim as an input. It is **not** a proof of linear
independence as the paper means it; it is agreement between two counts.

## 4. Q10 rank contributions, cumulatively

Adding candidates in index order, the running Q10 rank is:

    P10_01  0    P10_05  0    P10_09  1
    P10_02  0    P10_06  0    P10_10  2
    P10_03  0    P10_07  0    P10_11  3
    P10_04  0    P10_08  0    P10_12  3

`P10_12` adds nothing beyond `P10_09`, `P10_10`, `P10_11` — but it is not
redundant in the basis sense; see the Level-B selection, where it is preferred
over `P10_11`.

## 5. Why the projector's zeros are trustworthy

The positive control (`degree10_positive_control.json`) pushes the Level-A
representatives Q10_A/B/C — established independently as explicit F-index
contraction graphs — through the identical `closure_span` → `rref` → `project`
chain and recovers rank **3 of 3** at all six primes.

Without that control, a projector stuck at zero would print the same zeros for
P10_01…P10_08 and nothing would distinguish the two situations.

## 6. Scope

**Permitted**: "the twelve degree-10 structures displayed in equation (4.24)
span the quotient Q10, of dimension 3, under the stated formal flow
definition, at two primes with fit/holdout agreement; eight of the twelve lie
inside the reachable closure D10."

**Forbidden**:

1. **Not** an all-orders or degree-12 statement. Q12 is untouched by this.
2. **Not** independent of the AMB-01 / AMB-02 source readings, except where
   robustness is explicitly measured — see `intrinsic_degree10_levelB.md` §3.
3. **Not** a claim that these are the *simplest* structures reaching Q10. No
   larger class has been enumerated.
