# A compact Level-B basis for Q10

## 1. The basis

    Q10 Level-B basis  =  { P10_10, P10_11, P10_12 }

Preferred under the simplicity criterion in
`scripts/build_Q10_levelB_basis.py`, and minimal among the twelve published
equation-(4.24) candidates. Artifacts:
`intrinsic_Q10_levelB_basis.json`, `Q10_basis_search.json`.

This closes Level B for degree 10. Level A gave an explicit ten-F-index
contraction graph per class; Level B expresses the same three quotient
directions in the `M / N^(1050) / N^(4125)` block language.

## 2. The search, in full

Four of the twelve candidates have nonzero image in Q10: P10_09, P10_10,
P10_11, P10_12. Of the four three-element subsets, **three** are independent,
and both primes agree on which:

    { P10_09, P10_10, P10_11 }
    { P10_09, P10_10, P10_12 }
    { P10_10, P10_11, P10_12 }     <- preferred

`{P10_09, P10_11, P10_12}` is dependent.

**Removal minimality.** `P10_10` appears in every independent triple, so it
cannot be dropped from any basis. The reason is visible in the vectors: at
32749 the images of P10_09, P10_11 and P10_12 all have third component 0, so
those three span at most a two-dimensional subspace. `P10_10` is the only
candidate reaching the third quotient direction. The same holds at 32717.

Any two of `{P10_09, P10_11, P10_12}` are pairwise independent in the first two
coordinates — pairwise determinants 19461, 1674, 9775 at 32749, all nonzero —
so the other two members are genuinely a choice, settled by the criterion
rather than by necessity.

## 3. Source robustness — and a selection that changed

This is where the criterion does real work, and the answer moved once the
evidence was complete.

Every candidate carrying an unresolved bracket ambiguity had **both** readings
implemented and projected separately. What matters is not whether the raw
values differ but whether the **quotient image** does:

| candidate | ambiguity | Q10 image under the two readings |
|---|---|---|
| P10_10 | AMB-02 | **identical** at both primes |
| P10_12 | AMB-02 | **identical** at both primes |
| P10_11 | AMB-02 | **differs** at both primes |
| P10_09 | AMB-01 | **differs** at both primes |

An earlier state of this document preferred `{P10_09, P10_10, P10_12}` on the
grounds that P10_09's robustness was merely *unmeasured* while P10_11's failure
was measured. That gap was then closed: the AMB-01 alternative reading of
P10_09 was implemented, and P10_09's quotient image **does** move —
`[30992, 15284, 0]` against `[29055, 18687, 0]` at 32717.

So P10_09 and P10_11 are both ambiguity-sensitive and take the same penalty.
With that tie, the rest of the score decides, and P10_11 wins because P10_09
carries a RED symmetrisation expanding to 6 permutation terms plus an extra `M`
block:

    222.2   P10_10, P10_11*, P10_12     <- preferred
    233.1   P10_09*, P10_10, P10_12
    283.7   P10_09*, P10_10, P10_11*

(`*` marks a member whose quotient image moves under its source ambiguity.)

**Consequence, stated plainly**: no independent triple is free of source
ambiguity, because P10_10 is forced and the two remaining slots must draw at
least one member from the sensitive set `{P10_09, P10_11}`. The basis is
therefore certified *conditionally on the AMB-01/AMB-02 readings named in the
registry*, and resolving those readings is a genuine open item, not a
formality. Only `P10_10` and `P10_12` are unconditional.

## 4. Scoring criterion

Lower is simpler. Weights are stated, not tuned:

    10 x (number of M + N1050 + N4125 blocks)
     5 x (explicit BLACK bracket operations)
     8 x (RED bracket operations)
     1 x (expanded permutation terms; a RED sym over k slots gives k!)
     1 x (measured single-evaluation seconds)
    50   if the quotient image is NOT invariant under an unresolved
         source ambiguity

Per-member figures, measured at prime 32749, seed 5:

| candidate | M | N1050 | N4125 | black | red | perm terms | seconds | robust |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| P10_09 | 1 | 4 | 0 | 0 | 1 | 6 | 5.1 | **no** |
| P10_10 | 0 | 5 | 0 | 0 | 0 | 1 | 5.4 | **yes** |
| P10_11 | 0 | 5 | 0 | 0 | 0 | 1 | 7.2 | **no** |
| P10_12 | 0 | 5 | 0 | 0 | 0 | 1 | 6.6 | **yes** |

Peak RSS for the full two-prime projection was 531 MB; the resumed run,
reusing checkpoints, was 267 MB.

## 5. What is and is not claimed

**Permitted**: "`{P10_10, P10_11, P10_12}` is a three-element Level-B basis for
Q10, preferred under the documented simplicity criterion and minimal among the
twelve published equation-(4.24) candidates, verified at one fit and one
holdout prime, conditional on the AMB-02 reading recorded for P10_11."

**Forbidden**:

1. **Not** canonical. No class larger than the twelve published candidates has
   been enumerated, so "the simplest structure reaching Q10" is unsupported.
2. **Not** unconditional. Every independent triple contains an
   ambiguity-sensitive member; see §3.
3. **Not** a statement about Q12, which remains at rank 0 from the three
   published degree-12 structures.
4. **Not** rationally reconstructed. Certified modularly at two primes; no CRT
   reconstruction to rationals has been performed.
5. **Not** validated on five fitting primes. Two primes, one fit and one
   holdout, agreeing. The remaining four are incremental thanks to
   checkpointing but have not been run.

## 6. Open

- Resolve AMB-01 and AMB-02 from a colour render of journal page 17. This is
  now the binding constraint on an unconditional basis, not a tidiness item.
- Extend from two primes to six; the projection is checkpointed.
- Rational reconstruction of the change-of-basis matrices.
