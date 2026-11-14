# Exactly what is claimed for the Q10 Level-B basis

## 1. The claim, in the only permitted wording

> **Preferred Level-B basis among the twelve published degree-10 candidates
> under the documented deterministic simplicity score, verified at one fit and
> one holdout prime.**

    { P10_10, P10_11, P10_12 }

**Not** universally canonical. **Not** the simplest structure reaching Q10.
**Not** ambiguity-free — see §3, which is the most important section here.

## 2. Why "canonical" is unavailable

Canonicality would require a proof over a search space that has not been
enumerated. Only twelve structures were tested — the ones equation (4.24)
happens to display. The space of compact `M / N^(1050) / N^(4125)` contractions
of degree 10 is very much larger, and nothing here bounds it. A simpler
representative could exist outside the published list and this work would not
have seen it.

That is what the reverse graph-to-block benchmark exists to probe.

## 3. The basis is NOT ambiguity-robust, and cannot be

This correction matters more than the selection itself, so it is stated before
the score rather than after it.

Bracket **colour** does not survive PDF text extraction, and colour is what
fixes operation order in equation (4.24). Both readings of every ambiguous
candidate were implemented and projected. Measured:

| candidate | ambiguity | quotient image under the two readings |
|---|---|---|
| P10_10 | AMB-02 | **identical** |
| P10_12 | AMB-02 | **identical** |
| P10_11 | AMB-02 | **differs** |
| P10_09 | AMB-01 | **differs** |

Now count. `P10_10` is forced — it is the only candidate reaching the third
quotient coordinate, so every independent triple contains it. The remaining two
slots must be filled from `{P10_09, P10_11, P10_12}`, and only **one** of those
three is robust.

> **Therefore no ambiguity-robust three-element basis exists within the twelve
> published candidates.** Every independent triple contains at least one member
> whose quotient image moves when the source reading is resolved.

The selected basis contains exactly one such member, `P10_11`, which is the
minimum achievable. Calling the result "the robust basis" would be wrong; it is
the *ambiguity-minimal* basis. Only `P10_10` and `P10_12` are unconditional.

### How much this actually costs — measured, and less than stated above

An earlier version of this section left the impression that the basis itself is
in doubt. It is not, and the distinction is worth being exact about.

Measured at both primes (`Q10_ambiguity_relationship.json`):

| quantity | reading A | reading B | union |
|---|---:|---:|---:|
| rank of `{P10_10, P10_11, P10_12}` | 3 | 3 | **3** |

Under *either* reading the triple has rank 3 and therefore spans **all** of a
3-dimensional Q10. The union of both readings still has rank 3, so the two
readings span the same subspace.

> **AMB-02 changes the COORDINATES of `P10_11`, not the SPAN of the basis.**

So a *span-level* statement — "these three elements are a basis for Q10" — is
unaffected by the ambiguity. What is reading-dependent is any statement about
the specific coordinate vector of `P10_11`, and hence the change-of-basis
matrices. `ambiguity-minimal` remains the correct label for coordinate-level
claims; it should not be read as doubt about the basis property itself.

**Consequence.** Resolving AMB-01 and AMB-02 from a colour render of journal
page 17 is a binding scientific prerequisite for an unconditional basis, not a
tidiness item.

## 4. The three independent triples, and why this one wins

| triple | score | non-robust members |
|---|---:|---|
| **{P10_10, P10_11, P10_12}** | **222.2** | P10_11 |
| {P10_09, P10_10, P10_12} | 233.1 | P10_09 |
| {P10_09, P10_10, P10_11} | 283.7 | P10_09, P10_11 |

`{P10_09, P10_11, P10_12}` is dependent and is not a basis.

The top two triples each carry exactly one non-robust member, so the ambiguity
penalty does not separate them; the remaining score does. `P10_09` loses to
`P10_11` because it carries a RED symmetrisation expanding to 6 permutation
terms plus an extra `M` block, against `P10_11`'s single unbracketed
five-`N^(1050)` contraction.

The third triple carries two non-robust members and is eliminated by the
penalty alone.

## 5. Score definition

Lower is simpler. Weights are stated, not tuned.

| term | weight | rationale |
|---|---:|---|
| source-ambiguity penalty | 50 | dominates; a basis element that moves when the source is resolved is not a stable answer |
| block count (M + N1050 + N4125) | 10 each | the dominant structural cost |
| BLACK bracket operations | 5 each | explicit symmetrisation work |
| RED bracket operations | 8 each | staged, and strictly harder to transcribe correctly |
| expanded permutation terms | 1 each | a RED sym over k slots expands to k! |
| measured runtime, seconds | 1 each | tie-breaker only |

**Ambiguity robustness outranks evaluation cost by construction**: the penalty
is 50, larger than any achievable difference in the remaining terms across
these candidates. That ordering is a deliberate scientific judgement — a
cheaper formula whose meaning is uncertain is worth less than a costlier one
whose meaning is fixed.

### Terms declared but not separating

`coefficient height` and `manifest symmetry` appear in the specification and
are **not** implemented as score terms. They would not change the ranking here:
all three triples are drawn from the same four candidates, every basis
coefficient is a generic residue of comparable height at both primes, and no
candidate carries a manifest symmetry the others lack. They are recorded as
absent rather than silently folded in at weight zero.

`peak RAM` is likewise not a separating term — the four candidates evaluate
within the same allocation profile, dominated by the shared `composite_n1050`
tensor rather than by anything candidate-specific.

## 6. Stability of the selection

The independent-triple set is identical at both primes, so the selection is not
an artifact of one modulus. It **is** conditional on the AMB-02 reading
recorded for `P10_11`: if that reading is resolved against the implemented one,
`P10_11`'s quotient vector changes and the triple must be re-scored. The
selection machinery is deterministic and re-runnable, so this is a rerun rather
than a re-derivation.
