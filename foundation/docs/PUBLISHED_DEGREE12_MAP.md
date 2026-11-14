# The published degree-12 structures in Q12

## The certified claim, in exactly the permitted wording

> **The three degree-12 structures displayed in equation (4.25) lie inside the
> computed reachable closure D12 and have zero image in Q12 under the stated
> formal flow definition.**

Nothing stronger is claimed. See §4 for what this explicitly does *not* say.

## 1. Structures

From equation (4.25), *Some remarks on invariants*, J. Phys. A 59 (2026)
065203 (journal PDF p18; arXiv:2509.14350v2 p26), transcribed from the
rendered source:

    P12_01 = tr M^6
    P12_02 = (M^3)^{mu nu} (N^(4125) M M)_{mu nu}
    P12_03 = (N^(4125) M M)^{mu nu} (N^(4125) M M)_{mu nu}

Implemented in `src/sdinv/published_degree12_invariants.py` on the
repository's already-tested `n4125_mm`, `symmetric_inner` and
`matrix_trace_power`. No convention was re-chosen.

## 2. Result

| prime | dim Q12 | P12_01 | P12_02 | P12_03 | rank | seconds |
|---|---:|---|---|---|---:|---:|
| 32749 | 4 | `[0,0,0,0]` | `[0,0,0,0]` | `[0,0,0,0]` | **0** | 2925 |
| 32717 | 4 | `[0,0,0,0]` | `[0,0,0,0]` | `[0,0,0,0]` | **0** | 3249 |

`consistent_across_primes: true`. Artifact
`results/intrinsic_candidates/published_degree12_map.json`,
sha256 `9a8f4e32627ee5bef06324169dc823a4…`.

**The result is not vacuous.** Every structure has `status = "solved"` — each
lies in the atlas span and its coordinates were obtained exactly — with a
quotient vector of the correct length 4 that is identically zero. A structure
that merely failed to solve would also have contributed rank 0, for an
entirely uninformative reason; the `status` field separates the two cases and
the tests assert it.

## 3. Method

Per prime: 72 atlas columns evaluated at 80 generic samples, exact modular
linear solve for atlas coordinates, then projection through the echelon form
of D12. No floating point anywhere. The degree-12 registry stores
`graph_record` rather than `graph`, and products resolve recursively through
their factors with a cycle guard.

`tr(M^6)` retains **no** certified rational graph-basis coordinate vector —
29 of 72 columns exceeded the CRT uniqueness bound at 15 primes. Only the
modular span and quotient statements above are certified, which is sufficient
for a rank claim and is all that is asserted.

## 4. What this does NOT establish

Explicitly forbidden inferences, listed because each is a plausible misreading:

1. **Not** that the complete M/N degree-12 space has zero image in Q12. Three
   specific structures were tested.
2. **Not** that Q12 has no compact tensor representation. Q12 has dimension 4;
   these three structures simply do not reach it.
3. **Not** that the three invariants are redundant. They are nonzero elements
   of the degree-12 space; they lie in D12.
4. **Not** that no generalized flow can generate Q12. That is a different
   question about generators, not about these three scalars.

## 5. Consequence

Combined with the degree-10 results — `P10_01 = tr M^5` and `P10_02` both
projecting to zero, and the complete M-only family having quotient rank 0 —
every published compact structure tested so far lies inside the reachable
closure at its degree.

So all four compact Q12 directions, and all three compact Q10 directions,
remain **unknown**. That is the honest state of the Level-B problem, and it
is what the reverse graph-to-block engine exists to attack.
