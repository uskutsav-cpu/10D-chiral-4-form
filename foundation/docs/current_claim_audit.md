# Claim audit — Phase 0

Every scientific statement in the existing reports, audited against raw
artifacts and code. Statuses feed `docs/CLAIM_LEDGER.md`.

## 1. Claims requiring CORRECTION

None found that are currently *asserted* in the repository. The four
correction targets named in the objective were checked individually:

| target | finding |
|---|---|
| static span vs. closure | **already correct.** `stress_flow_classification.md` §1 states the closure (1,1,3,11,67) strictly exceeds the static span (1,1,2,2,4) and that a no-go built on the static span "would be wrong." `stress_flow_definition_audit.md` records the universal-`I6_2` obstruction as "contradicted as an inference from the static map alone; not established." |
| 3-prime failure vs. 5-prime success | **was stale, now fixed.** `assumptions_limitations_and_open_questions.md` §2.2 previously asserted the flow equations had no rational lift. Corrected in place to record that three primes fail and five succeed with holdout validation. |
| arbitrary `I6_2` vs. intrinsic K6 | **already correct.** `I6_obstruction_proof.md` explicitly warns `I6_1` is not `Tr(M³)` (it is `3J6/32`) and `I6_2` is not `K6`. Intrinsic quotient `q6` is used throughout. |
| universal-obstruction vs. transported-not-created | **already correct.** All documents use "transported, never created." Ledger entry F-1 forbids "K6 must vanish in every pure stress flow." |

## 2. Claims requiring WEAKENING

| claim | current wording risk | required wording |
|---|---|---|
| C-ATLAS-03 | rank 81 saturation could be read as "we have all 81 invariants" | "the cumulative Jacobian rank saturates 81"; **never** "the atlas is complete." Rank saturation is not a generating set. |
| C-ATLAS-04 | I12_61/I12_62 "do not raise the rank" sits close to "syzygy" | must never be called a syzygy until PO-03 produces a verified polynomial identity |
| C-FLOW-02 | "valid for all λ" | must carry *formal power series* and *at field degree six*; convergence in λ is unestablished |
| C-MIN-01 | "minimal" | minimal **under removal, in a fixed basis**; basis-change minimality is PO-08 |
| C-SEXTIC-* | K6 presented as the intrinsic sextic | true as our definition; must not be identified with the published Σ₂ (PO-04) |
| all `MOD-CERT` | prime agreement reads as proof | rank mod p is a characteristic-zero *lower bound*; exceptional primes are PO-09 |

## 3. Claims that may be STRENGTHENED

| claim | basis |
|---|---|
| C-GEN-01 | upgraded to `EXACT-CA-THM`: generator enumeration verified complete, 18 = 18, additivity holds for every generator |
| C-ATLAS-05 | PO-01 discharged — WL removed, canonicalisation exact or raising |
| C-FLOW-01 | PO-06 discharged — two independent holdouts, 192/192 identical reconstructions across two fit sets |
| C-FLOW-02 | the degree-6 exhaustiveness argument is now *proved* complete rather than observed, since only generators with leading degree ≤ 6 can contribute and all three produce rows |

## 4. Statements that are NOT yet claims

Listed to prevent accidental promotion:

- any all-orders statement (PO-10 open);
- any physical or causal statement (Phase 7 not started);
- any Type IIB statement (Phase 8 not started; PO-07 gates it);
- "independently reproduced" (Phase 10 not started — the dense N^(1050)
  contraction is a second implementation *inside* this repository, which is
  useful but is not clean-room);
- "new theory" (Phase 6 not started).

## 5. Numerical baseline re-verified against artifacts

| reported | artifact | agrees |
|---|---|---|
| dims 1, 2, 7, 14, 72 | `interacting_flow_equations.json:basis_dimensions` | yes |
| static 1, 1, 2, 2, 4 | `dimension_table.json:dimension_rows` | yes |
| new forcing 1, 1, 3, 5, 21 | `interacting_flow_equations.json` | yes |
| complement 0, 1, 4, 9, 51 | same; and forcing + complement = full at every degree | yes |
| closure 1, 1, 3, 11, 67 | `closure_and_minimality.json:closures.free` | yes |
| minimal set {K6, I8_3..I8_6} | `closure_and_minimality.json:minimal_completion_through_degree8` | yes |
| deficits 3 at deg 10, 4 at deg 12 | `closure_and_minimality.json:still_open` | yes |
| 96 tests | fresh-clone run | yes |

One discrepancy against the objective's stated baseline, noted for the record:
the objective lists primitive candidate counts `1 + 2 + 6 + 12 + 62 = 83`,
whereas the *homogeneous space* dimensions are `1, 2, 7, 14, 72`. These are
different quantities (primitives vs. full spaces including products) and are
not in conflict; the ledger keeps them as separate claims C-ATLAS-01 and
C-ATLAS-02.
