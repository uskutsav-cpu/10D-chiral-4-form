# Maximal program — status

Resume point for any future session. Read this first, then
`docs/CLAIM_LEDGER.md` and `docs/PROOF_OBLIGATIONS.md`.

- **Branch**: `research/maximal-chiral-four-form-program`
- **Baseline**: `3ed32805b38ce34216b34888f6539e3538e90fb9`
- **Current phase**: 0 complete → 1 starting
- **Nothing pushed** beyond the already-public classification branch.

## Phase status

| phase | title | status |
|---|---|---|
| 0 | baseline reproduction + claim audit | **COMPLETE** |
| 1 | finish generalized flow through degree 12 | **IN PROGRESS** |
| 2 | resolve Tr(M⁶) / intrinsic stress basis | not started |
| 3 | trace–spinor equivalence, real syzygies | not started (PO-04 external) |
| 4 | geometric structure | not started |
| 5 | all-orders theorem | not started |
| 6 | new nonlinear theory | not started |
| 7 | causality / stability | not started |
| 8 | Type IIB | not started |
| 9 | algorithm generalization | not started |
| 10 | clean-room reproduction | not started |
| 11 | adversarial referee review | not started |
| 12 | mentor package | not started |
| 13 | manuscripts | **BLOCKED** by design until 1–12 |
| 14 | reproducible release | not started |

## Pillar status

| pillar | state |
|---|---|
| 1 general structure | partial — degrees 10/12 closure deficits open (C-MIN-02) |
| 2 proof | **not started** — no induction exists (PO-10) |
| 3 new theory | not started |
| 4 string theory | not started |
| 5 independent reproducibility | not started (no clean-room) |

## Established at Phase 0

- 96/96 tests from a fresh clone at the verified commit; fingerprint
  `26b61c44…` reproduces; 6D gate PASS; closure artifact regenerates
  byte-identically.
- **PO-01 discharged**: WL heuristic removed; canonicalisation is exact
  (pynauty) or raises.
- **PO-06 discharged**: two independent holdout primes, 192/192 identical
  reconstructions across two fit sets.
- **PO-02 partially discharged**: generator enumeration verified complete
  (18 = 18, additive leading degrees).

## Open, in priority order

1. **C-MIN-02** — degree-10 deficit 3, degree-12 deficit 4. Phase 1 target.
2. **PO-05** — Tr(M⁶) has no certified rational coordinates. Route (B),
   basis-independent formulation, is the likely resolution.
3. **PO-03** — I12_61/I12_62 are *not* known to be polynomial syzygies.
4. **PO-07** — K6 statement not yet tested under field redefinitions or EOM.
   **Gates all physical and Type IIB reading.**
5. **PO-10** — no induction for any all-orders claim.
6. **PO-04** — external: K6 ↔ Σ₂ needs the authors' change of basis.

## Measured resource data (Apple M1, 8 GiB)

| operation | measured cost |
|---|---|
| full test suite | 171 s clean / 210–233 s under memory pressure |
| static degree-12 certificate, 1 prime | 154 s |
| interacting-flow certificate, 1 prime | ~530 s |
| flow assembly (5 fit + 1 holdout) | seconds |
| closure fixed point (4 primes) | ~95 s; converges in 2 sweeps |
| degree-6 10D order-6 enumeration | 65 s |

**Machine constraint.** 8 GiB total and frequently < 100 MB free. Three pytest
runs were killed by memory pressure during Phase 0. Any long computation must
be checkpointed and run with other applications closed.

**iCloud hazard.** The canonical tree is in an iCloud-synced folder. Do not put
active checkpoints there; judge activity by mtimes, not `ps`.

## Honest scope note

Phases 2–14 are months of work, and several are gated on things a computation
cannot supply: PO-04 needs the authors; Phase 8 needs primary-source Type IIB
expertise; Phase 10 needs a genuinely independent implementation. The program
should be expected to *pause* at those, not to power through them.
