# Research decision log

Append-only. Each entry: date, decision, alternatives, reason.

## 2026-07-30 — Phase 0

**D-001 — Reproduce in a non-synced directory.**
The canonical tree is under an iCloud-synced path. During the prior session
this produced writes with no local process and a torn copy that ran stale
bytecode, yielding one spurious test failure. Rejected: reproducing in place.
Chosen: fresh clone in `~/Downloads`, caches excluded, activity judged by
mtimes over a ≥60 s window.

**D-002 — Treat the reported baseline as claims, not axioms.**
Every number in the objective was re-derived from committed artifacts before
being entered in the ledger. This surfaced that the objective's "83 primitive
candidates" and the "1,2,7,14,72" dimensions are different quantities; both
retained as separate claims rather than reconciled by assumption.

**D-003 — Discharge PO-06 immediately rather than in Phase 1.**
Two independent holdouts cost minutes because six certificates already
existed. Doing it in Phase 0 means C-FLOW-01 enters Phase 1 already at its
strongest supportable status.

**D-004 — Do not upgrade any MOD-CERT claim on prime agreement alone.**
Fifteen primes agreeing is strong evidence, not proof. PO-09 records the
exceptional-prime gap explicitly rather than treating agreement as closure.

**D-005 — Gate all physical/Type IIB reading behind PO-07.**
The K6 transport result is off-shell and convention-fixed. Until the action of
field redefinitions and EOM on q6 is known, no physical consequence may be
drawn. This is recorded as a hard gate, not a caveat.

**D-006 — Phase 1 attacks the deficits via seeding, and will say so.**
The available certificates support seed-enlargement closure exactly. Adding a
*generator* is a different computation requiring rows that do not exist. Phase 1
will compute what it can compute exactly and label it precisely, rather than
blur the two (C-MIN-03).
