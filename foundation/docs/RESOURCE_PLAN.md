# Resource plan for the remaining phases

All figures below are either **measured** on this machine or **extrapolated
from a measured figure**, and each is labelled. No estimate is invented; where
the cost cannot be bounded from existing data the entry says so.

## 0. Machine and measured primitives

| | |
|---|---|
| CPU | Apple M1, 8 cores |
| RAM | 8 GiB total; frequently **< 100 MB free** in practice |
| BLAS | Accelerate (irrelevant to the integer paths, which are scalar loops) |
| disk | repo ≈ 4 MB tracked + ≈ 3.4 MB certificates |

**Measured primitives** (from run logs and artifact fields, this session):

| operation | measured | source |
|---|---|---|
| full test suite | 168–174 s clean; 210–233 s under memory pressure | four runs |
| degree-12 atlas validation | 679 s | `10d_order12.json:validation_seconds` |
| degree-12 atlas discovery | 830 s, peak RSS **1.18 GB** | `discovery.seconds`, `peak_rss_bytes` |
| static degree-12 cert, 1 prime | 154 s | `static_degree12_32749.json:seconds` |
| interacting-flow cert, 1 prime | 530 s | `interacting_degree12_*.json:seconds` |
| flow assembly (5 fit + holdout) | ~5 s | observed |
| closure fixed point, 4 primes | 96 s; converges in **2 sweeps** | observed |
| degree-10 deficit scan (14 dirs × 4 primes) | ~5 min | observed |
| degree-12 deficit scan (72 dirs, 1 prime + 4-prime confirm) | ~13 min | observed |
| 10D order-6 graph enumeration | 65 s | observed |

**The binding constraint is RAM, not CPU.** Peak RSS of 1.18 GB was already
recorded for atlas discovery, and three pytest runs were killed by memory
pressure during Phase 0. Any step projected above ~1.5 GB should be assumed to
fail on this machine unless sharded.

## 1. Phase-by-phase

### Phase 1 — finish generalized closure (IN PROGRESS)

| item | estimate | basis |
|---|---|---|
| deficit location | **done**, ~18 min | measured |
| basis-change minimality (PO-08), 20 random bases | ~35 min | 20 × 96 s closure, measured |
| reordering falsification (test 5) | ~15 min | one extra scan |
| **intrinsic expressions for 7 directions** | **cannot be estimated** | see below |

M1 verdict: the mechanical parts are comfortable. The intrinsic-expression
work is **not compute-bound** — it is the same task that produced
`K6 = N1050…`, which required choosing a candidate tensor contraction, then
verifying it. Search over candidate contractions is unbounded without a
structural idea; this is a **mathematical-insight step**, and no amount of CPU
substitutes. Honest statement: unknown, plausibly the largest single item in
the whole program.

Embarrassingly parallel: per-direction, per-prime, per-basis-sample.
Checkpoint: one JSON per (direction, candidate) with a verdict.

### Phase 2 — Tr(M⁶)

Route (A), analytic identity: **insight-bound, not compute-bound.**
Route (B), basis-independent formulation: mostly writing, ~1 day.

Extrapolated cost of the brute-force alternative (adding primes): the lift
failed at 15 primes with 29/72 columns over the bound. The *flow* coefficients
lifted at 5 primes, so Tr(M⁶) has genuinely larger height. Doubling the
modulus needs ~15 more primes at 154 s = ~40 min compute — **but there is no
bound saying 30 primes suffice**, and computing that bound (a Hadamard-type
estimate) is itself the useful step. Recommendation: do not buy primes blind.

### Phase 3 — trace–spinor equivalence, real syzygies

Gated on **PO-04, external** — the authors' change of basis is unpublished.
The syzygy search (PO-03) is estimable: enumerating degree-2 monomials in 83
candidates gives ~3.5k columns; an exact modular nullspace at that size is
minutes. Degree-3 is ~100k columns and would need sharding. Feasible on M1 at
degree 2; degree 3 borderline.

Degree 14/16 extension: **not feasible here.** Degree 12 discovery alone was
830 s at 1.18 GB peak with 376 candidates evaluated. Degree 14 grows the
candidate count superexponentially; the README already records order 12 at
~10⁷ candidates for the enumeration problem. Needs a larger machine.

### Phase 4 — geometric structure

Lie brackets of the stress-generated vector fields on a 96-coefficient space:
the objects are the same size as the flow targets already handled, so
**hours, not days**, and it fits in RAM. This is the phase most likely to pay
off per unit compute, because it is where CJ-01 and the coordinate-alignment
observation either become structure or dissolve.

### Phase 5 — all-orders theorem

**Pure insight.** Computer algebra can supply finite certificates for named
lemmas, but PO-10 requires an induction on degree that no computation
produces. Unestimable, and the honest failure mode is "no theorem", which is
an acceptable outcome (option F in the objective).

### Phase 6 — new theory

Depends entirely on Phase 5's outcome. Cannot be estimated before it.

### Phase 7 — causality/stability

Principal symbol and characteristic analysis in D=10 on rational backgrounds:
symbolic work on 10×10 systems, tractable — days of work, light compute.
Explicitly requires **rederiving** the higher-form analogue rather than
importing 4D nonlinear-electrodynamics formulas.

### Phase 8 — Type IIB

**Gated by PO-07 and by expertise, not compute.** Requires primary-source
reading with exact conventions, derivative orders, and field-redefinition
bookkeeping. This is where fabrication risk is highest and where I should stop
and involve a person.

### Phase 9 — algorithm generalisation

A third representation at low degree is cheap (the 6D benchmark runs in 27 s).
Choosing a *scientifically informative* representation is the hard part.
Estimate: days, low compute.

### Phase 10 — clean-room reproduction

**Cannot be discharged by me alone in the same session.** A second Python
implementation written from my own memory of the first is not independent, and
certifying it as such would be false. Needs either a different stack driven
from published definitions only, or a different person/agent with no access to
this code. Compute cost is small; the integrity requirement is the obstacle.

### Phases 11–14

Review, mentor package, manuscripts, release: writing and literature work.
Light compute. Phase 13 is blocked by design until 1–12 land.

## 2. Summary table

| phase | compute-bound? | fits on M1? | main blocker |
|---|---|---|---|
| 1 mechanical | yes | yes | — |
| 1 intrinsic | **no** | n/a | mathematical insight |
| 2 | no | yes | insight; height bound |
| 3 syzygies deg 2 | yes | yes | — |
| 3 spinor / deg 14+ | — | **no** | PO-04 external; hardware |
| 4 | yes | yes | — |
| 5 | **no** | n/a | induction (PO-10) |
| 6 | — | — | depends on 5 |
| 7 | no | yes | rederivation effort |
| 8 | **no** | n/a | PO-07 + human expertise |
| 9 | yes | yes | choosing the case study |
| 10 | no | yes | independence integrity |
| 11–14 | no | yes | prerequisites |

## 3. Checkpoint and stopping strategy

- One immutable JSON per (task, prime, sample); atomic write via temp +
  fsync + rename, as the existing scripts already do.
- **Never** place active checkpoints in the iCloud-synced tree. Completed
  immutable shards only.
- Stop condition for every scan: exhaustive over the basis, or an explicitly
  stated searched range if not exhaustive (standard 15).
- Pilot before any unbounded run: time one unit, multiply, compare against the
  1.5 GB RAM ceiling.
- Memory guard: if free RAM < 200 MB, do not launch; three pytest kills in
  Phase 0 were all attributable to this.

## 4. What more compute cannot buy

Recorded explicitly, because the program's remaining risk is concentrated
here: **Phase 1's intrinsic expressions, Phase 2's analytic identity, Phase 5's
induction, Phase 8's Type IIB reading, and Phase 10's independence** are not
compute problems. Four of the five ultimate pillars depend on at least one of
them. A session that reports progress by burning CPU on the estimable parts
while these stay open is not advancing the program.
