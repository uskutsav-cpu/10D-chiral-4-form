# Claim ledger

Every scientific assertion this repository is permitted to make, with its
evidence class, its allowed manuscript wording, and the stronger wording that
is **forbidden**.

Baseline commit: `3ed32805b38ce34216b34888f6539e3538e90fb9`
Branch: `research/maximal-chiral-four-form-program`
Last updated: Phase 1 Step 1 (quotient spaces constructed)

## Evidence classes

| class | meaning |
|---|---|
| `EXACT-THM` | proved by hand, no computer step |
| `EXACT-CA-THM` | proved, with a finite computer-checked certificate |
| `FINITE-ORDER` | exact but only up to a stated degree cutoff |
| `MOD-CERT` | modular (finite-field) certificate; characteristic-zero inference needs stated hypotheses |
| `RAT-RECON` | rational reconstruction, CRT-bounded, holdout-validated |
| `NUM-EV` | high-precision numerical evidence only |
| `CONJ` | conjecture |
| `HEUR` | heuristic |
| `BLOCKED` | external blocker |

**Global caveat applying to every `MOD-CERT` row.** A rank computed mod p is a
lower bound for the characteristic-zero rank, and equals it unless p divides a
relevant minor. We use ≥4 primes and observe identical ranks; that makes an
exceptional-prime coincidence unlikely but does **not** prove it. No row below
may be upgraded to `EXACT-THM` on the strength of prime agreement alone.

---

## C-ATLAS-01 — homogeneous invariant dimensions through degree 12

- **Wording**: "The spaces of homogeneous Lorentz invariants of a self-dual
  five-form in D=10 have dimensions 1, 2, 7, 14, 72 at degrees 4, 6, 8, 10, 12."
- **Assumptions**: generic self-dual five-form; the graph-contraction candidate
  set is complete at each degree; exact canonicalisation.
- **Range**: degrees 4–12.
- **Class**: `MOD-CERT` + `FINITE-ORDER`
- **Sources**: `results/10d_order12.json`, `src/sdinv/graphs.py`
- **Verification**: 15 primes, 4 samples/prime, identical ranks.
- **Independent verification**: **NOT DONE** (clean-room, Phase 10).
- **Permitted**: "exact ranks over several finite fields, identical at every
  prime tested."
- **Forbidden**: "proved dimensions"; "the invariant ring has these dimensions"
  (ring vs. space of homogeneous polynomials is not established here).
- **Caveats**: completeness of the candidate set relies on the multigraph
  argument, which is sound, plus exact canonicalisation for n ≤ 6 and pynauty
  above; see C-ATLAS-05.

## C-ATLAS-02 — degree-12 decomposition

- **Wording**: "At degree 12 there are 10 lower-degree product directions and
  62 connected primitive directions, total value-space rank 72."
- **Class**: `MOD-CERT` + `FINITE-ORDER`
- **Independent verification**: **NOT DONE**
- **Forbidden**: calling the 62 a basis of primitives of the invariant *ring*.

## C-ATLAS-03 — cumulative generic functional rank 81

- **Wording**: "The cumulative generic Jacobian rank through degree 12 is 81."
- **Assumptions**: genericity of the sample point; Jacobian rank at a generic
  point equals the functional rank.
- **Class**: `MOD-CERT`
- **Permitted**: "the cumulative Jacobian rank saturates the expected count 81
  = 126 − 45 at degree 12."
- **Forbidden**: "we have found all 81 invariants"; "the atlas is complete."
  Saturating a rank is not exhibiting a generating set.
- **Caveats**: 126 − 45 = 81 is a dimension count for a generic orbit, not a
  proof that these 81 generate the invariant ring.
- **Independent verification**: **DONE**, by a second implementation in spinor
  variables that shares no code with the trace side. See
  `results/rank81/certificate.json`: the exact analytic Jacobian of the complete
  83-candidate selection has modular rank 81 with zero evaluation errors and
  zero zero-rows, and `results/rank81/minor81_certificate.json` exhibits an
  explicit 81×81 minor with non-vanishing determinant under two independent
  determinant routines.
- **Strengthening this permits**: the bound `rank_Q ≥ 81` becomes
  *unconditional* rather than probabilistic, because the coordinate basis is
  integral and the Jacobian is therefore the reduction of an integer matrix. The
  genericity assumption above is still needed for the *upper* direction, which
  remains analytic. Nothing here licenses "proved 81 computationally".

## C-ATLAS-04 — I12_61 and I12_62: polynomially independent, functionally redundant

- **Wording**: "I12_61 and I12_62 are **polynomially independent** degree-12
  invariants — they are members of the 72-element basis whose polynomial rank
  is 72 — but they are **not new generic functional directions**: adjoining
  them leaves the cumulative generic functional rank at 81."
- **Assumptions**: genericity of the evaluation point; Jacobian rank at a
  generic point equals functional rank.
- **Range**: degree 12.
- **Class**: `MOD-CERT`
- **Sources**: `results/10d_order12.json` →
  `discovery.degree12_polynomial_rank = 72`,
  `discovery.cumulative_functional_rank = 81`,
  `discovery.functional_dependencies = ["I12_61","I12_62"]`
- **Holdout status**: covered by the atlas's three-prime / four-sample
  validation; no dedicated holdout for this sub-claim.
- **Independent verification**: **NOT DONE**.
- **Permitted**: "polynomially independent but functionally redundant at a
  generic point."
- **Forbidden**: **"these are polynomial syzygies"** — a Jacobian dependence
  is not a polynomial identity (PO-03). Equally forbidden: "these are
  redundant invariants" without the qualifier *functionally*, since they are
  linearly independent polynomials and removing them would shrink the
  degree-12 space from 72 to 70.
- **Caveats**: PO-03. Note also C-MIN-02 — both are in the degree-12
  unreachable set, which is conjecture CJ-01, not a theorem.

## C-SCOPE-01 — finite-order versus all-orders language

- **Wording**: "Every quantitative result in this repository is established at
  field degree ≤ 12 and nowhere above it."
- **Class**: `EXACT-THM` (a statement about our own evidence)
- **Permitted**: "verified through degree 12."
- **Forbidden**: "all orders", "in general", "for every degree", "the pattern
  continues", or any unqualified present tense that implies unbounded degree.
  No induction on degree exists (PO-10).
- **Caveats**: the single exception is C-FLOW-02, which is all-orders **in λ**
  at **fixed field degree 6**. Those are different variables and the
  distinction must never be blurred.

## C-SCOPE-02 — formal, local, generic, global

- **Wording**: "Flow statements are **formal** in λ (power series, no
  convergence proved), **generic** in the field configuration (evaluated at
  deterministic generic samples), and **local** in theory space (no global
  orbit statement)."
- **Class**: `EXACT-THM` (scope declaration)
- **Assumptions**: samples used are generic; non-generic strata are untested.
- **Permitted**: "a formal power-series flow, at generic configurations."
- **Forbidden**: "convergent", "analytic in λ", "globally defined", "for all
  field configurations", "the orbit is". Special/singular strata have not been
  examined at all.
- **Caveats**: falsification test 7 (special non-generic configurations) is
  **not yet run**. Until it is, every result is a generic-point statement.

## C-SCOPE-03 — Type IIB relevance versus a Type IIB result

- **Wording**: "This work concerns the algebraic invariant theory of a
  self-dual five-form. **No Type IIB result exists.**"
- **Class**: `EXACT-THM` (statement of absence)
- **Permitted**: nothing beyond noting that the same field content appears in
  Type IIB supergravity, explicitly flagged as motivation, not result.
- **Forbidden**: any claim that a Type IIB correction is constrained,
  excluded, matched, or predicted. Sharing a self-dual five-form is not a
  connection: Type IIB corrections involve curvature, derivatives, other
  fields, self-duality prescriptions and field-redefinition ambiguities, none
  of which are modelled here.
- **Caveats**: **PO-07 gates this absolutely.** The K6 statement is off-shell
  and convention-fixed; until its behaviour under field redefinitions and the
  equations of motion is known, no physical reading is permitted at all.

## C-ATLAS-05 — candidate-set completeness above degree 6

- **Wording**: "Candidate multigraphs are deduplicated by exact canonical form
  for n ≤ 6 and by pynauty above."
- **Class**: `EXACT-CA-THM` (given pynauty correctness)
- **Caveats**: PO-01 **discharged** in Phase 0 — `_canonical_wl` no longer
  exists, `canonical()` uses pynauty's exact certificate or raises above six
  vertices, and order 8 comes from nauty `geng | multig`. The historical WL
  collision (49 classes → 39 keys at order 6) cannot occur at this commit.

## C-SEXTIC-01/02/03 — intrinsic sextic basis

- **Wording**: "J6 = Tr(M³) = (32/3)·I6_1; K6 = −(1/1125)·I6_1 + (3/125)·I6_2;
  the change of basis has determinant 32/125 and quotient coordinate
  q6(c1,c2) = (125/3)·c2."
- **Class**: `RAT-RECON` + `EXACT-CA-THM`
- **Sources**: `results/stress_flow/change_of_basis/sextic_intrinsic.json`
- **Verification**: primes 32749/32719/32693, 4 samples each; an independently
  implemented dense N^(1050) cubic contraction agrees exactly on every sample.
- **Independent verification**: partial — the dense contraction is a genuine
  second implementation *within* this repository, not clean-room.
- **Permitted**: "exact change of basis, verified against an independent dense
  tensor contraction."
- **Forbidden**: claiming K6 is the published spinor invariant Σ₂. The source
  proves (Σ₁,Σ₂) is a sextic basis but does **not** publish the change of basis;
  our identification is inferred. See PO-04.

## C-STATIC-01 — static free-stress span

- **Wording**: "The static free-stress span has dimensions 1, 1, 2, 2, 4 at
  degrees 4, 6, 8, 10, 12."
- **Class**: `MOD-CERT` + `FINITE-ORDER`
- **Forbidden**: presenting this as the set a flow can reach. See C-CLOSURE-02.

## C-TRM6-01 — Tr(M⁶) coordinates are NOT certified

- **Wording**: "The degree-12 free-stress span is 4-dimensional. Tr(M⁶) is
  retained as a stress-adapted basis element; its coordinates in the 72-element
  graph basis are not certified."
- **Class**: rank is `MOD-CERT`; coordinates are `BLOCKED`
- **Evidence**: 15 primes, CRT modulus ≈ 5.2e67; **29 of 72 columns** exceed
  the uniqueness bound.
- **Permitted**: "the rank is certified; no certified rational coordinate
  vector exists at the current prime budget."
- **Forbidden**: quoting any rational coordinate vector for Tr(M⁶), including
  the entries under `rational_lift_audit.successfully_reconstructed_columns`
  (CRT terminating is not CRT being correct).
- **Caveats**: PO-05; Phase 2 must either find an analytic reduction or prove
  the theorems basis-independently.

## C-FLOW-01 — new-forcing dimensions

- **Wording**: "The new-forcing span of the interacting stress flow has
  dimensions 1, 1, 3, 5, 21 against full 1, 2, 7, 14, 72, leaving obstruction
  dimensions 0, 1, 4, 9, 51."
- **Assumptions**: analytic polynomial interactions; scalar trace generators;
  `Tr(tau)` excluded from forcing because it is the homogeneity operator.
- **Class**: `RAT-RECON` + `FINITE-ORDER`
- **Sources**: `results/stress_flow/interacting_flow_equations.json`
- **Verification**: fit primes 32749/32719/32693/32771/32713; independent
  holdout 32717; `all_modular_and_rational_holdouts_passed: true`.
- **Independent verification**: **NOT DONE**
- **Caveats**: PO-06 **discharged** in Phase 0 — two independent holdouts
  (32713 and 32717) both pass, and the two fit sets give 192/192 identical
  reconstructions.

## C-FLOW-02 — the K6 transport equation

- **Wording**: "The degree-six component of the flow obeys
  dq6/dλ = 10·(6−2)·a(λ)·q6 = 40·a(λ)·q6, where a is the coefficient of the
  Tr(τ) generator."
- **Assumptions**: the degree-6 sector closes (no higher-degree feedback);
  the generator set is complete at degree 6.
- **Range**: degree 6; **all orders in λ as a formal series**.
- **Class**: `EXACT-CA-THM`
- **Supporting facts now verified**: the generator enumeration is complete —
  18 expected multisets, 18 present, additivity of leading degrees holds for
  every generator (verified Phase 0). Only `tr_tau`, `tr_tau2`, `tr_tau3` have
  leading degree ≤ 6, and all three produce rows.
- **Permitted**: "an exact linear homogeneous transport equation, valid to all
  orders in λ as a formal power series at field degree six."
- **Forbidden**: "valid to all orders" without the words *formal* and *at
  degree six*. Convergence in λ is not established.

## C-FLOW-03 — K6 is transported, never created

- **Wording**: "q6(λ) = q6(0)·exp(40∫a). Hence q6(0)=0 ⟹ q6≡0, and q6(0)≠0 ⟹
  q6≠0 for all λ."
- **Class**: `EXACT-CA-THM` (immediate from C-FLOW-02)
- **Permitted**: "K6 is transported, never created."
- **Forbidden**: **"K6 must vanish in every pure stress flow"** — false for a
  seeded flow. Also forbidden: calling this a universal obstruction, or
  inferring it from the static span alone.
- **Caveats**: this is an off-shell, field-redefinition-fixed statement.
  Whether it survives field redefinitions and use of equations of motion is
  PO-07 and is **not** established.

## C-CLOSURE-01 — free-seed iterative closure

- **Wording**: "The iterative closure from the free seed has dimensions
  1, 1, 3, 11, 67."
- **Class**: `MOD-CERT` + `FINITE-ORDER`
- **Sources**: `results/stress_flow/closure_and_minimality.json`
- **Caveats**: computed under *seed* enlargement, which is not the same as
  adding a flow generator; see C-MIN-03.

## C-CLOSURE-02 — closure strictly exceeds the static span

- **Wording**: "The reachable set (1,1,3,11,67) strictly exceeds the static
  span (1,1,2,2,4) at degrees 8, 10 and 12."
- **Class**: `MOD-CERT`
- **Permitted**: "a no-go argument built on the static span alone is invalid."
- **Forbidden**: any claim that the static complement is a dynamical
  obstruction.

## C-MIN-01 — minimal completion through degree 8

- **Wording**: "Through degree 8 the completion requires exactly five
  directions — K6, I8_3, I8_4, I8_5, I8_6 — each non-redundant under removal."
- **Class**: `MOD-CERT` + `FINITE-ORDER`
- **Forbidden**: claiming minimality under *arbitrary* basis change; only
  removal-minimality within the fixed basis is shown. PO-08.

## C-MIN-02 — residual deficit — RESOLVED (Phase 1)

- **Wording**: "The degree-10 deficit of 3 is spanned by I10_6, I10_7, I10_12;
  the degree-12 deficit of 4 by I12_59, I12_60, I12_61, I12_62. Adjoining all
  seven closes both degrees. Each is non-redundant under removal."
- **Class**: `MOD-CERT` + `FINITE-ORDER`
- **Range**: degrees 10, 12.
- **Sources**: `results/generalized_flow/degree10_missing_directions.json`,
  `results/generalized_flow/degree12_missing_directions.json`
- **Verification**: exhaustive scan over every basis direction at each degree;
  confirmed on primes 32749, 32719, 32717, 32693.
- **Independent verification**: **NOT DONE**.
- **Permitted**: "the deficits are located exactly and each element is
  necessary."
- **Forbidden**: calling these directions *intrinsic* — they are graph labels
  in one basis. Also forbidden: presenting them as flow **generators**
  (C-MIN-03).
- **Caveats**: Phase 1 gate NOT met pending intrinsic expressions; PO-08
  still PARTIAL, but its two halves are now separated:
  - *cardinality* minimality is **basis-independent and proved**. Any set that
    closes degree `d` has at least `dim Q_d` elements, because the quotient map
    is linear and a span of `k` vectors has dimension at most `k`; the exhibited
    sets have exactly 3 and 4 elements against deficits of 3 and 4, so they are
    of minimum cardinality **in any basis**. See PO-08 in
    `docs/PROOF_OBLIGATIONS.md` and `tests/test_quotient_cardinality_bound.py`.
  - *removal* minimality — that no single exhibited element may be dropped —
    remains shown only in the fixed basis and under the permutation subgroup,
    and is **still open** under general `GL`.
- **Newly permitted by the above**: "no fewer than three degree-10 and four
  degree-12 directions can close the respective deficits, in any basis."
- **Still forbidden**: "this is the unique minimal set", and any basis-free
  claim about *which* directions are required.

## C-MIN-04 — the closure is coordinate-aligned

- **Wording**: "At degree 12 the 68-dimensional closure contains exactly 68 of
  the 72 basis vectors; at degree 10, 11 of 14."
- **Class**: `MOD-CERT`
- **Permitted**: "the closure coincides with a coordinate subspace of the graph
  basis, which a generic subspace of that dimension would not."
- **Forbidden**: treating the coordinate alignment as basis-independent. Under
  a generic change of basis it disappears; only the quotient dimension is
  intrinsic.

## CJ-01 — CONJECTURE: functional dependence implies unreachability

- **Wording**: "Every functionally dependent degree-12 candidate is dynamically
  unreachable."
- **Class**: `CONJ`
- **Evidence**: 2 of 2 (`I12_61`, `I12_62` are both the atlas's recorded
  `functional_dependencies` and members of the unreachable set). The converse
  is **false**: `I12_59`, `I12_60` are unreachable and functionally
  independent.
- **Forbidden**: any manuscript appearance. Two data points, no mechanism.
- **Target**: Phase 4 must explain it or expose it as coincidence.

## C-MIN-03 — seeding vs. generators

- **Wording**: "The closure analysis enlarges the seed, not the generator set."
- **Class**: `EXACT-THM` (definitional)
- **Forbidden**: presenting seeding results as `f(T,S,λ)` generator results.

## C-GEN-01 — generator enumeration complete

- **Wording**: "The 18 trace generators are exactly the multisets of Tr(τ^k)
  whose leading field degrees sum to at most 12."
- **Class**: `EXACT-CA-THM` (verified Phase 0: 18 expected = 18 present,
  none missing, none extra, additivity holds for all)
- **Caveats**: relies on `leading_field_degree` being correct per generator,
  which is itself computed; PO-02.

## C-REPRO-01 — test suite

- **Wording**: "96 tests pass at commit 3ed3280."
- **Class**: `EXACT-CA-THM` (software fact)
- **Status**: Phase 0 fresh-clone reproduction — see
  `docs/baseline_reproduction_maximal_program.md`.

---

## Claims that must be WEAKENED or are FORBIDDEN outright

| # | forbidden statement | why |
|---|---|---|
| F-1 | "K6 must vanish in every pure stress flow" | false for a seeded flow (C-FLOW-03) |
| F-2 | "I12_61/I12_62 are syzygies" | Jacobian dependence ≠ polynomial identity (C-ATLAS-04) |
| F-3 | any rational coordinate vector for Tr(M⁶) | uncertified (C-TRM6-01) |
| F-4 | "we have all 81 invariants" | rank saturation ≠ generating set (C-ATLAS-03) |
| F-5 | "all-orders theorem" from degree-12 data | no induction established |
| F-6 | "relevant to Type IIB" as a result | not a result; Phase 8 gate |
| F-7 | "independently reproduced" | no clean-room implementation yet (Phase 10) |
| F-8 | "K6 is the published spinor invariant Σ₂" | change of basis not in the source |
| F-9 | "the static complement is an obstruction" | contradicted by C-CLOSURE-02 |
| F-10 | "new nonlinear theory" | no theory constructed yet (Phase 6) |


## C-QUOT-01 — the intrinsic quotient spaces

- **Wording**: "The generalized-flow quotients have dimensions dim Q_10 = 3 and
  dim Q_12 = 4, where Q_d = A_d / D_d and D_d is the reachable seed closure."
- **Assumptions**: seed closure as computed by the fixed point; certificate
  target rows complete at each degree.
- **Range**: degrees 10 and 12.
- **Class**: `MOD-CERT`
- **Sources**: `results/generalized_flow/quotient_degree10.json`,
  `results/generalized_flow/quotient_degree12.json`
- **Holdout status**: identical at all six available primes (32749, 32719,
  32717, 32693, 32771, 32713).
- **Independent verification**: **NOT DONE**.
- **Permitted**: "the quotient dimensions are 3 and 4"; "dim Q_d is the
  basis-independent content."
- **Forbidden**: calling `I10_6, I10_7, I10_12` or `I12_59..I12_62`
  **intrinsic**. They are coordinate labels of one graph basis, selected by an
  echelon form. No intrinsic tensor expression exists for any of them.
- **Caveats**: PO-09 (modular); PO-08 (GL basis change partial).

## C-QUOT-02 — the quotient classes are primitive, not products

- **Wording**: "No quotient direction at degree 10 or 12 is a product of
  lower-degree invariants: every product entry of both bases lies inside the
  closure."
- **Range**: degrees 10 and 12.
- **Class**: `MOD-CERT`
- **Evidence**: at degree 10 both product entries (`I4_1*I6_1`, `I4_1*I6_2`)
  and at degree 12 all ten (`I4_1^3`, `I6_1^2`, `I6_1*I6_2`, `I6_2^2`,
  `I4_1*I8_1..6`) are in D_d; none is in a missing set.
- **Permitted**: "the quotient classes are genuinely primitive."
- **Forbidden**: inferring from this what the classes *are*. It is an exclusion,
  not an identification.

## C-QUOT-03 — SUPERSEDED, and the error is recorded

The earlier wording — "intrinsic quotient rank 0 of 3 and 0 of 4" — was
**wrong**. It conflated *no compact named formula* (true) with *no intrinsic
representative* (false). A verified contraction graph already defines a
coordinate-independent Lorentz scalar; the missing step was the deterministic
translation into explicit index form, not new mathematics. Superseded by
C-LEVELA-01. Retained here because a ledger that silently deletes its own
mistakes is not an audit trail.

## C-LEVELA-01 — Level-A intrinsic representatives, all seven classes

- **Wording**: "Each of the seven generalized-flow quotient classes has an
  explicit Einstein-index representative
  `I = prod_v F_{s(v,1..5)} x prod_pairs eta^{s_a s_b}`, with 25 metric
  factors at degree 10 and 30 at degree 12."
- **Assumptions**: the committed contraction graphs are correct (C-ATLAS-05,
  PO-01 discharged).
- **Range**: degrees 10 and 12.
- **Class**: `EXACT-CA-THM`
- **Sources**: `results/intrinsic_candidates/explicit_F_contractions.json`,
  `src/sdinv/graph_to_tensor.py`
- **Holdout status**: 3 primes x 6 samples (4 fitting + 2 fresh) per class =
  **126 dense-vs-graph checks, 126 agreements**, no sign discrepancy.
  Homogeneity `F -> cF ~ c^degree` verified on every class.
- **Independent verification**: the dense evaluator is a second code path
  within this repository (einsum built from the dummy-index assignment, metric
  applied by raising one slot). Not clean-room.
- **Permitted**: "explicit, coordinate-independent intrinsic representatives";
  "Level A complete for all seven classes."
- **Forbidden**: calling these **compact**, **canonical**, **simplest**, or
  **explanatory**. They are 10- and 12-fold contractions. Also forbidden:
  claiming the quotients are *understood* — that needs Level B.
- **Caveats**: intrinsic IDs `Q10_A/B/C`, `Q12_A/B/C/D` are provisional; their
  correspondence to graph labels is recorded, not an identification.

## C-LEVELB-01 — NEGATIVE: no M/N expression for any class

- **Wording**: "Level B (expression in M^(54), N^(1050), N^(4125)) is derived
  for **0 of 7** classes. Level C is not attempted."
- **Class**: bounded negative result
- **Families searched**: products of lower intrinsics — excluded by C-QUOT-02
  (all product basis entries lie inside the closure, so the classes are
  genuinely primitive). **Not searched**: quadratic-block multisets
  `M^a (N1050)^b (N4125)^c` with a+b+c = 5 and 6; epsilon/parity-odd channels;
  flow-generated structures. Static-span-versus-D_d relation not computed.
- **Next input needed**: the quadratic-block contraction generator (Stage 2),
  or the reverse route of Stage 5 — pairing the F's of a known Level-A graph
  into blocks and reducing.
- **Forbidden**: any compact tensor formula not validated against fitting
  primes, holdout primes and fresh samples.


## C-MONLY-01 — the M-only quotient rank at degree 10 is zero

- **Wording**: "The complete family of M-trace monomials at five-form degree 10
  — `tr(M^5)` and `tr(M^3)tr(M^2)` — lies entirely inside the reachable closure
  D_10. The M-only quotient rank is 0."
- **Assumptions**: scalars built purely from M are generated by traces of its
  powers; `tr(M) = 0` removes parts of size 1; generic samples.
- **Range**: degree 10.
- **Class**: `MOD-CERT`
- **Sources**: `results/intrinsic_candidates/M_only_quotient_test_deg10.json`,
  `scripts/test_M_only_quotients.py`
- **Holdout status**: primes 32749, 32719, 32693, 32771, 32713 (fit) and 32717
  (holdout); identical at all six.
- **Non-vacuity**: both monomials have `status = "solved"` — they lie in the
  atlas span — and quotient vectors exactly `[0,0,0]`. A monomial that merely
  failed to solve would also yield rank 0 for an uninformative reason; the
  status field distinguishes these and is asserted in the tests.
- **Independent verification**: **NOT DONE**.
- **Permitted**: "every degree-10 quotient class requires at least one N
  channel"; "M-only structures do not reach the quotient."
- **Forbidden**: extending this to degree 12, which is **not computed**. Also
  forbidden: concluding anything about which N channel is needed — this is an
  exclusion, not an identification.
- **Consequence**: combined with C-QUOT-02 (products excluded), the block
  multiset enumeration `M^a (N1050)^b (N4125)^c` with `a+b+c = 5` narrows to
  `b + c >= 1`.


## C-P12-01 — the published degree-12 structures have zero image in Q12

- **Wording (the only permitted form)**: "The three degree-12 structures
  displayed in equation (4.25) lie inside the computed reachable closure D12
  and have zero image in Q12 under the stated formal flow definition."
- **Assumptions**: the seed closure D12 as computed by the fixed point; the
  certificate target rows complete at degree 12; generic samples.
- **Range**: degree 12.
- **Scope**: generic, formal, modular.
- **Class**: `MOD-CERT`
- **Sources**: `src/sdinv/published_degree12_invariants.py`,
  `results/intrinsic_candidates/published_degree12_map.json`
  (sha256 `9a8f4e32627ee5bef06324169dc823a4…`)
- **Primes**: 32749 and 32717, independently, both giving rank 0 with all
  vectors of length 4 identically zero.
- **Non-vacuity**: every structure has `status = "solved"`. A structure that
  failed to solve would also contribute rank 0 for an uninformative reason.
- **Independent verification**: **NOT DONE**.
- **Forbidden stronger wording**, each an available misreading:
  1. that the complete M/N degree-12 space has zero image in Q12 — only three
     structures were tested;
  2. that Q12 has no compact tensor representation — Q12 has dimension 4;
  3. that the three invariants are redundant — they are nonzero elements of
     the degree-12 space that happen to lie in D12;
  4. that no generalized flow can generate Q12 — a different question, about
     generators rather than these three scalars.
- **Caveats**: `tr(M^6)` has no certified rational graph-basis coordinates
  (29 of 72 columns exceed the CRT bound at 15 primes); only the modular span
  and quotient statements are certified, which suffices for a rank claim.


## C-P10-BUG-01 — the shipped P10_07 was not a Lorentz scalar

- **Status**: DEFECT FOUND AND FIXED. Recorded because the *detection route*
  matters more than the bug.
- **What was wrong**: `p10_07_n1050_mm_n1050_n1050` raised all six axes on both
  inner `N^(1050)` factors, so the three alpha edges each joined two raised
  slots and contracted with `delta` rather than `eta`.
- **Observable**: `not_in_atlas_span` on all six primes. A consistent failure
  across independent moduli is structural, never modular bad luck.
- **Confirmation**: rotation invariant (22769 -> 22769) but boost violating
  (22769 -> 19807) at 32749, and likewise at 32719. That asymmetry is the
  signature of a metric misplacement: `delta` and `eta` agree on the spatial
  block, so a rotation cannot see it.
- **Fix**: `mm = mixed @ mixed` is `(MM)_{a}{}^{b}` — slot 0 DOWN, slot 1 UP —
  the opposite of what the docstring asserted. Raise `(5,)` on the outer `N`,
  `(0,1,2)` and `(0,1,3,4,5)` on the inner pair.
- **Independent check that the fix is not accidental**: moving the metric to
  the other end of every alpha edge gives bit-identical values (10384 at
  32749, 25537 at 32719). A genuine tensor contraction cannot depend on which
  end carries the metric.
- **Why it survived review**: no boost test existed for any degree-10
  candidate. The suite tested homogeneity, which catches int64 overflow, and
  nothing that constrains index placement.
- **Permitted**: "a metric-placement defect in P10_07 was found by the atlas
  projection and confirmed by a boost test."
- **Forbidden**: treating the earlier `Q10 rank 0` reading as having covered
  P10_07 — it never contributed a vector, because it never solved.
- **Regression cover**: `test_every_published_candidate_is_boost_invariant`,
  `test_p10_07_alpha_edge_placement_is_free`,
  `test_the_original_p10_07_placement_really_was_broken`.


## C-P10-METHOD-01 — index placement is derived, not read

- **Wording**: "Up/down index placement for the equation-(4.24) candidates is
  determined by the rule that every contracted edge carries exactly one raised
  end, and validated by a boost test; it is not read from the PDF text
  stream."
- **Justification**: PDF extraction returns stacked super/subscripts in
  glyph-position order, which does not reliably distinguish `M_{a}{}^{b}` from
  `M^{a}{}_{b}`. The edge rule has a unique answer up to which end carries the
  metric, and that choice provably cannot change the value.
- **Validation**: the rule independently reproduces the placement that was
  verified for P10_07 by boost test.
- **Class**: `METHOD`
- **Forbidden**: claiming the *bracket colour* is similarly recoverable. It is
  not — see AMB-01 and AMB-02 in `PUBLISHED_DEGREE10_INDEX_AUDIT.md`.


## C-P10-01 — the twelve published degree-10 candidates span Q10

- **Wording (the only permitted form)**: "The twelve degree-10 structures
  displayed in equation (4.24) span the quotient Q10, of dimension 3, under the
  stated formal flow definition; eight of the twelve lie inside the reachable
  closure D10."
- **Class**: `MOD-CERT`
- **Primes**: 32749 (fit) and 32717 (holdout), independently, both rank 3 of 3
  and published atlas rank 12 of 14.
- **Non-vacuity**: all thirty-two projections have `status = "solved"`. A rank
  assembled from unsolved rows would be uninformative.
- **Positive control**: the Level-A representatives recover rank 3 of 3 through
  the identical projector at six primes, so a zero row is a statement about the
  candidate, not about the pipeline.
- **Which candidates reach Q10**: P10_09, P10_10, P10_11, P10_12 only.
- **SUPERSEDES** the earlier reading in this branch that every published
  degree-10 candidate had zero image in Q10. That was drawn from the five
  simplest candidates, selected for ease of implementation — which correlates
  with structural simplicity, which is what decides whether a candidate reaches
  the quotient. A null result over a subset chosen for convenience is not
  evidence about the whole set.
- **Forbidden**: extending to Q12; claiming independence of the AMB-01/AMB-02
  readings; claiming five fitting primes — two primes were used.


## C-P10-02 — a preferred Level-B basis for Q10, conditional on source readings

- **Wording**: "{P10_10, P10_11, P10_12} is a three-element Level-B basis for
  Q10, preferred under the simplicity criterion documented in
  `scripts/build_Q10_levelB_basis.py` and minimal among the twelve published
  candidates, conditional on the AMB-02 reading recorded for P10_11."
- **Removal minimality**: P10_10 appears in every independent triple and cannot
  be dropped — it is the only candidate reaching the third quotient direction.
  Established at both primes.
- **Level-A map**: Q10_A, Q10_B, Q10_C expressed in the basis by exact modular
  solve; all six round-trips verified.
- **Conditionality, stated rather than buried**: NO independent triple is free
  of source ambiguity. P10_10 is forced and the remaining two slots must draw
  at least one member from {P10_09, P10_11}, both of which move under their
  unresolved bracket reading. Only P10_10 and P10_12 are unconditional.
- **Forbidden**: "canonical" — no class larger than the twelve was enumerated;
  "rational" — the certification is modular, with no CRT reconstruction.


## C-REV-01 — Q10 recovered by a formula-independent bounded reverse search

- **Wording (the only permitted form)**: "A formula-independent bounded reverse
  search independently recovered the full three-dimensional quotient Q10."
- **Class**: `MOD-CERT`
- **Independence**: enforced, not promised.
  `src/sdinv/reverse_block_decomposition.py` does not import
  `published_degree10_invariants`, and
  `test_reverse_engine_does_not_import_the_published_formulas` parses its AST
  and fails if it ever does. The published basis is loaded only in the
  post-hoc comparison script, after the search has finished.
- **Result**: rank 3 of 3, from 588 evaluated candidates over 22 samples at
  prime 32749; validated on holdout prime 32717 and on a fresh seed base.
- **The recovered basis is NOT the published one.** All three directions come
  from the `N^(4125)`x5 sector — pure `N^(4125)` contractions with no
  `N^(1050)` and no `M`. Every published Level-B basis element is
  `N^(1050)`-based. This is a stronger outcome than rediscovery: Q10 is
  reachable from a block sector the published list does not use.
- **Span equality**: both spans have rank 3 inside a 3-dimensional quotient, so
  each *is* Q10 and they are equal; the union rank is the direct check.
- **Topology provenance**: each recovered candidate is regenerated from its
  recorded einsum by re-enumerating its sector, not deserialised, so the
  recorded specification is proved producible by the declared search.
- **FORBIDDEN wording**, each an available misreading:
  1. "Complete enumeration of every M/N contraction." **5 of 21 sectors are
     capped** at 30 000 raw topologies and are not exhausted; even the
     exhausted sectors were sampled at 40 candidates for evaluation.
  2. "The reverse search proves the published basis minimal." It does not
     address minimality at all.
  3. "The reverse search reproduces the published formulas." It does not, and
     for `P10_04/05/08/09` it *cannot* — the generator applies no explicit
     BLACK or RED bracket, so those are outside the declared class.
- **Scope document**: `docs/DEGREE10_REVERSE_SEARCH_SCOPE.md` separates the
  RECOVERY goal (met) from the EXHAUSTION goal (not met).


## C-REV-02 — the corrected degree-10 reverse run

- **Permitted wording (exact)**: "A three-element compact basis for Q10 was
  independently recovered in the pure `N^(4125)`^5 sector and shown to span the
  same quotient space as the graph-derived and published-formula bases."
- **Corrected accounting**, all terminal statuses recorded:

      planned 259 = evaluated 259 + structurally_rejected 0
                  + evaluation_error 0 + interrupted 0

  Reconciles exactly. **0 candidates outside the atlas span**, so every
  evaluated candidate is a genuine Lorentz scalar. 22 of 259 have nonzero
  quotient image; recovered rank **3 of 3**; 926 s; peak 2811 MB.
- **SUPERSEDES** the original pilot's statistics. That run silently skipped 9
  of 588 candidates with a bare `continue`, and its M-containing sectors were
  computed with a mixed-variance `M` block that was not a Lorentz scalar. Those
  numbers must not be quoted.
- **The rank-generating candidates CHANGED under the fix.** With a valid `M`
  block, the first two rank directions now come from `M+N1050+N4125^3` and
  `M+N4125^4`; the third is the `N4125^5` candidate. In the defective run all
  three came from `N4125^5`, because the M-containing candidates were not
  scalars and could not lie in the span.
- **The pure-`N4125` basis is unaffected** by that change: `N4125^5` contains
  no `M`, that sector is fully enumerated (15 of 15 canonical), and the triple
  was validated at two primes, on fresh samples, with determinants 8930 and
  29506, mutual inverses against both Level-A and the published basis, and
  removal minimality 3 → 2.

### Required claim distinctions, each machine-tested where possible

| distinction | status |
|---|---|
| recovery ≠ exhaustion | 5 of 21 sectors capped; sectors sampled at 15 |
| pilot coverage ≠ complete sector enumeration | 259 of ≥16 522 canonical |
| source independence ≠ clean-room independence | AST test proves no import; the author had read the paper |
| alternative basis ≠ canonical | no universal minimality theorem |
| minimal within a searched class ≠ universally minimal | `N4125^5` swept; the class was not |
| modular agreement ≠ characteristic-zero proof | rational reconstruction UNCERTIFIED at 2 primes |
| ambiguity removal ≠ ambiguity resolution | the pure-N basis avoids AMB-02; it does not settle the source reading |

- **AMB-02, measured**: changes the **coordinates** of `P10_11`, not the
  **span**. Both readings give rank 3, union rank 3. A span-level basis claim
  is unaffected; coordinate-level and change-of-basis claims are not.
