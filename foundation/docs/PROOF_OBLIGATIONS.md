# Proof obligations

Open mathematical debts. Each has an ID, the claim it supports, what would
discharge it, and its current state. Nothing in `docs/CLAIM_LEDGER.md` may be
strengthened while its supporting obligation is open.

Baseline: `3ed32805b38ce34216b34888f6539e3538e90fb9`

Final classification for the JHEP submission, with one status per
obligation: `audit/PROOF_OBLIGATIONS_FINAL.md`. This file remains the
working ledger; that one is the decision record.

| ID | supports | state |
|---|---|---|
| PO-01 | C-ATLAS-05 | **DISCHARGED** (Phase 0) |
| PO-02 | C-GEN-01 | **CERTIFIED** (rule checked independently; derivation is mentor item G-10) |
| PO-03 | C-ATLAS-04 | OPEN |
| PO-04 | C-SEXTIC-* | OPEN (external) |
| PO-05 | C-TRM6-01 | OPEN |
| PO-06 | C-FLOW-01 | **DISCHARGED** (Phase 0) |
| PO-07 | C-FLOW-03 | OPEN |
| PO-08 | C-MIN-01 | **PARTIAL** — cardinality half PROVED; removal / general GL open |
| PO-09 | all `MOD-CERT` | **CERTIFIED** for degree ten except `B10 cap P10` |
| PO-10 | Phase 5 | NOT STARTED |

---

## PO-01 — canonicalisation is exact at every degree actually used — **DISCHARGED**

**Supports** C-ATLAS-01, C-ATLAS-05.

**Original concern.** An earlier Weisfeiler–Leman fallback was measured to
collide on precisely the graphs in play: at order 6, 49 true isomorphism
classes collapse to 39 WL keys. Because the enumerator keeps the first graph
per key, a collision **drops a genuine candidate**, and a dropped candidate can
only lower a rank — a confident undercount.

**Discharged, Phase 0.** Verified at `3ed3280`:

- `_canonical_wl` **no longer exists** in `src/sdinv/graphs.py`
  (`hasattr(sdinv.graphs, "_canonical_wl") == False`);
- `canonical()` dispatches to `_canonical_nauty` whenever pynauty is present,
  falls back to the exact n! form only for n ≤ 6, and otherwise **raises**
  `RuntimeError` rather than degrading silently;
- pynauty is importable at runtime (`pynauty is not None == True`) and pinned
  at `2.8.8.1` in `requirements-lock.txt`;
- order 8 is generated directly by nauty's `geng | multig` toolchain, so
  labelled duplicates never arise.

The failure mode is now loud rather than silent, which is the property that
matters. No residual undercount risk.

## PO-02 — `leading_field_degree` is correct, not merely self-consistent

**Supports** C-GEN-01, and through it the exhaustiveness of the degree-6
argument behind C-FLOW-02/03.

Phase 0 verified: 18 expected multisets = 18 present, none missing, none extra,
and leading degrees are additive across products for every generator. That
establishes *internal* completeness given the per-generator
`leading_field_degree` values.

**Still open**: those values are computed by the same code they validate. A
derivation — Tr(τ) has leading degree 4 because the free stress tensor is
traceless, Tr(τ^k) has leading degree 2k for k ≥ 2 — should be written out and
independently checked, ideally clean-room.

## PO-03 — are I12_61 and I12_62 polynomial syzygies?

**Supports** C-ATLAS-04.

Currently only a Jacobian dependence is established. A Jacobian dependence at a
generic point means functional dependence; it does **not** exhibit a polynomial
identity, and the manuscript may not say "syzygy" until one is constructed.

**Discharged by**: bounding the relation degree, enumerating candidate
monomials, computing the exact evaluation matrix, taking a modular nullspace,
reconstructing rationals, and verifying on fresh samples at independent primes
— or by proving no such identity exists in the searched range, with the range
stated.

**Target**: Phase 3.

## PO-04 — the K6 ↔ Σ₂ identification (EXTERNAL)

**Supports** C-SEXTIC-01/02/03 wording.

arXiv:2509.14350v2 proves (Σ₁, Σ₂) is a sextic basis in the spinor formalism
but does **not** publish the change of basis to (Tr(M³), K6). Our
identification is inferred.

**Discharged by**: obtaining the spinor implementation or the explicit map from
the authors, or deriving it independently and verifying numerically. Marked
external because it may need a person, not a computation.

## PO-05 — Tr(M⁶): analytic reduction or basis-independent formulation

**Supports** C-TRM6-01.

29 of 72 columns exceed the CRT uniqueness bound at 15 primes (modulus ≈
5.2e67). Informative contrast: the *flow* coefficients lifted cleanly at five
primes, so this is a height problem specific to Tr(M⁶), not a shortage of
effort. More primes is the wrong instrument.

**Discharged by** either
(A) an analytic identity via Cayley–Hamilton / trace identities / self-duality
    / Schouten identities, verified exactly; or
(B) a proof that every downstream theorem can be stated with Tr(M⁶) as a
    primitive stress-adapted element, making graph coordinates unnecessary.

**Target**: Phase 2. (B) is the more likely route and is sufficient.

## PO-06 — second independent holdout prime — **DISCHARGED**

**Supports** C-FLOW-01.

The committed artifact used five fit primes and **one** holdout (32717); the
standard requires at least two.

**Discharged, Phase 0.** Two assemblies were run over the six available
interacting certificates, each excluding a different prime:

| fit primes | holdout | result |
|---|---|---|
| 32749, 32719, 32693, 32771, 32717 | **32713** | passed |
| 32749, 32719, 32693, 32771, 32713 | **32717** | passed |

Both report `all_modular_and_rational_holdouts_passed: true`. The stronger
check also holds: the two fits agree on the new-forcing dimensions
`{4:1, 6:1, 8:3, 10:5, 12:21}` and produce **192/192 identical** reconstructed
coordinate vectors. Two different fit sets converging on byte-identical
rationals, each validated on the prime it excluded, is materially stronger
than a single holdout.

Reproduce: `scripts/assemble_interacting_stress_adapted.py` with
`--certificates` omitting one prime and `--validation-certificate` set to it.

## PO-07 — does the K6 statement survive field redefinitions and EOM?

**Supports** C-FLOW-03, and gates any physical reading.

The transport equation is an off-shell statement in fixed conventions. A
field-redefinition-dependent obstruction is not a physical obstruction.

**Discharged by**: classifying allowed local field redefinitions at each
degree, computing their action on q6, and determining whether q6 = 0 is
preserved; then repeating modulo the leading equations of motion.

**Until discharged**, no physical or Type IIB consequence may be drawn from
C-FLOW-03.

## PO-08 — minimality under basis change — **PARTIAL**

**Supports** C-MIN-01, C-MIN-02.

**Done (Phase 1): the permutation subgroup.** Random relabellings of the basis
at a degree, applied consistently to `basis`, `coordinates` **and**
`coefficient_monomial`, with the missing set tracked back through the inverse
relabelling. Four trials per degree, all invariant:

| degree | closure dim | missing count | recovered set |
|---|---|---|---|
| 10 | 11 | 3 | I10_6, I10_7, I10_12 |
| 12 | 68 | 4 | I12_59, I12_60, I12_61, I12_62 |

Certificate: `results/generalized_flow/minimality_certificates/basis_permutation.json`

**Why this is not the whole obligation, and why the obvious test is invalid.**
Each certificate target row is indexed by a `coefficient_monomial` that is
itself expressed in the *same* basis as the output `coordinates`. A general
invertible change of basis therefore re-expresses **both**. Multiplying the
coordinate rows by a random matrix `B` and re-running would *look* like a
basis-change test and would not be one: the resulting rows correspond to no
actual flow problem, and would very likely still return 3 and 4 while meaning
nothing.

**Still open: general GL.** Requires regenerating the certificates in the new
basis (~530 s/prime, plus changes to the generator machinery so the monomial
index is re-expressed consistently). Not attempted; scope recorded rather than
claimed.

### The cardinality half — **DISCHARGED**, analytically

The obligation "an argument that the intrinsic quotient dimension bounds the
minimal cardinality from below" is discharged by the following, which needs no
computation and no choice of basis.

**Proposition.** Let `A` be the degree-`d` atlas, `D ⊆ A` the reachable subspace,
and `π : A → Q = A/D` the quotient map. If a finite set `S ⊆ A` closes the
degree, i.e. `D + span(S) = A`, then

        |S| ≥ dim Q .

*Proof.* `D + span(S) = A` gives `π(span(S)) = π(A) = Q`, since `π` is surjective
and kills `D`. But `π` is linear, so `π(span(S)) = span(π(S))`, and a span of
`|S|` vectors has dimension at most `|S|`. Hence `dim Q ≤ |S|`. ∎

**Why this is basis-independent.** Every object in the statement is defined
without reference to a basis: `A` is a space of invariants, `D` is the set
reachable by the flow, `Q` is their quotient, and `dim Q` is an invariant of that
quotient. A change of basis of `A` induces an isomorphism of `Q` and leaves
`dim Q` unchanged, so the bound transports unchanged. In particular the observed
deficits — 3 at degree 10 and 4 at degree 12 — cannot be closed by fewer than 3
and 4 elements respectively **in any basis**, and the exhibited sets meet the
bound and are therefore of minimum cardinality.

Checked in `tests/test_quotient_cardinality_bound.py`, which verifies the bound
numerically against the recorded certificates and also verifies the hypothesis
the proposition needs (that the exhibited sets do close the degree). The test is
a guard on the *inputs*, not evidence for the proposition, which is proved above.

**What this does not discharge.** Minimality of *cardinality* is now
basis-independent. Minimality under *removal* — the statement that no single
element of the exhibited set can be dropped — is a stronger, set-specific
property, and it remains shown only in the fixed basis and under the permutation
subgroup. A general `GL` test still requires regenerating the certificates in the
new basis as described above, and is still not attempted. So PO-08 stays
**PARTIAL**, with its two halves now separated:

| half | status |
|---|---|
| cardinality bound `|S| ≥ dim Q`, any basis | **DISCHARGED** (analytic, above) |
| removal-minimality of the exhibited set under general `GL` | **OPEN** |

## PO-09 — exceptional primes

**Supports** every `MOD-CERT` row.

Rank mod p is a lower bound for the characteristic-zero rank, with equality
unless p divides a relevant minor. Agreement across 15 primes makes a
coincidence unlikely but does not prove it.

**Discharged by**: a Hadamard-type bound on the relevant minors giving a prime
count beyond which agreement forces equality, or one exact
characteristic-zero computation at a decisive rank.

### Which claims are actually exposed — **PARTIAL DISCHARGE**

PO-09 has been applied to every degree-10 headline number as a blanket caveat.
That is imprecise, and the imprecision cuts both ways: some of those numbers are
not exposed at all, and for the ones that are, the *direction* of a possible
failure is determined and worth stating. Only `rank_{F_p} <= rank_Q` is used
below.

**Not exposed. Unconditional over `Q`.** A space spanned by exactly `k` explicit
elements has `dim_Q <= k` for free; if its modular rank is also `k` then
`dim_Q >= k`, so `dim_Q = k` exactly and no prime can be bad:

| space | spanning set | modular rank | conclusion |
|---|---|---|---|
| `A10` | the 14 atlas elements | 14 | `dim_Q = 14` |
| `P10` | `I4_1*I6_1`, `I4_1*I6_2` | 2 | `dim_Q = 2` |
| `G10` | `I10_1 .. I10_12` | 12 | `dim_Q = 12` |
| `B10` | the 12 published candidates | 12 | `dim_Q = 12` |

The exact Jacobian bound `rank_Q >= 81` is likewise unexposed, for the separate
reason that the matrix is an integer reduction and the claim is one-sided.

**Was exposed; `D10` is now settled.** `D10` is built by admitting a generated
row only when it raises the rank *modulo p*, so the recorded 11 was a lower bound
over `Q` and the quotient was bounded above:

    dim_Q D10 >= 11        hence   dim_Q Q10 <= 3      <- SUPERSEDED
    dim_Q(B10 + P10) >= 13 hence   dim_Q(B10 cap P10) <= 1

The first line is superseded. Re-running the closure in exact rational
arithmetic, after lifting its non-integral targets by CRT and rational
reconstruction validated at a held-out prime, gives `dim_Q D10 = 11` exactly and
hence `dim_Q Q10 = 3`, with an explicit non-vanishing integer minor. See
`docs/D10_Q10_CHARACTERISTIC_ZERO_STATUS.md`. The second line stands.

Both consequences use the exact values of `dim_Q B10` and `dim_Q P10` from the
table above, which is why isolating the unexposed cases first was worth doing.

So a bad prime could only make the flow reach **more** than recorded and the
quotient **smaller** than 3, and could only make the published span meet the
products in **fewer** dimensions than 1 — that is, it would restore the `12 = 12`
coincidence that this project spent a session refuting. The refutation is
therefore the claim most worth hardening, not the atlas dimension.

**What is not yet done.** Neither exposed statement has been checked in
characteristic zero. Both would be discharged by exact evaluation over `Z` at the
same sample points — the samples are integer five-forms and the invariants are
integer contractions, so the obstruction is arithmetic width, not principle:
the pipeline is `int64`-modular throughout and would need arbitrary-precision or
CRT reconstruction with a height bound. Running further primes lowers the
probability but discharges nothing, and is not counted here as progress.

Recorded in the manuscript limitations rather than left in this file.

## PO-10 — the induction for any all-orders claim

**Supports** Phase 5.

No induction on degree currently exists. Degree-12 data cannot imply an
all-orders theorem. Required: base cases, a recursion step, preservation of
hypotheses, and an argument that no exceptional degree exists.

**Until discharged**, "all orders" is forbidden in every document (ledger F-5).

## PO-11 — bracket colour in equation (4.24)

**Supports** the Q10 Level-B basis.

`AMB-01` and `AMB-02` are unresolved because bracket **colour** does not
survive PDF text extraction, and colour is what fixes operation order. Both
readings of every ambiguous candidate are implemented and projected, and the
measurement is decisive rather than academic:

| candidate | Q10 image under the two readings |
|---|---|
| P10_10, P10_12 | identical |
| P10_09, P10_11 | **differ** |

Because `P10_10` is forced and only `P10_12` of the remainder is robust, **no
ambiguity-robust triple exists** among the twelve published candidates. The
selected basis attains the minimum of one source-reading-dependent member.

**Required**: a colour render of journal page 17 / arXiv page 25, or a
statement from the authors. This is a binding prerequisite for an
unconditional basis, not a tidiness item.

**Until discharged**: the basis is "preferred ambiguity-minimal", never
"ambiguity-robust". Only `P10_10` and `P10_12` are unconditional.

## PO-12 — exhaustion of the degree-10 block class

**Supports** any minimality theorem over compact block contractions.

The reverse benchmark met its RECOVERY goal — Q10 rank 3, independently, span
equal to the published Level-B span on fit and holdout primes. It did **not**
meet the exhaustion goal: 5 of 21 sectors are capped at 30 000 raw topologies,
and even the exhausted sectors were sampled at 40 candidates for evaluation.

**Required to discharge**: sweep every canonical candidate in the declared
class. Lower bound 16 522 candidates at 22 samples and ~150 ms is ≈ 15 h on one
worker, and the capped sectors add an unknown multiple. Streamed generation
with `--sector`/`--shard-*` makes this parallelisable.

**Until discharged**: forbidden to claim "complete enumeration of every M/N
contraction", or that the reverse search establishes minimality of anything.

## PO-13 — rational reconstruction of the Q10 change-of-basis matrices

**Supports** any characteristic-zero statement about the Level-A/Level-B map.

Both directions are certified modularly at two primes and verified mutually
inverse. Rational reconstruction was **attempted and is not certified**: two
primes give a CRT modulus of ~1.07e9, so a lift is unique only when numerator
and denominator are both below ~2.3e4, and the entries are generic residues of
that magnitude.

**Required**: more primes. The projection is checkpointed, so the four
remaining primes are incremental rather than a rerun.

**Until discharged**: the maps are modular certificates, not rational
identities.
