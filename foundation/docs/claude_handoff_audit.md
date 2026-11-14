# Claude handoff audit

Audit of the Codex working state at takeover. Every claim below was checked
against repository evidence, not against the handoff narrative. Where the two
disagree, the repository wins and the disagreement is recorded.

- **Tree**: `~/Documents/Codex/2026-07-29/now/work/selfdual-5form-invariants`
- **Branch**: `stress-flow/classification-through-degree12` (already existed)
- **HEAD**: `123cde497a350b31cc7b0b4a49af8cb2af9ca060`
- **Preservation**: `.handoff/` — see `claude_takeover_inventory.md`. No
  destructive git command was run at any point.

## 1. Degree-12 certificate isolation — VERIFIED

The handoff reported that modifying the degree-12 engine had perturbed the
pinned certificate hash, and that Codex restored the engine byte-for-byte and
moved new machinery into separate modules.

Checked by hashing every file three ways (public commit / HEAD / working tree):

| file | status |
|---|---|
| `src/sdinv/__init__.py` | identical |
| `src/sdinv/catalog.py` | identical |
| `src/sdinv/checkpoint.py` | identical |
| `src/sdinv/contract.py` | identical |
| `src/sdinv/forms.py` | identical |
| `src/sdinv/graphs.py` | identical |
| `src/sdinv/modp.py` | identical |
| `src/sdinv/spinor_adapter.py` | identical |
| `results/10d_order12.json` | identical |
| `results/degree12_benchmarks.json` | identical |

`git diff --stat de696ca HEAD -- src/sdinv/` is **1702 insertions, 0
deletions**, entirely in five new modules: `exactmap.py`, `interaction.py`,
`invariant_registry.py`, `sextic.py`, `stress.py`.

The published semantic fingerprint
`26b61c440b64b1beefe5d19f8eba4fcb9f299502a58796e82900ed2080605c25`
is present in `results/degree12_benchmarks.json`, and
`test_semantic_fingerprint_ignores_only_runtime_measurements` passes.

**Conclusion: the claim holds.** The atlas is isolated from the new work.

## 2. Test suite

- Full suite: **70 passed**, 0 failed, 43.7 s.
- The handoff reported 61/61. The extra tests come from commits `b789ec6`
  and `123cde4` plus the uncommitted `test_degree12.py` addition. Not a
  discrepancy in substance.

Milestone tests re-run explicitly (24 passed):

| milestone | test |
|---|---|
| anti-self-dual projection | `test_registry_gradients_are_anti_selfdual_and_obey_euler` |
| eq (3.3)/(3.4) equivalence | `test_equations_3_3_and_3_4_are_exactly_equivalent` |
| eq (2.33) → V(I4) | `test_general_interacting_stress_reduces_to_v_i4_formula` |
| derivative normalization | `test_registry_gradient_has_paper_normalization` |
| intrinsic sextic | `tests/test_sextic.py` (4) |
| degree-12 atlas order | `tests/test_degree12.py` (4) |

## 3. Intrinsic sextic basis — VERIFIED against raw artifact

Source: `results/stress_flow/change_of_basis/sextic_intrinsic.json`.

The stored `intrinsic_to_registry` matrix is `[[32/3, -1/1125], [0, 3/125]]`,
read **column-wise** against `registry_basis = [I6_1, I6_2]`:

    Tr(M^3) = (32/3) I6_1
    K_1050  = -(1/1125) I6_1 + (3/125) I6_2

- determinant `32/125` — matches the stored value and recomputes correctly
  as `(32/3)(3/125) - (-1/1125)(0)`
- inverse `registry_to_intrinsic = [[3/32, 1/288], [0, 125/3]]`, and
  `(3/32)(125/3) = 125/32 = (32/125)^-1` ✓
- quotient coordinate `q(c1 I6_1 + c2 I6_2) = 125 c2 / 3` ✓
- validated over primes 32749 / 32719 / 32693, four seeds each, with the
  independent dense 1050-tensor contraction agreeing exactly
  (`direct_1050_value == K_1050` on every sample)

**All coefficients quoted in the handoff are correct.** The only nuance is
matrix orientation: the handoff presents the relations as formulas, the
artifact stores them column-wise. No numerical disagreement.

Intrinsic definitions recorded in the artifact:

    Tr(M^3) : M_mu^nu M_nu^rho M_rho^mu
    N1050   : Lambda^mn_[abc Lambda_de]fmn
    K_1050  : N1050_[abc,de]f N1050^[abc,]_[ghi] N1050^[def,gi]h

A caveat is already recorded in the artifact and is worth repeating: the
source proves `(Sigma_1, Sigma_2)` is a sextic basis in the spinor formalism
but **does not publish the change of basis** to `(Tr(M^3), K_1050)`. The
identification here is ours, not the paper's.

## 4. Uncommitted diff — understood and reproducible

`src/sdinv/invariant_registry.py` (+54 −12)

1. imports `graph_to_record`;
2. reorders `degree12_product_items()` so `I4_1^3` precedes the six
   `I4_1*I8_k` entries — this makes the registry's degree-12 ordering match
   the committed artifact order;
3. adds `load_verified_registry_through_degree12(repository_root)`, which
   loads the 62 committed generators, extends the registry, and **raises**
   if the resulting basis order differs from `results/10d_order12.json`.

`tests/test_degree12.py` (+14) adds
`test_committed_degree12_basis_loads_in_exact_artifact_order`, asserting the
72-element basis loads in exactly the artifact order.

The reordering is a genuine behavioural change to the *registry*, but it is
guarded by an equality assertion against the committed artifact and does not
touch the atlas engine or its outputs (section 1). It is the intended
"canonical order" milestone.

## 5. Tr(M^6) — PARTIALLY COMPLETE at takeover

`results/stress_flow/certificates/static_degree12_32749.json`:

    prime 32749 (one prime only)   basis_rank 72   stress_rank 4
    quotient_dimension 68          stacked_columns 504   seconds 154.3
    seeds [20260729, 20260730, 20260731, 20260732]
    target n12[0-1^4,0-11^1,1-2^1,2-3^4,...,10-11^4]

Four target rows, each a 72-vector over F_32749:

| target | first entries | note |
|---|---|---|
| `tr_M2^3` | `[8, 0, 0, 0, …]` | analytic product |
| `tr_M3^2` | `[0, 29224, 0, 0, …]` | analytic product |
| `tr_M4*tr_M2` | `[0, 0, 0, 0, 2, …]` | analytic product |
| `tr_M6` | `[23343, 25757, 31060, 0, 2083, 31128, …]` | full atlas reduction |

So the handoff's "remaining row" **was in fact computed**, for one prime. What
is missing relative to the stated exactness bar:

- only 1 of ≥3 required validation primes;
- no independent holdout prime;
- no rational reconstruction (residues only);
- no fresh-sample verification of a reconstructed rational vector.

**Resumption action taken:** running `compute` for 32719 and 32693 (completing
the three validation primes) and 32717 (holdout), using the project's
canonical `DEFAULT_PRIMES` and the same default seeds. The script already
supports `assemble --certificates … --validation-certificate …`, which is the
CRT + reconstruction + holdout step.

## 6. Recorded discrepancy: eq (2.36) trace sign

The handoff reports that eq (2.36) appears to print a trace sign opposite to
the displayed eq (2.33), with the conformal zero-trace condition unaffected.

**Status: not yet independently re-derived.** No formula has been altered to
force agreement. This is carried forward as an open item for
`docs/published_formula_reproduction.md`, where both literal forms must be
preserved and the convention that follows from the derivation identified.

## 7. Open items

1. Finish Tr(M^6): 3 primes + holdout, rational reconstruction, fresh-sample
   check. **In progress.**
2. Commit the uncommitted milestone (currently only protected by `.handoff/`).
3. Complete the static stress subalgebra table for degrees 4, 6, 8, 10, 12.
4. Everything dynamical — flow reduction, role of K6, closed families,
   minimal generalized flow — is **not started**.
5. Document the eq (2.36) sign question with evidence on both sides.

## 8. Audit verdict

The diff is understood and reproducible; the atlas is provably untouched; the
test suite is green at 70; the intrinsic sextic coefficients check out against
raw artifacts over three primes with an independent contraction. There is no
reason to reject or rewrite any of the existing work, and none was rewritten.
