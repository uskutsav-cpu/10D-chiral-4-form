# Projection checkpoint validity semantics

## 1. The rule

> **Block on what changes the computed values. Record everything else as
> provenance.**

Two failure modes bracket the design, and both were reached in practice:

- **too strict** — `source_commit` was a blocking key, so committing a
  documentation change refused to resume and would have discarded a completed
  two-prime projection. That is the opposite of what a checkpoint is for.
- **too lax** — `evaluator_version` was a hand-maintained string. A string only
  invalidates when someone remembers to bump it, and the forgotten case serves
  values computed by the *old* code that look entirely plausible.

## 2. Two levels of key

### Global identity — blocks the whole store

| key | blocks? | why |
|---|---|---|
| `atlas_sha256` | **yes** | coordinates are reported against a column order |
| `quotient_sha256` | **yes** | a different D10 closure gives different quotient vectors from identical atlas coordinates |
| `evaluator_version` | **yes** | explicit manual override, retained as an escape hatch |
| `modular_backend` | **yes** | a contraction-path change can alter results |
| `prime`, `seed_base`, `degree` | **yes** | different inputs entirely |
| `source_commit` | **no** | provenance; recorded, and drift surfaced as `store.commit_drift` |
| `block_fingerprint` | **no** | enforced per-unit instead — see below |

The atlas hash is deliberately **order-sensitive**: `json.dumps(names)` without
`sort_keys`. Two runs with the same column *set* in a different order are not
interchangeable, because every stored coordinate vector is indexed by position.
Sorting the hash input would silently accept a permuted basis.

### Per-unit fingerprint — blocks exactly the affected units

`load_unit(..., fingerprint=...)` rejects a unit whose stored fingerprint
differs. This is what makes invalidation surgical:

- **atlas units** are keyed by `block_fingerprint()` — the hash of the shared
  tensor machinery (`stress.py`, `modp.py`, `index_symmetry_ops.py`);
- **formula units** are keyed by `evaluator_fingerprint(evaluator)` — that
  evaluator's own source *plus* the shared blocks.

Reimplementing one formula therefore invalidates that formula's 22 samples and
leaves the expensive atlas cache intact. Putting `block_fingerprint` in the
global identity instead would detonate the whole store for the same change.

## 3. Why the evaluator fingerprint includes the shared modules

Almost every equation-(4.24) evaluator is a thin wrapper over
`composite_n1050`, `_raise_axes`, `mod_einsum` and the bracket engine. Hashing
only the evaluator's own text would miss a change to `composite_n1050`, which
changes every cached value while leaving each evaluator's source untouched.

Including the shared modules over-invalidates on comment edits. That is the
right way to be wrong: a needless recomputation costs minutes, a stale value
costs a wrong scientific claim.

## 4. Measured effect

| run | wall time | peak RSS |
|---|---:|---:|
| cold, two primes, 16 evaluators | 912 s | 531 MB |
| resumed after adding one evaluator | **55 s** | 267 MB |

The resumed run recomputed only the new evaluator's units.

## 5. Test cover

`tests/test_projection_checkpoint_semantics.py` (19) and
`tests/test_projection_checkpoint_resume.py` (9):

1. documentation-only commit preserves the atlas cache;
2. `source_commit` still recorded, drift surfaced;
3. formula fingerprint change invalidates that formula;
4. formula change does **not** invalidate unrelated atlas units;
5. each semantic identity key refuses resume (7 parametrised cases);
6. atlas hashing is order-sensitive;
7. corrupted checksum rejected;
8. partially written record ignored;
9. missing unit reads as absent, not an error;
10. duplicate save is idempotent and not double-counted;
11. no two evaluators share a fingerprint;
12. the shared machinery enters the fingerprint.

## 6. Not covered

- **Concurrent writers.** The design assumes one heavy worker at a time, which
  is also the machine's memory constraint. Two processes writing the same unit
  path would race; `os.replace` keeps each file internally consistent, but no
  locking exists.
- **Reordered assembly.** Units are read back by explicit `(prime, sample,
  column)` key rather than directory order, so assembly order cannot scramble
  a matrix. This is a property of the caller, not something the store enforces.
- **Worker failure mid-flush.** The manifest is atomic, but a crash between
  `save_unit` and `flush_manifest` leaves a unit on disk that the manifest does
  not list. It is still found by `load_unit`, which checks the filesystem, so
  the effect is a stale count in `completed()`, not lost or wrong data.
