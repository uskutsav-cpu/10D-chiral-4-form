# Live state discovery

Discovered by direct inspection on 2026-08-01, before any new work. Machine-readable
form in `results/live_state.json`.

## Repositories

| name | path | role | branch | HEAD | dirty | tags |
|---|---|---|---|---|---:|---|
| trace+bridge | `work/selfdual-5form-invariants` | **user-authored**: tensor implementation, bridge package, manuscripts, certificates | `research/maximal-chiral-four-form-program` | `ba19944` | 0 | `q10-freeze-v1`, `order10-verified-30fd0f6` |
| mentor archive work copy | `work/spinor-work` | **THIRD PARTY** — not redistributable | `main` | `2a53370` | 2 | none |
| clean clone | `work/clean-clone` | disposable reproduction clone | same branch | `8f91a16` | 16 | inherited |

The trace and bridge implementations live in one repository but are separate
packages (`src/sdinv`, `spinor_trace_bridge/src/sdbridge`) with no shared code.

## Reconciliation against the specification's stated starting point

| specification says | actually found |
|---|---|
| trace frozen at `5c6a883`, tag `q10-freeze-v1` | confirmed; tag resolves to `5c6a883`, which is an ancestor of HEAD |
| 198 passing tests | **199** — a product-decomposition regression test was added after the freeze |
| mentor baseline frozen at `b70bf33` | **that commit does not exist on this machine.** Every repository was searched with `git cat-file`. Only the raw archive survived; it was re-imported as `2a53370` and the baseline re-verified. No result is attributed to `b70bf33`. |
| 83 selected candidates = 70 port graphs + 13 tensor words, degrees 1/2/6/12/62 | confirmed by direct inspection of `selected_graphs.json` |
| float64 rank-81 gap rested on one point and one seed | confirmed, and superseded: see below |
| low-degree no-Hilbert-stop reached 1, 2, 7 through degree 8 | confirmed and extended — degree 8 terminated on `candidate_exhaustion` at rank 7 |
| bridge work may already have corrected the real form | confirmed: the frame is split `(5,5)`, not Euclidean. Verified from anticommutators, not from directory names. |

## Work already complete, verified rather than repeated

- Exact modular Clifford algebra, frame congruence, forward map and left inverse,
  with kernel and image identified by **span equality**; 72 tests at two primes.
- Common-sample comparison at degrees 4, 6, 8, 10.
- Exact analytic Jacobian by amputation for port graphs, validated by Euler's identity.
- Both archived Jacobians re-analysed under an explicit noise-floor rule.

## The gap this execution must close

The 13 structured tensor-word candidates were **never implemented in exact
arithmetic**. Three required results depend on them and cannot be reached without
them:

1. all 83 candidates having a terminal record;
2. exact rank 81 — the port-graph subset alone reaches 59, and 70 port graphs
   cannot exceed 70 in any case;
3. degree-8 spinor rank 7 — the port-graph family reaches 6.

That is the critical path and it is where this execution starts.

## Running processes

The degree-10 no-Hilbert-stop scan was killed (SIGKILL) after roughly four hours
of CPU time at rank 12 of 14. No healthy computation was terminated by this
execution.
