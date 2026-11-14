# Exact checkpointed degree-10 pipeline

## Result

The order-10 computation found twelve explicit connected primitive directions.
Together with the two products

```text
I4_1 * I6_1
I4_1 * I6_2
```

they have exact finite-field value rank 14. Adding the twelve connected
Jacobian rows to the nine generators through order 8 gives cumulative rank 21.
This holds at seeds 20260729, 20260730, and 20260731 under both primes 32749
and 32719. The full formulas, canonical graph IDs, discovery locations, values,
costs, and validation records are in `results/10d_order10.json`.

The statement “all twelve” uses two facts with different scopes:

1. the twelve saved graph contractions are independently verified exact
   directions; and
2. the published Hilbert series supplies the upper bound of twelve new
   primitive directions at degree 10.

The graph catalog itself was also generated exhaustively, but discovery stopped
after attaining the Hilbert target. No floating-point SVD is used.

## Exact catalog

A connected contraction of ten 5-forms is a connected loop-free 5-regular
multigraph on ten vertices. Self-duality removes multiplicity five, so edge
multiplicity is at most four. Its simple support consequently has degree from
two through five.

The exact generator is:

```text
geng -cq -d2 -D5 10 RESIDUE/64 | multig -q -T -m4 -r5
```

The 64 canonical support-graph shards contain 187,392 weighted isomorphism
classes. Exact incidence-graph certificates find 187,392 distinct canonical
IDs. The logical catalog fingerprint is:

```text
17b112b6aa4e52ad9a2786802397d03675d9fd5910d59f21208aaa540b01d9e7
```

Each gzip JSONL shard is written atomically with a deterministic logical
SHA-256 manifest. Interrupted generation verifies and skips completed shards.
The generated catalog is 11,214,833 compressed bytes and is deliberately not
committed; it regenerates from the commands below.

## Reproduction

Create the environment and pass the nauty 2.9.3 executables explicitly if they
are not on `PATH`:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements-lock.txt

.venv/bin/python scripts/degree10_pipeline.py generate \
  --catalog-dir work/degree10-catalog \
  --geng /path/to/nauty2_9_3/geng \
  --multig /path/to/nauty2_9_3/multig \
  --shards 64 --resume

.venv/bin/python scripts/degree10_pipeline.py schedule \
  --catalog-dir work/degree10-catalog \
  --schedule work/degree10-schedule.sqlite3

.venv/bin/python scripts/degree10_pipeline.py discover \
  --catalog-dir work/degree10-catalog \
  --schedule work/degree10-schedule.sqlite3 \
  --checkpoint work/degree10-final.checkpoint.json \
  --lower-result results/10d_order8.json \
  --prime 32749 --seed 20260729 \
  --schedule-mode hash-bounded --max-greedy-peak 100000000 \
  --checkpoint-every 5

.venv/bin/python scripts/degree10_pipeline.py validate \
  --catalog-dir work/degree10-catalog \
  --schedule work/degree10-schedule.sqlite3 \
  --discovery-checkpoint work/degree10-final.checkpoint.json \
  --lower-result results/10d_order8.json \
  --primes 32749 32719 \
  --jacobian-seeds 20260729 20260730 20260731 \
  --value-seed-start 20260801 --value-samples 16 \
  --out results/10d_order10.json
```

Discovery is deterministic. With the saved settings it reaches rank 21 after
132 candidates and selects the same twelve canonical graph IDs. Every
checkpoint contains the run identity, engine fingerprint, cursor, fully
reduced modular basis, and selected generators under a payload checksum.

A clean checkout can revalidate the explicit committed formulas without the
generated catalog:

```bash
.venv/bin/python scripts/degree10_pipeline.py validate \
  --selection-result results/10d_order10.json \
  --skip-catalog-check \
  --primes 32749 32719 \
  --jacobian-seeds 20260729 20260730 20260731 \
  --value-seed-start 20260801 --value-samples 16 \
  --out /tmp/10d_order10_revalidated.json
```

For a cluster, run `generate --residue "$SLURM_ARRAY_TASK_ID" --shards 64` as
an array job into one shared catalog directory, then run `generate --resume`
once on the coordinator to verify every shard and create the global index.

The pinned Python dependencies have explicit licenses: NumPy 2.5.1 is
BSD-3-Clause, pynauty 2.8.8.1 is GPL-3.0-or-later, and pytest 9.1.1 is MIT.
The externally built nauty 2.9.3 gtools are Apache-2.0. No nauty binary or
mentor code is vendored in this repository.

## Profiling and optimizations

Measurements are from an arm64 Apple M1 host with 8 logical cores, 8 GiB RAM,
macOS 26.5.1, Python 3.13.13, and NumPy 2.5.1. The benchmark command detects
the current host rather than embedding these values. Raw distributions are in
`results/degree10_benchmarks.json`.

The exhaustive enriched catalog took 69.87 seconds across its shard manifests;
global checksum, canonical-ID, and cross-shard uniqueness verification took
34.18 seconds in the final run. The compressed catalog is 11.2 MB.

| Stage | Median | p90 | Sample |
|---|---:|---:|---:|
| Exact canonicalization | 0.079 ms | 0.146 ms | 500 |
| Greedy scheduling estimate | 0.068 ms | 0.124 ms | 500 |
| Globally optimal plan, final | 12.5 ms | 15.8 ms | 100 |
| Scalar value | 65.8 ms | 93.1 ms | 12 selected graphs |
| Optimized value + reverse Jacobian | 198 ms | 303 ms | 12 |
| Reference value + amputated Jacobian | 635 ms | 1,050 ms | 12 |
| Modular rank insertion | 0.106 ms | 1.87 ms | 12 |
| Atomic checkpoint | 1.36 ms | 6.34 ms | 20 |

Major before/after measurements:

| Change | Before | After | Effect |
|---|---:|---:|---:|
| Bitset subset DP + bounded plan cache, plan p50 | 95.2 ms | 12.5 ms | 7.6x faster |
| Same change, value/Jacobian p50 | 330 ms | 198 ms | 1.7x faster |
| Full two-prime validation | 116.3 s | 50.5 s | 2.3x faster |
| Cost-first discovery | 4/12 after 1,000 candidates, 225.3 s | — | correlated narrow graphs |
| Hash-stratified safe-cost discovery, same pre-plan engine | — | 12/12 after 132 candidates, 51.5 s | attains target |
| Hash-stratified discovery, final engine | — | 12/12 after 132 candidates, 27.3 s | final reproduction |

The compact derivative projection removes the dense 252 by 100,000 direction
matrix. Its direction matrix is 252 by 126; exact signed antisymmetric orbit
sums pull the dense gradient back to the self-dual tangent coordinates. The
optimized and dense projections are unit-tested against each other.

## Correctness boundaries

- Arithmetic and ranks are exact over finite fields.
- Both primes and all saved sample seeds independently give the target ranks.
- All twelve selected graphs give nonzero values in every Jacobian validation
  sample.
- All twelve selected graphs match the separate amputated backend at the first
  prime and seed.
- Exact catalog enumeration is claimed only for the defined connected
  contraction-graph class.
- Completeness of the twelve primitive directions uses the published
  Hilbert-series upper bound.
- No claim is made that the 81 functionally independent invariants constitute
  a complete polynomial generating set.

The future mentor spinor comparison interface is documented in
`docs/spinor_adapter.md`; no mentor source has been reconstructed or copied.
