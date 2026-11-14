# Exact bounded degree-12 pipeline

## Result and proof boundary

The committed certificate contains 72 linearly independent homogeneous
degree-12 scalar polynomials:

```text
I4_1^3
I6_1^2
I6_1*I6_2
I6_2^2
I4_1*I8_1, ..., I4_1*I8_6
I12_1, ..., I12_62
```

The 62 `I12_*` entries are explicit connected valence-5 multigraph
contractions in `results/10d_order12.json`. Together with the 21 primitive
generators through degree 10, the first 60 raise generic Jacobian rank to 81.
The last two, `I12_61` and `I12_62`, are independent degree-12 polynomials but
add no new cumulative functional direction.

There are two logically separate completeness statements:

1. Exact stacked-gradient rank 72 gives a lower bound of 72 on the homogeneous
   degree-12 scalar space. The supplied Hilbert-series dimension 72 is the
   upper bound, so the degree-12 list is complete.
2. Exact generic Jacobian rank 81 reaches the known dimension of the orbit
   quotient, so the selected functions are a complete functional coordinate
   set through degree 12.

Finite-field generic ranks are lower bounds in characteristic zero. Agreement
at four deterministic samples under each of three primes—32749, 32719, and
32693—guards against exceptional points and exceptional reductions. No
floating-point rank tolerance is used.

The supplied eight-page spinor calculation note was rendered and visually
audited as an independent target. Its SHA-256 is recorded in the result. No
mentor code or candidate list was copied.

## Why the degree-12 rank uses several points

At one point, every product gradient is a linear combination of lower-degree
gradients. Consequently the ten degree-12 product gradients cannot have rank
ten in one 126-column tangent space. This is not a polynomial relation.

For same-degree homogeneous polynomials, an identity

```text
c1 f1 + ... + cr fr = 0
```

would force the same linear combination of gradients to vanish at every
point. Concatenating each polynomial's gradients at several generic points
therefore provides an exact polynomial-independence test. Discovery uses two
points (252 columns); validation uses four (504 columns). In every validation
field all ten products and all 62 connected contractions are pivots.

## Exact graph shards

Order-12 candidates are produced directly without labelled duplicates:

```text
geng -cq -d2 -D5 12 RESIDUE/4096 | multig -q -T -m4 -r5
```

Every gzip JSONL record stores the unambiguous upper triangle, a canonical
incidence-graph SHA-256 ID, and scheduling metadata. Generation then rereads
the entire shard, checks its logical SHA-256, recomputes every canonical ID,
and proves within-shard uniqueness.

The discovery used these deterministic shards:

| Residue | Candidates | Logical SHA-256 |
|---:|---:|---|
| 63 | 23,480 | `163075c73ae68dc9e25788e80e3646d6921d1ea4603738c6b249df2c5848e53c` |
| 0 | 112 | `8eec0cb679495f2f3aea45fa05aa9a470bcd8bad02983e2c35a7dcb728509757` |
| 1 | 53 | `f94496d896dc9427342c67c382070cbdb82d02b2ceefbf083d29c24fc8a09036` |
| 2 | 375 | `f3b0a9e0ff192099ab23d62a2e6e343d6f5c823171204285da0e1530ef3a4ef5` |
| 3 | 217 | `cde856426cc133e8fc47be3801f6d638b67ba4d02a0c04dde60b6c9f11cd3061` |
| 4 | 28 | `e4e43ede84fdbd3f68153aad51c887cc30671625fd88d05ab28c65c2ca1141a3` |

These are exact non-isomorphic search shards. The result makes no exhaustive
global order-12 graph-count claim.

## Planning, memory, and checkpoints

Every evaluated graph first receives a globally optimal binary contraction
tree from subset dynamic programming. Candidates are rejected before tensor
allocation unless:

- maximum intermediate boundary rank is at most 7;
- maximum pair-union rank is at most 9;
- the conservative retained-forward/reverse plus workspace estimate is at
  most 2 GiB.

The process also stops if observed peak RSS crosses 3 GiB. Discovery peaked at
1.03 GB and clean-artifact validation at 1.18 GB. The exact plan profile, cost, width, and
memory estimate are stored with every selected graph.

The rank checkpoint contains a content-addressed run identity, schedule
cursor, two fully reduced exact rank sieves, selected records, timings, and a
payload checksum. It is atomically replaced after every accepted direction,
at configured candidate intervals, and at least once per minute. Candidate
updates are staged and committed as one state assignment. Schedule extensions
are accepted only when the old candidate list is a byte-identical,
hash-verified prefix; this path was exercised during discovery.

## Reproduce discovery

Create the environment and supply nauty 2.9.3 `geng` and `multig` paths. The
following records the staged search used for the committed basis:

```bash
.venv/bin/python scripts/degree12_pipeline.py generate \
  --geng /path/to/geng --multig /path/to/multig \
  --shards 4096 --residues 63 --resume

.venv/bin/python scripts/degree12_pipeline.py schedule \
  --shards 4096 --residues 63 \
  --candidate-limit 5000 --target-safe 160 --minimum-safe 80

# This three-candidate stop verifies interruption/resume.
.venv/bin/python scripts/degree12_pipeline.py discover --max-evaluated 3
.venv/bin/python scripts/degree12_pipeline.py discover

.venv/bin/python scripts/degree12_pipeline.py schedule \
  --shards 4096 --residues 63 \
  --candidate-limit 2500 --target-safe 320 --minimum-safe 200 \
  --extend-existing

.venv/bin/python scripts/degree12_pipeline.py generate \
  --geng /path/to/geng --multig /path/to/multig \
  --shards 4096 --residues 0 --resume
.venv/bin/python scripts/degree12_pipeline.py schedule \
  --shards 4096 --residues 63 0 \
  --candidate-limit 1500 --target-safe 400 --minimum-safe 320 \
  --extend-existing --extension-new-sources-only

.venv/bin/python scripts/degree12_pipeline.py generate \
  --geng /path/to/geng --multig /path/to/multig \
  --shards 4096 --residues 1 2 3 4 --resume
.venv/bin/python scripts/degree12_pipeline.py schedule \
  --shards 4096 --residues 63 0 1 2 3 4 \
  --candidate-limit 2200 --target-safe 600 --minimum-safe 400 \
  --extend-existing --extension-new-sources-only
.venv/bin/python scripts/degree12_pipeline.py discover
```

The final run examined 376 scheduled candidates: 62 were selected and 314
were exact polynomial dependencies. There were no memory rejections and no
early functional dependencies.

## Revalidate a clean checkout

The explicit formulas do not require the ignored `work/` shards or checkpoint:

```bash
.venv/bin/python scripts/degree12_pipeline.py validate \
  --selection-result results/10d_order12.json \
  --primes 32749 32719 32693 \
  --seeds 20260729 20260730 20260731 20260732 \
  --out /tmp/10d_order12_revalidated.json
```

This recomputes all 996 value/Jacobian contractions and 186 finite-field boost
checks. Compare the semantic result fingerprint in
`results/degree12_benchmarks.json`; timing and RSS fields are deliberately
excluded from that fingerprint.

## Benchmark summary

The recorded machine was an Apple M1 with 8 GiB RAM, Python 3.13.13, and NumPy
2.5.1. Full distributions, including minima and worst cases, are in
`results/degree12_benchmarks.json`.

| Stage | Median | p90 | p99 | Worst |
|---|---:|---:|---:|---:|
| Exact canonicalization | 0.154 ms | 0.174 ms | 0.259 ms | 0.306 ms |
| Exact global planner | 192 ms | 298 ms | 553 ms | 766 ms |
| Discovery value + Jacobian | 0.978 s | 1.84 s | 3.12 s | 5.10 s |
| Validation value + Jacobian | 0.519 s | 1.21 s | 2.48 s | 8.15 s |
| Rank insertion | 0.199 ms | 0.647 ms | 4.57 ms | 48.5 ms |
| Atomic checkpoint | 21.7 ms | 142 ms | 788 ms | 4.85 s |
| Lorentz scalar value | 0.206 s | 0.402 s | 0.763 s | 1.63 s |

The schedule contains 214 width-6 and 386 width-7 plans. Estimated peak bytes
have median 618 MB, p90 794 MB, p99 1.21 GB, and worst case 1.36 GB.
