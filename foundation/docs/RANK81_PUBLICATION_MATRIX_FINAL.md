# Rank matrix — final

Artifacts: `results/rank81/full_rank_matrix_publication_final.{json,csv,sha256}`.
Producer: `scripts/aggregate_rank_matrix.py`. Cells: `results/rank81/cells/`.

## Result

    cells planned      15   (5 primes x 3 seeds)
    cells complete     15
    cells invalidated   0
    cells missing       0
    distinct ranks    [81]  -- every cell
    stable pivot rows      81   unstable 0
    stable pivot columns   81   unstable 0

Every cell: 83 candidates planned, 83 evaluated, 0 evaluation errors, 0
interrupted, 0 structurally rejected, 0 zero rows, Euler homogeneity 83/83,
schedule complete.

Primes `32693, 32713, 32717, 32719, 32749`; seeds `11, 22, 33`.

The pivot structure being *identical* across all fifteen cells is the stronger
statement. Equal ranks alone would allow each cell to be reaching 81 through a
different set of candidates; the same 81 pivot rows and the same 81 pivot columns
everywhere says the rank is carried by the same directions at every prime and
seed tested.

## Why 32707 is not in the matrix

An earlier attempt used `32707` as a holdout and lost three cells to

    RuntimeError: image has dimension 0, not 126

That is **not** an exceptional prime. The null-frame congruence was fixed only up
to a square-root branch, which reverses orientation, which flips the sign of the
Hodge star, which swaps the eigenspace the gamma map annihilates. See
`docs/FRAME_ORIENTATION_FINDING.md`. With the orientation pinned, 32707 behaves
like every other prime.

The matrix was completed with `32693` in its place and has not been recomputed,
because 15/15 at five primes already establishes what it is for. 32707 is now
available if a sixth prime is ever wanted.

## The aggregator refuses more than it accepts

`--self-test` mutates a copy of the cells eleven ways and requires each to be
rejected: missing cell, duplicate cell, mixed candidate ordering, mixed
coordinate dimension, mixed flop budget, 82 of 83 candidates, an evaluation
error, a zero row, a failed Euler check, a malformed provenance hash, and a cell
not marked complete. It then re-reads the cells to confirm none was modified, and
re-renders from a reversed input list to confirm the output is
order-independent. All eleven fire; the output is byte-identical on rerun.

A gate that has never fired is not a gate, which is why the rejections are
exercised rather than merely written.

## Provenance, including what is not established

The cells were produced by a rank-matrix run in a **separate working tree** of
this project (`~/Downloads/sdinv-jhep`) and record no source-critical hash.
Rather than trust them, the two cells that overlap this repository's own
certificate — `(32749, 11)` and `(32749, 22)` — were compared field by field:
total rank, row and column counts, per-degree block ranks, cumulative ranks,
pivot rows, pivot columns and pivot row candidate ids. **All agree.** All fifteen
share a single candidate-ordering hash.

What is **not** established: the remaining thirteen cells were not independently
recomputed here. The evidence is internal consistency across fifteen cells plus
exact agreement with this repository on the two that overlap.

## Standing

No manuscript claim depends on the matrix. The rank-81 statement rests on the
committed certificate and its explicit `81 x 81` minor, and the matrix is
corroboration across primes and seeds.
