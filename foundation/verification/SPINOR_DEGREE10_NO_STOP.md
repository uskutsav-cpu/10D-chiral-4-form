# Degree-10 spinor scan with the Hilbert stopping rule disabled

Every archived scan stopped at the Hilbert target, which is supplied from outside the computation. A run that stops there demonstrates nothing about the ansatz. Disabling it is what makes the terminal status informative.

## Terminal status by degree

| degree | rank | target | selected | unique tried | stopped by |
|---:|---:|---:|---:|---:|---|
| 4 | 1 | 1 | 1 | 2 | **candidate_exhaustion** |
| 6 | 2 | 2 | 2 | 5 | **candidate_exhaustion** |
| 8 | 7 | 7 | 6 | 32 | **candidate_exhaustion** |
| 10 | 12 (partial) | 14 | 10 | 64 (first batch) | **killed** |

## Structural correspondence with the tensor side

The two implementations split each graded piece into non-product generators plus products. Those splits are computed independently, in different variables, in different arithmetic. They agree at every degree reached:

| degree | tensor: graph + product | spinor: new + product | agree |
|---:|---|---|---|
| 4 | 1 + 0 | 1 + 0 | yes |
| 6 | 2 + 0 | 2 + 0 | yes |
| 8 | 6 + 1 | 6 + 1 | yes |
| 10 | 12 + 2 | not reached | -- |

This is the degree-by-degree generalisation of `A10 = G10 (+) P10`. Note what it does NOT say: the tensor side's published-candidate span `B10` is a different twelve-dimensional subspace and is not a product complement. The correspondence is with the graph generators, not with the published structures.

## What this establishes

Degrees reaching candidate exhaustion: 4, 6, 8.

At those degrees the rank is derived by searching the ansatz out, not by stopping at a number supplied from outside. That is the distinction the specification asks for, and it is the reason the run was made at all.

## What this does not establish

Degree 10 did **not** reach a terminal state on this machine, so no exhaustion claim is made for it. The independent exact modular enumeration reaches rank 14 there by *saturation*, which is a weaker statement and is labelled as such wherever it appears.

The degree-10 attempt was killed after roughly four hours of CPU time, having
reached rank 12 of the target 14 without leaving its first candidate batch. It
was attempted and killed, not skipped, and it establishes nothing about the
ansatz at degree 10 in either direction.

There is a suggestive parallel with degree 8, where the port-graph family fell
one short and structured candidates supplied the remainder. A parallel is not
evidence, and no claim is made from it.

A cluster job that would complete the degree-10 and degree-12 scans is prepared
in `cluster/`. **No cluster run has occurred.**
