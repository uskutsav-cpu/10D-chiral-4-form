# Degree-12 reverse pilot — plan only, not launched

**Status: PREPARED, NOT LAUNCHED.**
`scripts/reverse_engineer_degree12_pilot.py` exits with code 2 unless invoked
with `--i-mean-it`, so it cannot start by accident.

## 1. Objective and starting point

Independently recover Q12 rank 4 from quadratic-block topology, without
importing or encoding the published equation-(4.25) structures.

The starting point is a **negative** result, which is what makes this safe to
state openly: `P12_01`, `P12_02` and `P12_03` have combined Q12 rank **0 of 4**
at primes 32749 and 32717. All four compact Q12 directions are unknown, so
there is no published span the search could accidentally copy.

## 2. Sector survey

28 block multisets of six quadratic blocks. Slot counts run 12 (`M`x6) to 36
(six 6-slot `N` blocks). Survey cap 4000 raw topologies per sector; **19 of 28
sectors hit it**, so every count below is a lower bound.

Canonical lower bound across all sectors: **13 499**.

Cheapest sectors first: `M`x6 (1), `M`x5+`N4125` (2), `M`x5+`N1050` (6),
`M`x4+`N4125`x2 (11), `M`x4+`N1050`+`N4125` (42), `M`x4+`N1050`x2 (106).

## 3. First shard, and why

    sector          N4125 x6      (36 slots, canonical >= 107)
    shard           residue 0 of modulus 4
    fit prime       32749
    holdout prime   32717

**Rationale.** The degree-10 recovery came *entirely* from the analogous
`N4125`x5 sector — all three recovered directions, and it was also the cheapest
N-only sector to enumerate. That is an empirical prior derived from this
project's own reverse search, **not** from the published formulas, so using it
does not compromise independence.

This is a prior about where to look first, not a restriction. If the shard
yields nothing, the cheapest-first sector order in
`degree12_reverse_pilot_plan.json` is the fallback sweep.

## 4. Hard limits

| control | value |
|---|---|
| heavy workers | 1 |
| `--max-rss-mb` | 1500 |
| `--stop-after-seconds` | 5400 |
| `--stop-after-candidates` | 400 |
| `--topology-cap` | 4000 |
| checkpoints | immutable, local temp, **never iCloud** |
| generation | streamed |

## 5. Cost must be measured, not extrapolated

`--measure-only` times a handful of contractions and exits. Use it first.

The degree-10 figure of ~150 ms per contraction is **not** transferable:
degree-12 candidates contract six 6-index operands rather than five, and the
intermediate sizes differ. Any duration quoted before running `--measure-only`
is a guess.

This repository has already recorded one wrong runtime forecast built from
three points on a rising prefix of a non-monotone cost curve
(`P12_RESOURCE_FAILURE_ANALYSIS.md`). The rule adopted then applies here: do
not extrapolate a completion time from fewer than ~5 points spread across the
full work range, and prefer "unknown, monitoring" to a number that will be read
as a forecast.

## 6. Known hazards carried over from degree 10

Two defects were found in the degree-10 work, both of the same kind, and both
would recur here unchecked:

1. **Metric placement.** `P10_07` shipped raising both ends of three edges;
   the reverse engine shipped handing `M_{a}{}^{b}` to a routine expecting
   all-lower operands. Both were invisible to homogeneity and to rotations and
   both were caught by a **boost** test. Any degree-12 evaluator must be boost
   tested before its numbers are believed.
2. **Sampling bias.** The degree-10 published-candidate result initially read
   "everything projects to zero" because only the five *simplest* formulas were
   implemented, and simplicity is exactly what determines whether a candidate
   reaches the quotient. A null result over a subset chosen for convenience is
   not evidence about the whole set.

## 7. Prerequisite

Do not launch until the degree-10 benchmark is frozen: reverse rank 3
recovered, span equality with the published Level-B basis proven on both fit
and holdout primes, and the clean-clone reproduction passing.
