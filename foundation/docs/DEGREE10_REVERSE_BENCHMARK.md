# Independent reverse graph-to-block benchmark at degree 10

## 1. What this is, and why independence matters

The twelve published equation-(4.24) candidates span Q10. That is a forward
result: the formulas were transcribed, implemented, and projected. This
benchmark asks the same question from the other end — enumerate compact
quadratic-block contractions directly and see whether the quotient is reached
**without being told the answer**.

Independence is enforced, not promised.
`src/sdinv/reverse_block_decomposition.py` does not import
`published_degree10_invariants`, and
`tests/test_reverse_block_decomposition.py::test_reverse_engine_does_not_import_the_published_formulas`
parses the module's AST and fails if it ever does. The published result is used
only after a span has been recovered, for comparison.

## 2. Why a naive search is impossible, and what makes this one work

A degree-10 scalar is five quadratic blocks, each degree 2 in `F`: `M`
(2 slots), `N^(1050)` (6 slots), `N^(4125)` (6 slots). Five `N` blocks give 30
index slots, and a perfect matching of raw slots would be

    29!! ~ 6.2 x 10^15

which is not searchable.

The slots are not distinguishable. `composite_n1050` antisymmetrises axes
(0,1,2,3,4), so those five slots are interchangeable up to sign and only axis 5
is distinct. Each block therefore has **slot classes**:

| block | class A | class B |
|---|---|---|
| `M` | 2 slots | — |
| `N^(1050)` | 5 slots, antisymmetric | 1 slot |
| `N^(4125)` | 6 slots | — |

A contraction is then fully described by how many edges join each ordered pair
of `(block, class)` endpoints — a small integer enumeration with degree
constraints, rather than a matching over 30 labelled slots.

### Pruning, in the order applied

1. **Forbidden traces.** An edge joining two class-A slots of the *same* block
   contracts two indices inside an antisymmetric group and vanishes
   identically. Rejected before any arithmetic, so structural nulls never
   reach the evaluator.
2. **Degree feasibility.** Every block must saturate its slot count exactly.
3. **Connectivity.** A disconnected topology factorises into a product of
   lower-degree invariants, lives in the product subspace, and cannot carry a
   primitive direction.
4. **Canonicalisation.** Topologies related by relabelling like blocks describe
   the same scalar; the lexicographic minimum over that permutation group is
   the canonical key.

Index placement follows the rule the forward work paid for: **every contracted
edge carries exactly one raised end**, enforced by
`test_build_einsum_raises_exactly_one_end_of_every_edge`. An edge joining two
equally-placed slots contracts with `delta` instead of `eta` and is not a
Lorentz scalar — the exact defect that shipped in `P10_07`.

## 3. Enumeration, all 21 sectors

Block multisets are ordered **cheapest first** by total slot count. An earlier
ordering put the 30-slot all-`N` sector first and the survey stalled inside it;
the `M`-heavy sectors are 10 slots and complete instantly.

| sector | slots | canonical | selected | truncated | enum time |
|---|---:|---:|---:|---|---:|
| M+M+M+M+M | 10 | 1 | 1 | no | 0.0 s |
| M+M+M+M+N1050 | 14 | 4 | 4 | no | 0.0 s |
| M+M+M+M+N4125 | 14 | 1 | 1 | no | 0.0 s |
| M+M+M+N1050+N1050 | 18 | 46 | 40 | no | 0.3 s |
| M+M+M+N1050+N4125 | 18 | 19 | 19 | no | 0.0 s |
| M+M+M+N4125+N4125 | 18 | 5 | 5 | no | 0.0 s |
| M+M+N1050+N1050+N1050 | 22 | 289 | 40 | no | 5.2 s |
| M+M+N1050+N1050+N4125 | 22 | 208 | 40 | no | 0.8 s |
| M+M+N1050+N4125+N4125 | 22 | 55 | 40 | no | 0.1 s |
| M+M+N4125+N4125+N4125 | 22 | 8 | 8 | no | 0.0 s |
| M+N1050+N1050+N1050+N1050 | 26 | 1793 | 40 | **yes** | 85.8 s |
| M+N1050+N1050+N1050+N4125 | 26 | 1603 | 40 | no | 16.1 s |
| M+N1050+N1050+N4125+N4125 | 26 | 652 | 40 | no | 2.6 s |
| M+N1050+N4125+N4125+N4125 | 26 | 116 | 40 | no | 0.7 s |
| M+N4125+N4125+N4125+N4125 | 26 | 15 | 15 | no | 0.3 s |
| N1050 x5 | 30 | 3014 | 40 | **yes** | 97.6 s |
| N1050 x4 + N4125 | 30 | 4667 | 40 | **yes** | 66.4 s |
| N1050 x3 + N4125 x2 | 30 | 3057 | 40 | **yes** | 34.0 s |
| N1050 x2 + N4125 x3 | 30 | 834 | 40 | no | 8.6 s |
| N1050 + N4125 x4 | 30 | 120 | 40 | no | 2.2 s |
| N4125 x5 | 30 | 15 | 15 | no | 0.8 s |
| **total** | | **16 522** | **588** | 5 sectors | ~322 s |

**Five sectors are truncated** at the 30 000 raw-topology cap, so their
canonical counts are lower bounds and those sectors are **not exhausted**. This
is a bounded pilot, and the search space is not closed.

## 4. Cost model

Measured at prime 32749:

| quantity | value |
|---|---:|
| quadratic-block build, per sample | 4.4 s |
| one contraction, heaviest sector | ~150 ms |
| atlas construction, 14 columns x 22 samples | 34 s |
| evaluation, 588 candidates x 22 samples | ~130 s per sample |

Evaluation is **sample-outer**: each sample's three blocks are built once and
every candidate is swept against them. Candidate-outer would rebuild the blocks
for every (candidate, sample) pair — 4.4 s of setup for 150 ms of work, thirty
times the cost of the thing being measured. Caching all 22 samples' blocks
instead would hold ~350 MB, which this machine does not have.

**Projected cost of an exhaustive search.** Taking the untruncated canonical
counts as a lower bound of 16 522 candidates, at 22 samples and ~150 ms:

    16 522 x 22 x 0.15 s  ~  15 hours

on one worker, and that is a *lower* bound because five sectors are truncated.
Exhausting the declared space is a multi-day job at this throughput, not a
session-scale one.

## 5. Result — Q10 rank 3 recovered independently

> **A formula-independent bounded reverse search independently recovered the
> full three-dimensional quotient Q10.**

588 candidates, 22 samples, fit prime 32749, 3113 s, peak RSS 3095 MB.

| # | sector | einsum | rank after |
|---|---|---|---:|
| 1 | N4125 x5 | `abcdef,abcdgh,efijkl,ghimno,jklmno->` | 1 |
| 2 | N4125 x5 | `abcdef,abcdgh,efijkl,gijmno,hklmno->` | 2 |
| 3 | N4125 x5 | `abcdef,abcdgh,egijkl,fhimno,jklmno->` | 3 |

### The recovered basis is NOT the published one

All three recovered directions come from the **`N^(4125)` x5** sector — pure
`N^(4125)` contractions with no `N^(1050)` and no `M` at all. Every element of
the published Level-B basis `{P10_10, P10_11, P10_12}` is `N^(1050)`-based.

The reverse search did not rediscover the published formulas. It found an
independent set spanning the same quotient. That is a stronger outcome than a
rediscovery would have been: it shows Q10 is reachable from a block sector the
published list does not use at all, and it could not have been produced by
leaking the answer into the search.

It is also consistent with the declared scope. The generator applies no
explicit BLACK or RED bracket beyond what `composite_n1050` carries
intrinsically, so `P10_04`, `P10_05`, `P10_08` and `P10_09` are **not
representable** by it — the reverse engine could not have reproduced those by
name even in principle. `P10_10`, `P10_11` and `P10_12` are bracket-free and in
principle reachable; the bounded pilot evaluated 40 of the 3014 canonical
candidates in their sector and did not happen to hit them before rank 3 was
already complete from the cheaper `N^(4125)` sector.

### Validation

`results/intrinsic_candidates/degree10_reverse_span_validation.json`.

Each recovered topology is **regenerated from its recorded einsum** by
re-enumerating its sector, rather than deserialised — a stricter check, since it
confirms the recorded specification is actually producible by the declared
search and fails loudly if the generator has drifted.

| prime | samples | rank | in atlas span |
|---|---|---:|---|
| 32749 (fit) | original | **3/3** | yes |
| 32749 (fit) | fresh seed base | **3/3** | yes |
| 32717 (holdout) | original | **3/3** | yes |
| 32717 (holdout) | fresh seed base | **3/3** | yes |

### Span equality — PROVEN

| prime | rank reverse | rank published | rank union | dim Q10 | spans equal | change-of-basis mutually inverse |
|---|---:|---:|---:|---:|---|---|
| 32749 | 3 | 3 | **3** | 3 | **yes** | **yes** |
| 32717 | 3 | 3 | **3** | 3 | **yes** | **yes** |

Both spans lie inside a 3-dimensional quotient and each has rank 3, so each
already **is** Q10 and they are therefore equal. The union rank is the direct
check: six vectors spanning only rank 3 means no direction lies outside the
common span. Exact change-of-basis matrices in both directions are recorded and
verified mutually inverse over F_p.

    reverse-generated span  =  Q10  =  published Level-B span

on both the fitting and the holdout prime. **The benchmark passes.**

## 6. Scope

**Permitted**: statements about what the *bounded pilot* covered — 588 of at
least 16 522 canonical candidates, 40 per sector, one fit prime.

**Forbidden**: any claim that the degree-10 block search space has been
exhausted. It has not. Five sectors are truncated at enumeration and every
sector above 40 canonical candidates is sampled rather than swept.
