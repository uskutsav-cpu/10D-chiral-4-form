# Generalized flow: the complete degree-12 deficit

**Phase 1, partial.** The deficits are located and shown minimal. The Phase 1
gate is **NOT met**, because the gate requires every added direction to be
*intrinsic* and these are graph labels. See §5.

Reproduce:

    .venv/bin/python scripts/find_missing_flow_directions.py --degree 10 \
        --out results/generalized_flow/degree10_missing_directions.json
    .venv/bin/python scripts/find_missing_flow_directions.py --degree 12 \
        --single-prime-scan \
        --out results/generalized_flow/degree12_missing_directions.json

## 1. The deficits, resolved

Starting from the free seed with K6 and the four degree-8 directions adjoined,
the closure is `(1, 2, 7, 11, 68)` against full `(1, 2, 7, 14, 72)`.

| degree | deficit | missing directions | closure with all |
|---:|---:|---|---:|
| 10 | 3 | `I10_6`, `I10_7`, `I10_12` | 14 (full) |
| 12 | 4 | `I12_59`, `I12_60`, `I12_61`, `I12_62` | 72 (full) |

Each is non-redundant under removal, on every prime tested
(32749, 32719, 32717, 32693):

- degree 10: dropping any one gives 13, not 14;
- degree 12: dropping any one gives 71, not 72.

The scan is exhaustive over the basis at each degree — every basis direction is
tested individually — so these are *the* missing sets, not a sufficient set
found by search.

## 2. The full minimal completion through degree 12

| degree | directions | count |
|---:|---|---:|
| 6 | K6 | 1 |
| 8 | I8_3, I8_4, I8_5, I8_6 | 4 |
| 10 | I10_6, I10_7, I10_12 | 3 |
| 12 | I12_59, I12_60, I12_61, I12_62 | 4 |
| | **total** | **12** |

## 3. The closure is coordinate-aligned — and that is not generic

At degree 12 the closure has dimension 68 inside a 72-dimensional space, and
**exactly 68 of the 72 basis vectors lie inside it**, with exactly 4 outside.

A generic 68-dimensional subspace of a 72-dimensional space contains **no**
basis vector at all. Containing 68 of them means the closure is precisely the
coordinate subspace spanned by all basis elements except those four. The same
holds at degree 10 (11 of 14).

This is a genuine structural fact about the graph basis, and it is the reason
the deficits can be named by labels at all. It is also **basis-dependent as
stated**: under a generic invertible change of basis the closure remains a
68-dimensional subspace but stops being coordinate-aligned. The invariant
content is the 4-dimensional quotient, not the labels.

## 4. An unexplained coincidence worth recording

The atlas records exactly two degree-12 candidates that are **functionally
dependent** — they do not raise the cumulative Jacobian rank of 81:

    results/10d_order12.json → discovery.functional_dependencies = ["I12_61", "I12_62"]

Both of them are in the unreachable set. The overlap:

| | |
|---|---|
| functionally dependent | I12_61, I12_62 |
| dynamically unreachable | I12_59, I12_60, I12_61, I12_62 |
| dependent **and** unreachable | I12_61, I12_62 (2 of 2) |
| unreachable but **not** dependent | I12_59, I12_60 |

So on the available data, *functional dependence implies dynamical
unreachability* (2 out of 2), while the converse fails (4 unreachable, only 2
dependent).

**This is an observation, not a theorem.** Two data points is weak evidence,
and no mechanism is proposed. It is recorded as conjecture CJ-01 for Phase 4,
where the geometric framework must either explain it or expose it as
coincidence. It must not appear in any manuscript in its current state.

## 5. Why the Phase 1 gate is NOT met

The gate requires: *"every added direction is intrinsic."*

`I10_6`, `I12_59` and the rest are **graph labels in a particular basis** —
exactly the situation `I6_2` was in before it was replaced by the intrinsic
`K6 = N1050…` contraction. A label is not a tensor.

Remaining Phase 1 work:

1. an intrinsic tensor expression for each of the seven new directions;
2. minimality under basis change (PO-08), which the falsification list also
   demands (test 4);
3. re-running the scan under a reordered candidate list (falsification test 5)
   to confirm the sets are not an ordering artifact;
4. the distinction in §6.

## 6. Standing caveat: seeds are not generators

Everything here enlarges the **seed**. Adjoining `S` as an independent
**generator** of `f(T, S, λ)` is a different operation requiring certificate
rows that do not exist. The seeding result bounds the generator result from one
side only. This caveat is embedded in both JSON artifacts, not only in prose,
because it is the distinction most easily lost when the data is read later.
