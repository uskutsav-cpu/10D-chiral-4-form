# Immediate Execution Checklist

This is the shortest route from repository setup to the first genuinely new Paper-2 result.

## Week 1: foundation freeze

- [ ] Clone/check out the predecessor repository beside this one.
- [ ] Pin its exact commit and branch.
- [ ] Run `python scripts/freeze_source_snapshot.py ../selfdual-5form-invariants`.
- [ ] Re-run predecessor tests and the interacting-flow closure artifact.
- [ ] Add a new certificate here for `(1,1,3,11,67)` using fresh primes/seeds.
- [ ] Update C01–C04 from `IMPORTED` to `REPRODUCED` only after the new artifacts exist.

## Weeks 1–3: K6 theorem

- [ ] Reproduce the exact `(Tr(M^3), K6)` change of basis.
- [ ] Write the degree-6 generator-exhaustion table.
- [ ] Derive the `K6` projection of each eligible generator analytically.
- [ ] Produce a symbolic derivation of `q6_dot = C q6`.
- [ ] Check the coefficient `40` against the pinned normalization.
- [ ] Add a negative-control test showing a deliberately illegal generator can create an inhomogeneous `K6` term.

**First paper-grade milestone:** C05 and C06 become `ESTABLISHED`.

## Weeks 3–6: higher obstruction certificates

- [ ] Export exact reachable-row coordinates for degrees 8, 10, 12.
- [ ] Run the annihilator pipeline.
- [ ] Verify quotient dimensions under unrelated basis changes.
- [ ] Search for intrinsic tensor representatives of `Q8` first.
- [ ] Record any failure to find compact representatives rather than overclaiming graph labels.

## Weeks 5–8: true generalized-flow completion

- [ ] Implement the metric/stress variation induced by a candidate `S_i` inside the flow law.
- [ ] Compare generator augmentation against seed augmentation explicitly.
- [ ] Compute `Delta_d(S)`.
- [ ] Run removal and smaller-cardinality searches.

**Second paper-grade milestone:** C08 becomes either `ESTABLISHED` in a declared class or is replaced by a weaker, correct theorem.

## Weeks 7–10: conformal / ModMax

- [ ] Introduce algebraic localization (`r^2=I4` or equivalent homogeneous coordinates).
- [ ] Verify the polynomial limit reproduces known closure.
- [ ] Encode the HLS ModMax-like trajectory.
- [ ] Project onto intrinsic obstruction coordinates.
- [ ] Decide C11 versus C12 with a certificate.

## Weeks 10–12: freeze core paper

- [ ] Repeat novelty search.
- [ ] Ask mentor to resolve convention/sign questions.
- [ ] Run untouched holdout inputs.
- [ ] Freeze machine artifacts and manuscript tables.
- [ ] Attempt degree-14 probe only after the core result is safe.
