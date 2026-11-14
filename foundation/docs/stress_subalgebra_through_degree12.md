# The static free-stress subalgebra through degree 12

Source artifact: `results/stress_flow/dimension_table.json`
Gate: `tests/test_static_stress_degree12.py`

Scope, quoting the artifact itself: *"Static scalar invariants of the free
traceless INZ stress tensor T = M/48; interacting formal stress maps are
separate."* That separation matters — see
`docs/stress_flow_classification.md`, where the dynamical closure turns out
to be much larger than the static span computed here.

## 1. Dimensions

| degree | full invariant space | free-stress span | quotient |
|---:|---:|---:|---:|
| 4 | 1 | 1 | 0 |
| 6 | 2 | 1 | 1 |
| 8 | 7 | 2 | 5 |
| 10 | 14 | 2 | 12 |
| 12 | 72 | 4 | 68 |

The degree-12 stress basis is `(tr_M6, tr_M4*tr_M2, tr_M3^2, tr_M2^3)`.
Three of those four are analytic products of lower-degree traces; only
`tr_M6` required a full reduction against the committed 72-element atlas.

## 2. Exactness record

| property | value |
|---|---|
| primes | 15: 32603 … 32771 |
| samples per prime | 4 (seeds 20260729–20260732) |
| basis rank at every prime | 72 |
| stress rank at every prime | 4 |
| same standard complement at every prime | yes |
| degree-12 atlas sha256 | `9a784dc5…5113` |
| engine sha256 | `a3360be6…38ef` |

The pivot columns of the stress rows inside the original atlas basis are
`['I4_1^3', 'I6_1^2', 'I4_1*I8_1', 'I12_58']` at indices `[0, 1, 4, 67]`.

## 3. The Tr(M^6) rational reconstruction FAILED — and that is the result

This is the most important methodological point in the static map, and it
should not be softened in any write-up.

The original plan was: fit `Tr(M^6)` modulo several primes, reconstruct
rational coefficients by CRT, then validate on an independent holdout prime.
With **15 primes** — CRT modulus

    51681124540149423589916523149247886541408185180838373465166625800747

roughly `5.2e67` — **29 of the 72 columns still fall outside the certified
uniqueness bound.** The columns that did lift came back with heights like

    -4639872259700851804555423085146127 / 719951961948631043048952431460182

which are reconstruction artifacts of a modulus that is still too small, not
plausible physical coefficients.

The recorded policy is the correct one:

> No rational coefficient is guessed outside the certified CRT uniqueness
> bound. Tr(M^6) is used as an intrinsic stress-adapted basis element
> instead.

So the degree-12 row is **not** expressed in the atlas basis with certified
rationals. Instead `Tr(M^6)` itself is adopted as a basis element of a
stress-adapted basis, which is exact, intrinsic, and requires no lifting.

Two consequences for anyone writing this up:

1. Do not quote rational coordinates for `Tr(M^6)` in the atlas basis. They
   are not certified, and the ones in
   `rational_lift_audit.successfully_reconstructed_columns` are certified
   only in the sense that CRT *terminated*, not that the answer is right.
2. The statement "the degree-12 free-stress span is 4-dimensional" **is**
   certified — it is a rank, computed identically at 15 primes, and ranks do
   not require lifting.

Closing the gap would need either many more primes (the height suggests a
substantially larger budget) or an analytic identity for `Tr(M^6)` in terms
of the atlas generators. The latter is the better route and is listed as an
open question.

## 4. Physical normalisation

The artifact records `physical_free_stress_coordinates` alongside the raw
coordinates, related by a factor `48^6` per the `tau = 48*T` convention. Any
comparison with published formulas must fix this normalisation first; see
`docs/interacting_stress_tensor.md`.
