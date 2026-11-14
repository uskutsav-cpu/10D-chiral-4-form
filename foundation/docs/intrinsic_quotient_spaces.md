# The intrinsic quotient spaces Q_10 and Q_12

**Step 1 of Phase 1: complete.** The quotient spaces are constructed exactly and
their dimensions are certified. Step 3–4 — exhibiting *intrinsic* spanning
classes — is **not** complete; see §4 for a bounded negative.

Reproduce:

    .venv/bin/python scripts/solve_intrinsic_quotients.py --degree 10 \
        --out results/generalized_flow/quotient_degree10.json
    .venv/bin/python scripts/solve_intrinsic_quotients.py --degree 12 \
        --out results/generalized_flow/quotient_degree12.json

## 1. Construction

For each degree d, with A_d the full homogeneous invariant space and D_d the
reachable seed closure:

    Q_d = A_d / D_d

Exactly, over F_p: reduce D_d to row echelon form recording pivot columns P;
the quotient coordinates are the non-pivot columns; the projection
`pi : A_d -> Q_d` reduces a vector against the echelon rows and reads the
non-pivot entries.

## 2. Results

| degree | dim A_d | dim D_d | **dim Q_d** | primes agreeing |
|---:|---:|---:|---:|---:|
| 10 | 14 | 11 | **3** | 6 / 6 |
| 12 | 72 | 68 | **4** | 6 / 6 |

Primes: 32749, 32719, 32717, 32693, 32771, 32713.

**dim Q_d is the basis-independent content.** The representatives that the
echelon form happens to select —

    Q_10 representatives: I10_6, I10_7, I10_12
    Q_12 representatives: I12_59, I12_60, I12_61, I12_62

— are **coordinate labels of one graph basis and are not intrinsic**. They are
convenient handles on the classes, nothing more. This is the same distinction
that separates `I6_2` from `K6`.

## 3. A structural constraint on any intrinsic representative

Established directly from the committed bases, and it sharply narrows the
search:

**No quotient direction is a product.** Every product entry in both bases lies
*inside* D_d:

| degree | product entries in the basis | in closure? |
|---|---|---|
| 10 | `I4_1*I6_1`, `I4_1*I6_2` | yes, both |
| 12 | `I4_1^3`, `I6_1^2`, `I6_1*I6_2`, `I6_2^2`, `I4_1*I8_1..6` | yes, all ten |

None appears in a missing set. Therefore **the quotient classes are genuinely
primitive**, and an intrinsic representative cannot be built as a product of
lower-degree intrinsics.

This eliminates the entire "lower-intrinsic products" family from the Step 2
candidate library *before* any evaluation — a real narrowing, obtained from
data already committed rather than from new computation.

## 4. Bounded negative: the quotients are NOT yet intrinsically spanned

Stating this precisely, per the completion gate's requirement for a bounded
negative rather than a fabricated formula.

**Families searched so far** (by structural argument, not evaluation):

| family | status at degrees 10 and 12 |
|---|---|
| products of lower intrinsics (I4, J6, K6, I8_k) | **excluded** — all product basis entries lie in D_d (§3) |
| stress-adapted traces `Tr(M^k)` and their products | **not excluded, not yet tested**; the static stress span has dimensions 2 and 4 at degrees 10 and 12 and its relation to D_d has not been computed |
| N^(1050)-type primitive contractions (the K6 family) | **not searched** — no degree-10/12 generator implemented |
| new M–F mixed contraction channels | **not searched** |
| epsilon / parity-odd contractions | **not searched** |
| flow-generated (Lie derivative, bracket) structures | **not searched** |

**Quotient rank reached intrinsically: 0 of 3 at degree 10, 0 of 4 at
degree 12.** No intrinsic tensor expression has been produced for any of the
seven classes.

**What is needed next.** A generator for primitive scalar contractions at total
field degree 10 and 12 built from N^(1050) and M, analogous to the eight-term
graph expansion that defines K6 at degree 6. That is new machinery: the K6
expansion required a specific tensor ansatz plus verification, and the degree
10/12 analogues have a substantially larger contraction-channel space.

**This is the mathematical-insight step flagged in `docs/RESOURCE_PLAN.md`.**
It is not compute-bound, and no amount of CPU on the existing machinery
advances it. Fabricating a compact formula that happens to fit four samples is
explicitly forbidden by the gate and would be worse than the honest negative.

## 5. Caveats carried

- D_d is a **seed** closure. The generator-extension problem
  `dV/dlambda = f(T, S, lambda)` is different and is untouched (Step 8).
- Quotient dimensions are `MOD-CERT`: identical at six primes, which is strong
  evidence but not a characteristic-zero proof (PO-09).
- General GL basis-change validation remains **open** (PO-08 partial); only the
  permutation subgroup has been verified.
