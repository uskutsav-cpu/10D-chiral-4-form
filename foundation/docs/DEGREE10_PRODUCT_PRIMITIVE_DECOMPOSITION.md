# The degree-10 product / primitive decomposition — computed, not assumed

## 1. Result

The 14-dimensional degree-10 atlas `A10` carries **two explicit product
columns**, present in the registry by construction:

    A10 = span( I10_1 ... I10_12,  I4_1*I6_1,  I4_1*I6_2 )
    P10 = span( I4_1*I6_1, I4_1*I6_2 )          dim P10 = 2

These are the only products of total field degree 10: the Hilbert series has no
degree-2 invariant and begins at `b_4 = 1`, so the only splitting is
`4 + 6`, giving `1 x 2 = 2` product monomials.

**The published-candidate span is NOT a primitive complement.** Exact modular
computation at both validation primes:

| quantity | 32749 | 32717 |
|---|---:|---:|
| `dim B10` (published span) | 12 | 12 |
| `dim P10` (products) | 2 | 2 |
| `dim (B10 + P10)` | **13** | **13** |
| `dim (B10 ∩ P10)` | **1** | **1** |

Since `dim(B10 + P10) = 13 < 14`, the published span together with the products
does **not** span the atlas, and since `dim(B10 ∩ P10) = 1`, the published span
**contains exactly one product direction**.

Therefore the primitive content of the published span is

    12 - 1 = 11        NOT 12

## 2. The coincidence that had to be checked

Two decompositions both read `12 + 2`:

    14 = 12 + 2     published span B10 + complement        (trace)
    14 = 12 + 2     primitive generators + products        (spinor Hilbert data)

It was tempting to identify them and claim the published candidates span exactly
the primitive part. **That claim is false.** The published span carries one
product direction and misses one atlas direction entirely; the agreement of the
two `12`s is **dimensional coincidence, not structural identity**.

This is why the dimension gate exists. Asserting the identification would have
put a wrong structural claim into the manuscript while every individual number
remained correct.

## 3. Permitted wording

> The twelve published equation-(4.24) candidates span a 12-dimensional
> subspace `B10` of the 14-dimensional degree-10 atlas. `B10` meets the
> 2-dimensional product subspace `P10` in exactly one dimension, so `B10` has
> primitive content 11 and is not a complement of `P10`.

**Forbidden**: calling `B10` "the primitive part"; calling its dimension "12
primitives"; identifying the trace `12` with the spinor primitive count `m_10`.

## 4. Relation to Q10 — a third, independent splitting

    14 = 11 + 3     reachable closure D10 + quotient Q10

`Q10` is a quotient by the **stress-flow reachable** subspace, not by products.
Nothing about `Q10` follows from the product/primitive decomposition, and the
numerical coincidence that `dim D10 = 11` equals the primitive content of `B10`
is again **not** a structural statement — the two 11s describe different
subspaces.

## 5. The correspondence IS real — but with G10, not B10

Completing the pairwise incidence table (`degree10_space_incidence.json`)
identifies the correct primitive complement. Writing `G10` for the span of the
twelve graph generators `I10_1 ... I10_12`:

| pair | dim A | dim B | sum | ∩ | relation |
|---|---:|---:|---:|---:|---|
| `P10`, `G10` | 2 | 12 | **14** | **0** | `G10 ⊕ P10 = A10` |
| `P10`, `B10` | 2 | 12 | 13 | 1 | neither |
| `P10`, `D10` | 2 | 11 | 11 | 2 | **`P10 ⊆ D10`** |
| `B10`, `G10` | 12 | 12 | 14 | 10 | distinct 12-spaces |
| `D10`, `B10` | 11 | 12 | 14 | 9 | neither |

So:

    A10 = G10 (+) P10          12 + 2, a genuine direct sum
    A10 = D10 (+) (Q10 lift)   11 + 3, and P10 is contained in D10
    B10                        a third 12-dimensional subspace, equal to neither

**The structural `12 + 2` correspondence with the spinor Hilbert decomposition
holds for `G10`, the graph generators — not for `B10`, the published span.**
`G10` meets the product subspace trivially and complements it exactly, which is
what "primitive" means. `B10` does not.

`P10 ⊆ D10` is the expected statement and is now verified rather than assumed:
products of reachable lower-degree invariants are themselves reachable, so the
products contribute nothing to `Q10`.

## 6. Final permitted wording

> The degree-10 atlas decomposes as `A10 = G10 ⊕ P10` with `dim G10 = 12` and
> `dim P10 = 2`, matching the spinor-side primitive/product counts. The
> published candidate span `B10` is a distinct 12-dimensional subspace which
> meets `P10` in one dimension and `G10` in ten, and is therefore not a
> primitive complement. The product subspace lies inside the reachable closure,
> `P10 ⊆ D10`, so products contribute nothing to `Q10`.
