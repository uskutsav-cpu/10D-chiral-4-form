# Exact scope of the degree-10 reverse search

## 1. Two goals that must not be conflated

**Recovery goal** — recover Q10 rank 3 independently of the published formulas.
Satisfied as soon as three independent quotient directions are produced by
formula-independent generation and validated on a holdout prime. This is the
goal the current work targets.

**Exhaustion goal** — enumerate every candidate in a precisely defined
block-contraction class. **Not reached.** Five sectors are truncated at their
raw-topology cap.

A recovery result says *the quotient is reachable from block topology alone*.
An exhaustion result would additionally say *and here is everything that
reaches it*, which is what a minimality theorem would need. Only the first is
claimed.

## 2. What is inside the declared class

| dimension | restriction |
|---|---|
| block families | `M`, `N^(1050)`, `N^(4125)` only |
| block count | exactly 5 (degree 2 each, so field degree 10) |
| parity sector | self-dual 5-form `F` in D=10, Lorentzian signature |
| free indices | none; every candidate is a full scalar contraction |
| index placement | exactly one raised end per contracted edge |
| arithmetic | exact modular, `mod_einsum`; no floating point anywhere |

### Contraction-topology rules

A contraction is an edge multiset over `(block, slot-class)` endpoints.
Slot classes are:

| block | class A | class B |
|---|---|---|
| `M` | 2 slots | — |
| `N^(1050)` | 5 slots, antisymmetric | 1 slot |
| `N^(4125)` | 6 slots | — |

Every endpoint must be saturated exactly. This is what makes the search
tractable: raw slot matchings for five `N` blocks number 29!! ≈ 6.2×10¹⁵, while
the class-level edge-count enumeration is a small integer problem.

### Reductions and prunings applied

1. **Forbidden traces** — an edge joining two class-A slots of the same block
   contracts inside an antisymmetric group and vanishes identically. Rejected
   before arithmetic.
2. **Young symmetry** — used implicitly, in that `composite_n1050` carries the
   equation-(2.15) antisymmetrisation intrinsically and the class structure is
   read off it. No additional Young projector is applied on top.
3. **Automorphism / relabelling reduction** — topologies related by permuting
   like blocks are identified by a canonical key (lexicographic minimum over
   the permutation group of equal-kind blocks).
4. **Connectivity** — disconnected topologies factorise into products of
   lower-degree invariants, live in the product subspace, and cannot carry a
   primitive direction.

## 3. What is OUTSIDE the declared class

Recorded so the boundary is not mistaken for a completeness claim:

- **Higher blocks.** Only quadratic blocks are considered. A degree-10 scalar
  built from, say, a quartic block that is not a product of two quadratics is
  not reachable here.
- **Explicit bracket structures.** The reverse engine applies no additional
  BLACK or RED (anti)symmetrisation beyond what `composite_n1050` carries
  intrinsically. Published candidates that need an explicit bracket program —
  `P10_04`, `P10_05`, `P10_08`, `P10_09` — are therefore **not representable**
  by this generator at all. This is a real limitation of the class, not an
  oversight, and it means the reverse search cannot be expected to reproduce
  those formulas by name.
- **Non-generic index orderings.** The generator hands out slots within a class
  in order. Because class A is antisymmetric, different orderings differ at
  most by sign, so this loses no span — but it does mean signed variants are
  not enumerated separately.
- **Sign/orientation variants.** Canonicalisation identifies topologies related
  by block relabelling without tracking the induced sign. Two topologies that
  differ only by an odd permutation of antisymmetric slots are treated as one
  candidate. This can only lose an overall sign, never a direction, but it is
  recorded because it is an identification the code makes.

## 4. Sector-by-sector exhaustion status

Cap: 30 000 raw topologies per sector.

| sector | canonical | exhausted? |
|---|---:|---|
| M+M+M+M+M | 1 | yes |
| M+M+M+M+N1050 | 4 | yes |
| M+M+M+M+N4125 | 1 | yes |
| M+M+M+N1050+N1050 | 46 | yes |
| M+M+M+N1050+N4125 | 19 | yes |
| M+M+M+N4125+N4125 | 5 | yes |
| M+M+N1050+N1050+N1050 | 289 | yes |
| M+M+N1050+N1050+N4125 | 208 | yes |
| M+M+N1050+N4125+N4125 | 55 | yes |
| M+M+N4125+N4125+N4125 | 8 | yes |
| M+N1050 x4 | 1793 | **NO — capped** |
| M+N1050 x3+N4125 | 1603 | yes |
| M+N1050 x2+N4125 x2 | 652 | yes |
| M+N1050+N4125 x3 | 116 | yes |
| M+N4125 x4 | 15 | yes |
| N1050 x5 | 3014 | **NO — capped** |
| N1050 x4+N4125 | 4667 | **NO — capped** |
| N1050 x3+N4125 x2 | 3057 | **NO — capped** |
| N1050 x2+N4125 x3 | 834 | yes |
| N1050+N4125 x4 | 120 | yes |
| N4125 x5 | 15 | yes |

**16 of 21 sectors exhausted; 5 capped.** Total canonical candidates retained:
16 522, which is a **lower bound** on the class size.

Within each sector the pilot evaluated at most 40 candidates, so even the
exhausted sectors are *sampled*, not swept, at the evaluation stage.

## 5. Claims

**Permitted**: "A formula-independent bounded reverse search independently
recovered the full three-dimensional quotient Q10."

**Forbidden**:

- "Complete enumeration of every M/N contraction."
- "The reverse search proves the published basis is minimal."
- "All degree-10 compact structures reaching Q10 are known."

## 6. Cost of closing the exhaustion goal

At ~150 ms per contraction and 22 samples, the retained 16 522 canonical
candidates alone cost ≈ 15 h on one worker, and the five capped sectors would
add an unknown multiple. Streamed generation with `--sector`, `--shard-residue`
and `--shard-modulus` makes this parallelisable across runs; the coordinator's
memory is now bounded by distinct canonical keys rather than by raw topology
count (33 MB against the 2611 MB the list form reached).
