# Lorentz invariants of a self-dual 5-form in 10D

**Science completion:** an independent D5 character calculation now derives the
Hilbert dimensions through degree 22, closes degree-12 completeness, and removes
the literature-count premise from the corrected source maps. See
[the completion ledger](SCIENCE_STATUS.md) and the minimal
[standalone rank/orbit verifier](science/rank81/README.md).

**September 2026 graph-certificate release:** the same explicit 81 graphs now
have fresh nonzero 81×81 Jacobian minors at two points under each of two
primes. See [the roadmap and proof scope](docs/CLASSIFICATION_ROADMAP.md),
[the ordered basis](results/rank81_basis.json),
[the complete witnesses](results/rank81_certificate.json), and
[the manuscript with all formulas](paper/manuscript.pdf).
The independent audit preserved a failure in the frozen literature
transcription and repaired the J10–J12 antisymmetrizers. See
[the audit report](INDEPENDENT_AUDIT.md) for the frozen failure, corrected
source maps, exact evidence and remaining assumptions. Independent fresh
checks verify every Jacobian entry at two further points, the degree-12
polynomial rank 72, and generic orbit dimension 45. Run
`python scripts/verify_rank81.py --recompute` to re-evaluate the four frozen
graph witnesses, and `python scripts/verify_independent_audit.py` to check
the saved independent certificates using integer and rational arithmetic.

Extension of Elamaran–Ferko–Scarlett, *Machine Learning Invariants of Tensors*
([arXiv:2512.23750](https://arxiv.org/abs/2512.23750), Phys. Rev. D).

## The question

How many functionally independent Lorentz scalars can be built from a
self-dual 5-form $F^+_{\mu_1\ldots\mu_5}$ in ten dimensions, and what are they
explicitly?

The **count is known: 81** (Hutomo–Lechner–Sorokin,
[arXiv:2509.14351](https://arxiv.org/abs/2509.14351), JHEP 02 (2026) 147;
structural analysis in Cederwall–Hutomo–Kuzenko–Lechner–Sorokin,
[arXiv:2509.14350](https://arxiv.org/abs/2509.14350)).

The degree-by-degree partition function is also known through order 22. In
particular,

> $P(t)=1+t^4+2t^6+7t^8+\cdots
> =(1-t^4)^{-1}(1-t^6)^{-2}(1-t^8)^{-6}\cdots$.

Thus there is one quartic generator, two sextic generators, and six octic
generators; the seventh degree-8 scalar is the product $I_4^2$. The literature
gives a tensor basis for the octic generators. This repo supplies an explicit
contraction-graph basis and independently checks its Jacobian rank.

On the regular local quotient, these scalars give local coordinates for
Lorentz-invariant functions of $F_5$. This supports local descriptions of
interactions without derivatives, including chiral 4-form theories. Global
orbit separation and a presentation of the full polynomial invariant ring
are not established here.

## Status

| Case | Result |
|---|---|
| 6D generic 3-form (reproduction) | **PASS** — 5 invariants, pattern 1,2,1,1 at orders 2,4,6,8, confirmed under two primes |
| 10D self-dual 5-form, order 4 | **1** independent invariant (all 4 candidate graphs enumerated — complete) |
| 10D, order 6 | **2** new independent invariants, running rank **3 / 81** (49 exact graph classes) |
| 10D, order 8 | **6** new independent invariants, running rank **9 / 81** (1,689 exact graph classes; complete under two primes) |
| 10D, order 10 | **12** new connected primitive directions, running rank **21 / 81**; degree-10 value rank **14** after adding $I_4I_{6,1}$ and $I_4I_{6,2}$ (187,392 exact graph classes; two primes, three Jacobian samples per prime) |
| 10D, order 12 | **62** connected primitive polynomial directions plus **10** lower products give degree-12 rank **72**; 60 add functional directions, giving the full cumulative rank **81 / 81**. The remaining `I12_61` and `I12_62` are linearly independent homogeneous directions modulo lower products but add no cumulative functional direction (three primes, four samples per prime). |
| Exact stress-flow map through degree 10 | **PASS** — free-stress dimensions **1,1,2,2** inside five-form value dimensions **1,2,7,14**; exact ModMax $I_8/I_{12}$ reproduction under three primes |

Artifacts:

- [`results/10d_order8.json`](results/10d_order8.json) — the nine graph
  generators, two-prime rank logs, and completeness claim.
- [`results/10d_graph_catalog.json`](results/10d_graph_catalog.json) — all
  4, 49, and 1,689 exact isomorphism classes at orders 4, 6, and 8.
- [`results/10d_baseline.json`](results/10d_baseline.json) — the original
  order-4/order-6 baseline.
- [`results/10d_order10.json`](results/10d_order10.json) — the twelve
  explicit degree-10 graph generators, the two product directions, and exact
  two-prime validation evidence.
- [`results/degree10_benchmarks.json`](results/degree10_benchmarks.json) —
  stage distributions and before/after optimization measurements.
- [`docs/degree10.md`](docs/degree10.md) — catalog, checkpoint, reproduction,
  benchmark, and limitation details.
- [`results/10d_order12.json`](results/10d_order12.json) — all 62 explicit
  order-12 graphs, ten product directions, 83-candidate primitive inventory,
  and the 3-prime by 4-sample exact certificate.
- [`results/degree12_benchmarks.json`](results/degree12_benchmarks.json) —
  generation, planning, contraction, rank, checkpoint, width, and RSS
  distributions.
- [`docs/degree12.md`](docs/degree12.md) — proof scope, exact shard hashes,
  staged reproduction, memory bounds, and clean-checkout revalidation.
- [`results/stress_flow_exact_low_degree.json`](results/stress_flow_exact_low_degree.json)
  — five-prime rational change-of-basis matrices, exact complement
  certificates, ModMax regressions, and the flow-closure pilot.
- [`docs/stress_flow.md`](docs/stress_flow.md) — paper conventions,
  normalization, maps, obstructions, perturbative scope, and degree-12 import
  contract.

### The degree-12 result

The homogeneous degree-12 scalar space has exact rank **72** on the committed
basis:

- $I_{4,1}^3$;
- $I_{6,1}^2$, $I_{6,1}I_{6,2}$, and $I_{6,2}^2$;
- $I_{4,1}I_{8,k}$ for $k=1,\ldots,6$; and
- 62 explicit connected contractions
  $I_{12,1},\ldots,I_{12,62}$.

Polynomial independence is certified by exact ranks of gradients stacked
across four independent generic points. This matters: at one point the ten
product gradients lie in the lower tangent span and cannot themselves have
rank ten. Across four points, the ten products and all 62 connected
directions are pivots, giving rank 72. Attaining the supplied Hilbert-series
upper bound proves degree-12 completeness.

For functional independence, the 21 lower primitive rows and the first 60
order-12 rows give rank **81**. `I12_61` and `I12_62` remain new homogeneous
polynomials but their Jacobian rows reduce to zero against that cumulative
basis. This exact pattern repeats at all four seeds under each of primes
32749, 32719, and 32693. All selected values are nonzero in those runs, and
all 62 graphs pass an exact SO(1,9) boost check under every prime.

### The degree-10 result

The exact degree-10 value space has rank **14** on the saved basis:
twelve connected contractions $I_{10,1},\ldots,I_{10,12}$ plus the two lower
products $I_{4,1}I_{6,1}$ and $I_{4,1}I_{6,2}$. The twelve graph formulas and
canonical SHA-256 IDs are explicit in
[`results/10d_order10.json`](results/10d_order10.json).

At each of seeds 20260729, 20260730, and 20260731, and under both primes 32749
and 32719, the lower generators have cumulative Jacobian rank 9 and the twelve
degree-10 rows raise it to **21**. Separately, sixteen value samples per prime
give exact rank **14** for the homogeneous degree-10 basis. Reaching the
published upper bound of twelve new primitives proves completeness at this
degree; it does not claim a polynomial generating set through all degrees.

### The six order-8 graph generators

`nij^m` means that tensor copies `i` and `j` share `m` contracted indices.
Every omitted pair has multiplicity zero.

| ID | Exact contraction graph |
|---|---|
| $I_{8,1}$ | `n8[04^4,05^1,14^1,16^4,25^4,27^1,36^1,37^4]` |
| $I_{8,2}$ | `n8[03^4,06^1,14^4,16^1,25^4,26^1,37^1,47^1,57^1,67^2]` |
| $I_{8,3}$ | `n8[03^4,04^1,14^1,15^1,16^2,17^1,25^4,26^1,37^1,46^1,47^2,67^1]` |
| $I_{8,4}$ | `n8[03^2,06^2,07^1,14^4,16^1,25^4,26^1,36^1,37^2,47^1,57^1]` |
| $I_{8,5}$ | `n8[03^2,04^1,06^2,14^2,15^1,17^2,25^4,26^1,36^2,37^1,47^2]` |
| $I_{8,6}$ | `n8[03^2,06^2,07^1,14^2,15^1,16^2,24^1,25^2,27^2,36^1,37^2,45^2]` |

These six Jacobian rows add six directions to the quartic/sextic basis under
both primes 32749 and 32719. This meets the published upper bound of six new
octic generators, so the list is complete. The remaining degree-8 scalar is
the disconnected product $I_4^2$.

### Exact order-8 enumeration

The catalog is generated with nauty 2.9.3:

```bash
geng -cq 8 | multig -q -T -m4 -r5
```

`geng` emits each connected underlying simple graph once. `multig` assigns
edge multiplicities, enforces degree 5 and maximum multiplicity 4, and
suppresses isomorphic weighted outputs. The result is exactly **1,689**
connected candidates. `pynauty` supplies an independent exact certificate in
the tests; no Weisfeiler-Lehman hash is used for correctness.

### Why the tests use a boost, not just a rotation

`test_10d_contractions_survive_a_lorentz_boost` applies a genuine SO(1,9) element.
This is deliberate: if raised/lowered indices are assigned **per tensor copy**
rather than **per contracted edge**, an edge joining two same-placement vertices
contracts with $\delta$ instead of $\eta$. A pure rotation cannot see that error —
$\delta$ and $\eta$ agree on the spatial block — so the quantity looks invariant
and is not. A boost mixes the timelike direction and exposes it immediately.
Any candidate that changes value under a boost is not an invariant, whatever a
rank count says.

## Method

**Candidate generation is a graph enumeration, not a search.** A fully
contracted scalar from $n$ copies of a $p$-form is a perfect matching on the
$np$ index slots, and the induced multigraph has every vertex of valence $p$.
So the complete candidate list at order $n$ is:

> symmetric $n\times n$ non-negative integer matrices, zero diagonal, every row sum $= p$

Two exact pruning rules: no self-loops (a self-contraction of an antisymmetric
form vanishes), and $m_{ij} \le p-1$ for the self-dual case (multiplicity $p$
forces a factor $F\cdot F = 0$). Disconnected graphs are products of lower
invariants and can never raise the rank, so they are dropped.

**Functional independence is a Jacobian rank.** Given candidates $I_a$, the generic rank
of $\partial I_a/\partial A^k$ is the number of functionally independent
invariants. A uniformly random point over a large finite field attains that
generic rank with high probability. The order-12 certificate is recomputed at
four points under three primes to guard against unlucky specializations.

**Homogeneous polynomial independence uses stacked gradients.** If a linear
combination of same-degree homogeneous polynomials vanishes identically, its
gradient vanishes at every point. Therefore a rank-$r$ matrix obtained by
concatenating their exact gradients at several points proves at least $r$
linearly independent polynomials. The order-12 run attains the known upper
bound 72, so the lower and upper bounds coincide.

$I$ is multilinear in its $n$ vertices, so
$\partial I/\partial A^k = \sum_v [\text{graph with } F \text{ at } v \to P(e_k)]$.
The production code contracts the tensor network once and reverse-
differentiates its globally optimized contraction tree, reusing intermediates
for all vertices. A separate amputated-tensor implementation remains as a
correctness oracle in the tests.

**Arithmetic is exact over $\mathbb{F}_p$.** Invariant values span an enormous
dynamic range, so a float SVD forces you to guess a rank tolerance — and that
guess *is* the answer. Over $\mathbb{F}_p$, zero means zero. Every result is
confirmed under a second prime.

## Read this before changing the core

`np.einsum` multiplies **all** operands before summing. At order $n$ the
products reach $p^n$, int64 wraps past $2^{63}$, and **numpy raises nothing** —
you get plausible, confident, wrong numbers. Asking for an optimised path does
not save you: for a symmetric contraction numpy will return a single step
contracting all $n$ operands at once.

This bug appeared three times during development. It produced rank 9 and then
rank 7 for a quantity whose mathematical ceiling is 5. It was only ever caught
by comparing against exact bigint arithmetic.

`mod_einsum` therefore contracts **strictly two operands at a time** with
reduction after every step, plus an explicit overflow guard. Large pairwise
contractions use float64/BLAS only when every unreduced integer sum is proven
to stay below $2^{53}$, so the result remains exact. Do not replace this with a
bare multi-operand `np.einsum` call.

## Run it

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

pytest tests/ -v              # gates. must pass before trusting anything.
python3 scripts/run_6d.py     # reproduction. must print PASS.
python3 scripts/run_10d.py    # orders 4, 6, 8 under both primes
```

`run_6d.py` is the gate. If it does not print PASS, every 10D number is
meaningless — fix that first.

The exact catalog is committed. To regenerate it, install
[nauty](https://pallini.di.uniroma1.it/) so `geng` and `multig` are on `PATH`,
then run:

```bash
python3 scripts/generate_graph_catalog.py
```

The checkpointed degree-10 workflow uses nauty 2.9.3 and is documented with
copy-paste commands in [`docs/degree10.md`](docs/degree10.md). A clean checkout
can revalidate the committed formulas without regenerating the discovery
catalog:

```bash
python3 scripts/degree10_pipeline.py validate \
  --selection-result results/10d_order10.json \
  --skip-catalog-check \
  --primes 32749 32719 \
  --jacobian-seeds 20260729 20260730 20260731 \
  --value-seed-start 20260801 --value-samples 16 \
  --out /tmp/10d_order10_revalidated.json
```

The committed order-12 formulas can likewise be revalidated without any
generated catalog or checkpoint:

```bash
.venv/bin/python scripts/degree12_pipeline.py validate \
  --selection-result results/10d_order12.json \
  --primes 32749 32719 32693 \
  --seeds 20260729 20260730 20260731 20260732 \
  --out /tmp/10d_order12_revalidated.json
```

See [`docs/degree12.md`](docs/degree12.md) for exact shard generation and
discovery commands.

## Scope

The repository now gives an explicit trace-contraction realization of all
**81 functionally independent invariants through order 12**, along with all
72 homogeneous degree-12 directions. It does not claim that these 81
functions form a freely generated polynomial ring, nor does it claim an
exhaustive count of every order-12 contraction graph. Degree-12 completeness
uses the known Hilbert-series upper bound together with the exact rank-72
lower bound.

## Layout

```
src/sdinv/modp.py      exact F_p linear algebra, overflow-safe pairwise einsum
src/sdinv/forms.py     p-forms, Hodge dual via index complement, self-dual projector
src/sdinv/graphs.py    exact multigraph certificates, nauty generation, catalogs
src/sdinv/contract.py  contraction evaluation, optimized reverse-mode Jacobian
src/sdinv/catalog.py   atomic checksummed streaming graph shards
src/sdinv/checkpoint.py durable identity-checked rank checkpoints
src/sdinv/spinor_adapter.py future exact trace/spinor column-space comparison
src/sdinv/stress.py      paper-normalized M, N, T, ModMax I8/I12 identities
src/sdinv/invariant_registry.py committed low-degree and degree-12 import registry
src/sdinv/exactmap.py    modular solves, CRT, rational maps, complements
scripts/run_6d.py      reproduction gate
scripts/run_10d.py     two-prime complete computation through order 8
scripts/generate_graph_catalog.py  exact nauty catalog generation
scripts/degree10_pipeline.py generation, scheduling, discovery, validation, benchmarks
scripts/degree12_pipeline.py bounded order-12 shards, discovery, validation, benchmarks
scripts/stress_flow_pipeline.py exact stress map and physics-stage artifact
tests/test_core.py     core correctness gates
tests/test_degree12.py committed order-12 certificate gates
tests/test_stress_flow.py paper, Lorentz, map, and closure-pilot gates
```

## Next questions

1. The octic map is resolved for the explicitly selected six source expressions
   plus the quartic square; the full seven-dimensional map and product
   correction are in `results/order8_change_of_basis.json`. Other source
   selections must be specified and their span checked independently.
2. What are compact tensor-word expressions for the 62 saved order-12 graph
   contractions?
3. What are explicit differential relations for `I12_61` and `I12_62`
   against the 81-row functional basis?
4. At what degree do the first nonlinear relations among the published
   generator counts appear?
