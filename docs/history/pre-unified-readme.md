# 10D chiral four-form: stress-flow research engine

Executable research infrastructure for nonlinear chiral four-forms in ten
Lorentzian dimensions. It combines an independent five-form evaluator,
constrained derivatives, the interacting stress tensor, exact polynomial
control fields, and proof-scoped orbit/obstruction tools.

**Important correction:** a linear invariant hull is not a nonlinear reachable
orbit. The predecessor's degree-wise closure numbers remain **imported claims**;
degree-ten/twelve orbit dimensions are **not independently** established here.

## Start on your Mac

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
make verify
make smoke
```

`make verify` runs the regression suite, compilation checks, baseline schema
checks, and arithmetic verification of the packaged degree-eight certificate.
It does not rerun every expensive tensor computation. `make smoke` actually
creates fresh 10D self-dual tensors and fits the degree-six stress generators
at three primes, with separate fit and holdout points.

For the complete implemented degree-eight experiment:

```bash
bash scripts/run_research.sh core
```

This runs the tests, degree-six/degree-eight tensor fits, genuine generalized
flow fits, bounded-residual map certificate, rational lifting, and the exact
finite-model orbit/completion analyses. It executes in the foreground, prints
progress, and resumes hash-valid checkpoints on rerun. The commands do not
commit, push, delete, or overwrite another source checkout.

## What is implemented

| Layer | Executable capability | Boundary |
|---|---|---|
| Exact algebra | Prime validation, RREF, nullspaces, determinants, rank witnesses, rational reconstruction | A sampled rank is not automatically a characteristic-zero identity |
| Tensor engine | 126 independent components, Lorentzian Hodge star, boosts, contraction graphs, reverse derivatives | Array/work budgets are estimates, not an operating-system RSS quota |
| Physics | General HLS equation (2.33), anti-self-dual derivative projection, `tau=48T`, trace products | Analytic, derivative-free, classical interaction class |
| Degree-wise maps | Explicit graph/product bases, multiple primes, disjoint holdouts, preserved coupling monomials | Finite fits are labelled as fits |
| Flow geometry | Minimal linear invariant hulls, finite-depth Lie brackets, exact autonomous time jets | Hull rank and Lie lower bounds have different meanings |
| Obstructions | Exact polynomial ideal tangency and homogeneous basis covariance | Ideals live in **coupling space**, not directly in field-component space |
| Generalized flows | All eligible mixed monomials of `f(tau,S)` and finite-catalogue searches | Adding `S` is not changing the seed |
| Degree-eight result | Explicit nonlinear obstruction, constructive four-dimensional orbit, six-extra catalogue completion | Conditional on the certified map's stated physics/Hilbert inputs |
| Conformal tools | Exact quadratic extension `r^2=I4` and published ModMax stress-square reproduction | Alternative pure-stress ModMax reachability is unresolved |
| Larger calculations | Automatic import of all predecessor degree-10/12 graph formulas, preflight, checkpointing, one-prime workers | Full degree-10/12 recomputation is not shipped as a finished result |

## The degree-eight finding

Write the interaction coefficients as

```
a, j, k, b1, b2, b3, b4, b5, b6, p
```

multiplying `I4_1, I6_1, I6_2, I8_1,...,I8_6,I4_1^2` respectively.
For the declared polynomial stress-flow class, the exact reduced model has

```
Omega8 = b2 + 16*a^3
Omega8_dot = 60*u_tr1(lambda)*Omega8
```

Its free-seed orbit is the smooth four-dimensional graph

```
k = b3 = b4 = b5 = b6 = 0,   b2 = -16*a^3.
```

Four explicit control segments reach any point on that graph. Its **linear
hull has dimension five**, which is why the linear rank alone overcounts the
orbit. The six extras `I6_2,I8_2,I8_3,I8_4,I8_5,I8_6` complete the ten-dimensional
truncated coefficient model; each is necessary **within that catalogue**.
This is not global minimality over every possible choice of additional scalars.

Read [the derivation](docs/DEGREE8_RESULT.md) and
[the map certificate argument](docs/DEGREE8_MAP_CERTIFICATE.md). Novelty and
physical conventions still require mentor/literature review. No all-orders,
Type-IIB completion, causality, or nonanalytic classification is claimed.

## Commands

```bash
python -m chiral4form plan --config configs/degree8.json
python -m chiral4form run --config configs/degree8.json --output runs/degree8
python -m chiral4form status --output runs/degree8
python -m chiral4form degree8 --output runs/degree8-orbit.json
python -m chiral4form completion8 --output runs/degree8-completion.json
python -m chiral4form prove-map8 --output runs/map8-proof
python -m chiral4form verify-map8 verification/map8/map_certificate.json
python -m chiral4form modmax --output runs/modmax.json
```

## Import the real degree-twelve atlas

```bash
python -m chiral4form source --fetch --checkout .cache/upstream --output data/imported
python -m chiral4form plan --source .cache/upstream --config configs/degree12.json
```

The source command pins commit
`3ed32805b38ce34216b34888f6539e3538e90fb9`, checks Git blobs and SHA-256 hashes,
and translates actual JSON formulas without importing upstream Python.
An existing dirty or wrong-commit checkout is refused, not reset. The bundled
fixture is an explicitly attributed degree-eight excerpt, not a fabricated
complete upstream snapshot. The larger job only runs when the real registry
and resource preflight pass. See [execution instructions](docs/EXECUTION.md).

## Where to look

- `src/chiral4form/`: implementations, not generated result assertions.
- `configs/`: explicit degree/prime/holdout/resource policies.
- `verification/`: committed outputs of computations actually run for this release.
- `CLAIMS.md`: imported, conditional, refuted, and unresolved statements.
- `THEORY.md`: conventions and the precise flow class.
- `docs/`: proof arguments, architecture, reproduction, and mentor review.
- `manuscript/`: scope-controlled outline and a results draft, not a finished paper.

Upstream formulas and published mathematics are attributed in
[NOTICE](NOTICE) and [the source notes](docs/LITERATURE.md).
