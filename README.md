# Stress-Flow Orbits and Obstructions for 10D Chiral 4-Forms

This repository is the **Paper 2** research program built on the explicit invariant classification of a self-dual five-form in ten dimensions.

## Scientific question

For a nonlinear ten-dimensional chiral four-form theory with self-dual five-form field strength \(F_5\), determine which interactions are reachable from a chosen seed by scalar stress-tensor flows

\[
\partial_\lambda \mathcal V = f\!\left(\operatorname{Tr}T,\operatorname{Tr}T^2,\ldots,\operatorname{Tr}T^{10};\lambda\right),
\]

identify **intrinsic obstruction classes** to reachability, and determine the **minimal genuine enlargement of the flow-generator algebra** needed to remove those obstructions.

The flagship physics endpoint is a convention-controlled test of the ten-dimensional ModMax-like/conformal sector.

## Why a separate repository?

The predecessor repository, [`uskutsav-cpu/selfdual-5form-invariants`](https://github.com/uskutsav-cpu/selfdual-5form-invariants), answers primarily:

> **What Lorentz invariants exist?**

This repository answers:

> **Which nonlinear theories can stress-tensor flows actually generate?**

That separation is deliberate. It keeps the invariant atlas immutable while allowing the flow problem to develop its own definitions, certificates, conjectures, and manuscript.

## Imported baseline

The current source project reports the following exact/provisional Paper-2 baseline through five-form degree 12:

| degree | full invariant space \(A_d\) | static stress span | free-seed dynamic closure \(R_d\) | quotient \(Q_d=A_d/R_d\) |
|---:|---:|---:|---:|---:|
| 4 | 1 | 1 | 1 | 0 |
| 6 | 2 | 1 | 1 | 1 |
| 8 | 7 | 2 | 3 | 4 |
| 10 | 14 | 2 | 11 | 3 |
| 12 | 72 | 4 | 67 | 5 |

These numbers are **imported claims until independently frozen and revalidated in this repository**. The machine-readable record is `data/baseline/paper2_baseline.json`.

The conceptual surprise is that nonlinear closure is dramatically larger than the static stress algebra: at degree 12 it reaches **67/72** invariant directions despite a static stress span of only **4**.

## First theorem target

The first intrinsic obstruction is a sextic direction \(K_6\). In the source calculation its quotient coordinate obeys

\[
\dot q_6 = 40\,a(\lambda)q_6,
\]

so \(q_6=0\) is flow invariant. The target here is to replace computational evidence with a short analytic proof and a convention-independent tensor/representation statement.

Crucially, the intended statement is:

> **transported, never created from a seed with \(q_6=0\)**,

not “\(K_6\) vanishes in every pure stress flow.”

## Research gates

A strong core paper does **not** depend on proving an all-orders theorem. Submission-ready core:

1. exact nonlinear reachability through degree 12;
2. analytic intrinsic sextic obstruction theorem;
3. basis-independent obstruction/annihilator certificates for \(Q_8,Q_{10},Q_{12}\);
4. genuine generalized-generator analysis (not seed augmentation);
5. conformal/ModMax reachability test with an explicit certificate.

The stretch result is an all-orders finite obstruction module / reachability criterion.

## Repository map

- `ROADMAP.md` — phase-by-phase research program and kill tests.
- `CLAIMS.md` — claim ledger; no headline result is “established” without a certificate.
- `THEORY.md` — precise definitions and distinctions used by the code.
- `docs/LITERATURE.md` — closest prior work and novelty boundaries.
- `docs/REPRODUCIBILITY.md` — finite-field, seed, holdout, and provenance policy.
- `docs/DATA_CONTRACT.md` — interface to the invariant-atlas repository.
- `docs/PAPER_PLAN.md` — manuscript architecture and publication gates.
- `src/chiral4form/` — exact linear-algebra and certificate infrastructure.
- `scripts/` — executable validation/import/obstruction tools.
- `data/baseline/` — frozen baseline metadata.
- `tests/` — falsification-oriented tests.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'
pytest -q
python scripts/validate_baseline.py
```

The baseline validator intentionally checks only internal consistency. It does **not** pretend to recompute the physics calculation imported from the predecessor repository.

## Research integrity rule

A result may be called **established** here only if:

- its scope is stated precisely;
- a reproducible certificate exists;
- at least one independent falsification path exists where feasible;
- modular evidence is not silently promoted to an unsupported characteristic-zero equality;
- basis-dependent complements are not described as intrinsic obstructions;
- seed augmentation is not described as generalized-generator completion.

See `CLAIMS.md` for the live status ledger.
