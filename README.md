# Constructive invariants and nonlinear chiral four-form stress flows

A unified research project connecting explicit Lorentz invariants of a self-dual
five-form in ten dimensions to the interactions reachable by stress-tensor flows.

The scientific question is: **which analytic, derivative-free interactions can
be generated from the free chiral four-form theory, and which additional
invariant generators enlarge that reachable set?**

Read [SCIENCE_STATUS.md](SCIENCE_STATUS.md) for current verification levels and
[CLAIMS.md](CLAIMS.md) for the claim boundaries. The editable unified research
draft is [manuscript/unified/main.tex](manuscript/unified/main.tex). Its compiled
PDF and the computation evidence live under `reports/` in the release bundle.

## Install and verify

Python 3.10+ is required for the flow package. Fresh foundation verification is
exercised on Python 3.13; fresh character generation requires a C++17 compiler.

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'
make combined-verify
```

Verification levels are explicit:

| Command | Work performed |
|---|---|
| `make combined-verify` | Regression/adversarial suite, foundation saved arithmetic, independent exact orbit algebra, packaged degree-eight map arithmetic |
| `make combined-fresh` | Above plus fresh character tables, independent rank/orbit tensors, cross-engine bridge, degree-eight tensor map, disjoint degree-twelve physical holdouts |
| `make combined-full` | Above plus complete bounded degree-ten and degree-twelve physics-map calculations |
| `make combined-endpoint` | Explicit ten-stage rational endpoint and all 68 degree-twelve constraints |

The full degree-twelve calculation is expensive and saves each tensor evaluation.
A partial run is not a proof. Resume the same compatible checkpoint directory:

```sh
python scripts/resume_degree12_proof.py --workers 2
python scripts/prove_combined_degree12_map.py --verify-only
```

Use one or two workers on an 8 GB machine. Work budgets bound planned tensor
arrays; they do not enforce an OS-wide memory quota. Checkpoint namespaces bind
code, the invariant registry, source coordinates, points, bounds, and primes.
Never delete checkpoints to hide a failed result.

## Generate and classify interactions

```sh
python independent/verify.py --target '1,0,0,0,0,0,0,0,0,0' --output reports/reachable-example.json
python independent/verify.py --classify reports/reachable-example.json --output reports/reachable-classification.json
python independent/verify.py --classify examples/combined/unreachable-octic.json --output reports/unreachable-example.json
```

The endpoint contains all 96 exact interaction coefficients. Classification
returns either signed controls or nonzero obstruction residuals in the declared
finite model. Use exact invariant IDs and rational strings or integers.

## Scientific components

- `foundation/`: canonical scientific files plus a complete pinned source archive of `selfdual-5form-invariants`, with
  original attribution, graph formulas, character derivation, and standalone
  rank/orbit verifier. The upstream repository remains separate and unchanged.
- `data/combined/`: the actual complete 96-coordinate homogeneous registry through
  degree twelve, with a source lock and hashes. No ignored local run is needed.
- `src/chiral4form/`: tensor physics, projected derivatives, stress traces,
  polynomial vector fields, exact orbit and obstruction infrastructure.
- `independent/`: a separate SymPy implementation that imports neither research
  engine. It checks exact orbit algebra and constructs explicit endpoints.
- `scripts/`: fresh map proofs, cross-engine verification, resumable computation,
  and the combined verification entry point.
- `verification/`: preserved predecessor certificates and historical evidence.
- `reports/`: new logs, manifests, fresh certificates, and independent results.

## Central distinctions

The atlas has 81 algebraically independent invariants. The complete homogeneous
spaces through degrees 4, 6, 8, 10, 12 have dimensions 1, 2, 7, 14, 72, giving
96 coefficient slots. These are different counts. No global polynomial-ring
presentation or global real-orbit separation is claimed.

In the degree-eight reduced flow model, `Omega8=b2+16*a^3` is transported
homogeneously. Its free-seed orbit has dimension four although its linear hull
has dimension five. Six specified extra generators complete that truncated
model, with necessity restricted to the stated catalogue.

The independently checked reduced-model orbit dimensions are 4, 6, 10 through
cutoffs 8, 10, 12. Relating those exact model calculations to tensor physics is
handled by separate bounded map certificates. Signed-control reachability does
not establish causality, unitarity, stability, nonanalytic ModMax reachability,
Type-IIB completion, or an all-orders orbit classification.

The original entry points remain available in the [historical README](docs/history/pre-unified-readme.md).
Legacy degree-wise closure numbers remain imported claims. Historical success
labels are not independently verified merely by being imported. A linear invariant
hull is not a nonlinear reachable orbit. Legacy documents and certificate status labels describe their own computation
stages; they do not override the current unified science status.
