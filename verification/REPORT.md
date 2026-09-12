# Computational-engine release verification

**Release:** 0.2.0. **Base:** `0dedbeb`. **Local regression result:** **142 passed, 0 failed**.

This report describes computations actually executed, not planned experiments.
The repository contains 30 Python source modules and 2,636 source lines,
plus tests, configurations, documentation, and recorded results.

## Executed experiments

| Experiment | Primes | Samples per prime | Outcome |
|---|---:|---:|---|
| Degree-six tensor-to-flow map | 3 | 4 fit + 2 holdout | Basis ranks 1,2; independent analytic sextic map agrees |
| Degree-eight tensor-to-flow map | 3 | 9 fit + 2 holdout | Basis ranks 1,2,7; hull rank 5, Lie lower bound 4 |
| Genuine generalized degree-eight flow | 3 | 9 fit + 2 holdout | Hull and Lie ranks both 10 with six specified extras |
| Bounded-residual degree-eight map proof | 7 | 9 fixed ternary tensors | All 63 tensor evaluations / 1,323 cleared residual checks pass |
| Exact degree-eight orbit analysis | Rational arithmetic | Symbolic | Four-dimensional invariant graph plus explicit four-segment construction |
| Exact finite-catalogue completion | Rational arithmetic | Symbolic | Six necessary extras within the stated catalogue; old five-extra choice misses one constraint |
| Finite-catalogue subset search | GF(30011) | All 64 subsets | Only the full six-element set meets the requested rank-10 test |
| Rational reconstruction | Two fit primes + one unused prime | 20 nonzero coefficients | Candidate agrees; not itself used as the physics identity proof |

Tensor-fit primes: 30011, 30013, 30029. The bounded proof additionally uses
30047, 30059, 30071, 30089, with the same bounded integer tensors across fields.
The final map-certificate hash is
`53f82d39cbed38d089cc7fba124f7d7f0fca31228cc2d4db9bbea88ed1afe73f`.

## Scientific scope

For interaction coefficients `(a,j,k,b1,b2,b3,b4,b5,b6,p)`, the reduced analytic
model has `Omega8=b2+16*a^3`, transported by `Omega8'=60*u_tr1*Omega8`.
The free-seed orbit satisfies `k=b3=b4=b5=b6=0` and `b2=-16*a^3`.
Its dimension is four although its linear hull has dimension five. The proof
uses unrestricted signed polynomial controls in the declared degree-eight
truncation, not physical positivity/causality restrictions.

The physics-to-basis certificate is conditional on the external characteristic-
zero invariant-space upper bounds 1,2,7 and the documented HLS stress and derivative
conventions. Those inputs are explicitly cited and require expert review; they
were not independently derived by a character calculation here. No priority or
publication claim is certified by these files.

## Engineering checks

Editable installation succeeded. Python compilation, shell syntax, adversarial
certificate tests, array-budget checks, fit/holdout separation, actual small-
pipeline integration, checkpoint interruption/resume, and full `make verify`
passed. The GitHub Actions matrix is configured but was not run remotely.
The detailed environment and hashes are in `RELEASE_MANIFEST.json`.

## Not represented as finished

Full degree-ten/twelve tensor and orbit reproduction; all-orders classification;
alternative pure-stress ModMax reachability; global minimality over arbitrary
extra invariants; causality/unitarity; Type-IIB supersymmetric completion; mentor
approval. The larger-degree importer/evaluator infrastructure exists, but these
research claims remain open.

## Reproduce

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
make verify
bash scripts/run_research.sh core
```

`verify` rechecks packaged certificate arithmetic but not all tensor samples.
`core` recomputes the complete implemented degree-eight program with progress
and content-addressed resume. Neither command pushes to GitHub.
