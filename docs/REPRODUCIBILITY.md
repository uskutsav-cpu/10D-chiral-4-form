# Reproducibility and proof policy

## Distinct evidence levels

- A nonzero exact modular minor certifies a characteristic-zero **lower** bound
  for an integral/rational matrix with invertible denominators.
- A matching analytic upper bound can turn that lower bound into an equality.
- Agreement across finite samples/primes alone is not a polynomial identity.
- Bounded CRT reconstruction with a held-out prime produces a validated
  **candidate** unless a bound or independent identity argument closes the gap.
- `map_proof.py` instead bounds cleared integer residuals at fixed bounded
  tensors and uses a modulus larger than twice the bound. The external Hilbert
  upper bounds make the evaluation map injective; its exactness argument is in
  `DEGREE8_MAP_CERTIFICATE.md`.
- Linear invariant hull rank, Lie-rank lower bound, and nonlinear orbit dimension
  are different quantities. Only the explicit degree-eight construction closes
  both orbit bounds here.

## Actual released inputs

The tensor runs use primes 30011, 30013, 30029, fit and holdout seeds specified
in each configuration, and exact modular arithmetic. The bounded-residual proof
uses seven primes and nine fixed ternary electric tensors (63 tensor evaluations),
with all electric vectors stored in the certificate.

The NumPy tensor backend accepts primes at most 65521 and rejects incompatible
memory/work requests before executing contractions. Do **not** send arbitrary
large primes to it. The independent Python-integer linear algebra supports a
larger validated domain; this does not widen the tensor backend's safe range.

The optional float64 BLAS fast path is used only when all integer products and
sums are below its exact-integer limit. It does not compute floating-rank/SVD
certificates. Tests compare it with the integer path.

## Provenance and resume

Each run records configuration, basis fingerprint, engine hash, environment,
primes, sample counts, and reports. Source import pins a Git commit and hashes
actual files. The bundled low-degree registry is attributed as an excerpt.
Checkpoint keys include engine/configuration/registry state. Atomic writes and
exclusive locks prevent incomplete success records. Changed data/code invalidates
resume keys. Never delete a lock before verifying its process is gone.

`verification/` contains only results actually computed for this release. Runtime
outputs go under ignored `runs/`; they are not silently committed. The configured
GitHub CI matrix has not been executed remotely by this packaging session.
