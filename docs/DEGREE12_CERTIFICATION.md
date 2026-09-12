# Degree-12 characteristic-zero certification layer

This layer is intentionally separate from the already-frozen modular result.
It never upgrades the degree-12 orbit claim to a characteristic-zero physics
theorem merely because modular ranks stabilize.

## Pipeline

1. **Normalization audit** across all modular models.
2. **Modular obstruction rebuild** in the 81-dimensional weight-100 ansatz.
3. **Pivot/rank stabilization** of the 13-dimensional equation space and the
   10-dimensional free-seed control space.
4. **Asymmetric bounded reconstruction** using separate numerator and
   denominator windows.
5. **Fault-tolerant fallback** with leave-k-out voting and a Gaussian lattice
   candidate when a small number of modular images may be unlucky.
6. **Subspace-level lift** of canonical row spaces, followed by exact modular
   reduction checks.
7. **QQ linear-algebra gate**: equation rank 13, obstruction nullity 68,
   degree-12 leading rank 68, control rank 10, and the 10-dimensional
   sandwich.
8. **Direct physics gate remains separate.** A successful modular/QQ lift is
   recorded as an algebraic candidate, not automatically as a direct
   characteristic-zero tensor proof.

## CLI

```bash
python -m chiral4form degree12-certify \
  --models runs/degree12-reduced \
  --holdout-prime 30293 \
  --max-bad-primes 1 \
  --output verification/degree12/qq/certificate.json
```

The command is fail-closed. If any normalization, rank, pivot, reconstruction,
or verification gate fails, `qq_theorem` remains false.

## Proof-bound helper

`modular_reconstruction.sage_style_echelon_height_bound()` implements the
sufficient height inequality used by multimodular echelon algorithms after a
valid characteristic-zero source-matrix height bound is supplied. This helper
is not invoked automatically because the current degree-12 physics pipeline
has not yet produced a rigorous global source-height bound.
