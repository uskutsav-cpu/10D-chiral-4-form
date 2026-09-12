# Reproducibility Policy

## Provenance

Every generated research artifact must record:
- source commit SHA(s);
- basis/invariant-registry fingerprint;
- conventions identifier;
- prime(s);
- random seed(s);
- sample count;
- exact script/CLI version;
- hashes of imported input artifacts.

## Finite-field policy

A rank seen at one prime or one random point is never enough for a headline equality.

Use two separate roles:
1. **reproduction primes** matching historical artifacts;
2. **fresh audit primes/seeds** chosen independently from development runs.

When rational reconstruction is required, record:
- CRT modulus;
- uniqueness/height bound;
- reconstructed fractions;
- at least one held-out prime not used for fitting.

## Baseline historical fields

The predecessor interacting-flow assembly reports fit primes
`32749, 32719, 32693, 32771, 32713`
and holdout `32717`.

These should be preserved as historical reproduction inputs, not used as the only evidence in this repository.

## Fresh audit inputs

Large-prime audit runs may use fields such as
`998244353`, `1000000007`, `1000000009`, or other backend-supported primes **only after checking every fixed denominator is invertible and arithmetic remains exact in the implementation**.

Do not blindly copy these values into a backend with a restricted integer range.

## Holdout discipline

Maintain at least one seed/prime set that is not used while developing a claimed identity. Open it only after the derivation/code has been frozen.

## Negative controls

Every central pipeline should include at least one mutation that must fail, for example:
- remove a required generalized generator;
- inject an invalid quotient direction;
- flip a convention sign;
- perturb a claimed rational coefficient;
- replace a self-dual sample with an unprojected control where appropriate.

## No floating-rank claims

Floating-point SVD may be used for debugging or visualization, but not for certificate-grade rank claims where exact modular/rational arithmetic is available.
