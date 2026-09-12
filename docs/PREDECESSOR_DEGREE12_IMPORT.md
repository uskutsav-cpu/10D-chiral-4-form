# Importing predecessor degree-12 modular certificates

The pinned predecessor commit already contains six exact interacting-flow
coefficient certificates at primes 32693, 32713, 32717, 32719, 32749 and
32771.  They use the same `tau=48*T` normalization and include independent
holdout validation.

This adapter converts those certificates into the current polynomial-vector-
field JSON schema, maps the predecessor generator IDs to the current generator
IDs, restricts the resulting fields to the frozen six-parameter degree-10
orbit, and combines them with the current nine reduced prime models.

Importing a predecessor certificate is **not** a characteristic-zero proof.
The normal `degree12-certify` fail-closed gate must still succeed.
