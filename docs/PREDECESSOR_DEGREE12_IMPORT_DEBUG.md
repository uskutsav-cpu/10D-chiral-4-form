# Degree-12 predecessor certificate importer — debug hardening

The importer now treats the existing native reduced models as the canonical
schema.  It no longer invents optional `generator_catalogue` or `restriction`
metadata for predecessor-derived models.

Before a predecessor prime is accepted, the importer also verifies the complete
degree-4/6/8/10 vector-field sector against the frozen characteristic-zero
degree-10 model after exact reduction modulo that prime.  Degree-12-only
generators are required to act trivially on all first-24 lower coordinates.

This separates three questions:

1. Did the predecessor certificate pass its own fit/holdout checks?
2. Is it normalized identically to the current degree-10 model?
3. Does its reduced degree-12 polynomial support match the native reduced
   degree-12 model?

Only if all three pass is the modular image added to the extended certificate
set.  This is deliberately stricter than simply matching metadata strings.
