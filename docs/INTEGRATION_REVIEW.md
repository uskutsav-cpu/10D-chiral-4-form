# Review guide for the unified project

The integration retains the invariant and flow engines as separately inspectable
implementations. It adds a complete pinned registry, an independent exact orbit
checker, exact generalized-completion blocks, explicit endpoint generation,
and self-contained fresh tensor-map proof commands.

Review in this order:

1. `data/combined/source-lock.json`: source commits, tracked file hashes,
   symlink target and complete registry binding.
2. `independent/verify.py`: independent reconstruction of degree-ten tangency,
   all 36 degree-twelve source rows, 68 constraints, vertical block and ten-stage
   endpoint algorithm. Its symbolic triangular conditions establish the
   construction for arbitrary targets, beyond the example endpoints.
3. `independent/completion.py`: exact leading blocks for all selected generalized
   generators. Minimum arbitrary-scalar/nonlinear-control cardinality is not
   asserted.
4. `scripts/prove_combined_degree12_map.py`: explicit derivative denominator
   clearing, absolute bounds, fresh tensor points, bounded CRT and injectivity.
   It does not depend on the ignored historical interpolation runs.
5. `scripts/verify_combined.py`: quick/fresh/full command scope and stage logs.
6. `docs/PHYSICS_CONVENTION_AUDIT.md`: auxiliary versus physical fields,
   quartic normalization and the trace sign fixed by HLS (2.33).
7. `manuscript/unified/main.tex`, `SCIENCE_STATUS.md`, and `CLAIMS.md`: mathematical
   narrative, exact evidence levels and unresolved research boundaries.

The source-map proof workers are checkpointed and may take substantial time.
Use the supervisor to resume compatible files and expose failures immediately.
A complete set of checkpoints still requires the final arithmetic verifier.

The unified manuscript is an editable research draft. The original graph atlas,
character derivation and source manuscripts remain available in `foundation/`.
No automatic merge, journal submission or external outreach is performed.

The inherited absolute-bound metadata retains its historical integrality-description
field for provenance. The new proof does not require that mathematical premise:
it multiplies each target by an explicit stress/derivative denominator clearer
before applying bounded integer CRT. Compatibility checks on old certificate
labels are not promoted to a proof of coefficient integrality.

The complete map has two replay paths: strict production input validation and
`independent/map12_arithmetic.py`, which imports neither tensor engine. The latter
checks all saved residues, CRT lifts, denominator scaling and source identities,
while retaining the declared absolute-bound derivation as an assumption.
