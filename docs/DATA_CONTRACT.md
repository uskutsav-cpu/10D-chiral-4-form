# Data contract with the predecessor atlas

The engine consumes graph/product definitions as data, not executable upstream
code. The target checkout is never substituted for the original atlas repository.

```bash
python -m chiral4form source --fetch --checkout .cache/upstream --output data/imported
python -m chiral4form plan --source .cache/upstream --config configs/degree12.json
```

The importer pins `3ed32805b38ce34216b34888f6539e3538e90fb9`, checks Git objects and
SHA-256 content, rejects dirty/wrong-commit existing checkouts, and exports an
explicit registry with ordered homogeneous bases through twelve. It never resets
an existing checkout. Network access is necessary for an absent source checkout.
The older wrapper `scripts/freeze_source_snapshot.py` delegates to this importer;
use `--help` for its current arguments rather than the old scaffold invocation.

The bundled fixture contains nine primitive graph formulas through degree eight
plus the quartic-square product and source attribution. It is intentionally not a
full predecessor snapshot. Full degree-ten/twelve computation requires the pinned
source import; missing formulas cause an error, never a fabricated slot.

Imported counts and source hashes only prove provenance. Reproduction requires a
new evaluator run; orbit classification additionally requires nonlinear control
analysis. `data/baseline/paper2_baseline.json` retains historical numbers with that
warning, while `verification/` contains this release's actual calculations.
