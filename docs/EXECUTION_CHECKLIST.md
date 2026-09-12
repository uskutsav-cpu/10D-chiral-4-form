# Execution checklist for the implemented engine

This supersedes the original foundation-only checklist. Do not promote the
imported degree-ten/twelve linear closure counts into orbit dimensions.

## Run the delivered degree-eight program

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
make verify
bash scripts/run_research.sh core
```

`verify` checks tests and packaged certificate arithmetic. `core` additionally
recomputes real tensors, fits, holdouts, generalized generators, the bounded map
certificate, rational candidate lifting, and symbolic orbit/completion results.
Check `runs/*/summary.json` and the output status; a successful fit is not an
all-orders theorem. No Git write or remote push is performed by these commands.

## Next research milestones

1. Have the mentor review the HLS normalization, constrained derivative, and
   identification of the auxiliary-interaction flow with the intended physical
   flow. Review `THEORY.md` and `docs/DEGREE8_MAP_CERTIFICATE.md`.
2. Independently verify the explicit four-segment orbit construction and the
   nonlinear obstruction `b2 + 16*a^3`. Test genuine generator additions, not
   altered seeds. The earlier five-extra claim is refuted in this model.
3. Import the pinned full atlas using `python -m chiral4form source --fetch
   --checkout .cache/upstream --output data/imported`, then inspect the degree-12
   preflight before choosing a compute budget. A source import is not a rerun.
4. Extend the full coupling-space vector fields and orbit constraints through
   degrees ten/twelve. Do not substitute coefficient-support activation for a
   nonlinear orbit computation.
5. Treat conformal/ModMax reachability and all-orders statements as separate
   unresolved research. The delivered localization algebra only reproduces the
   published ModMax stress-square identity.
6. Freeze paper claims only after source/convention review, independent result
   reproduction, and a fresh novelty review. Code completeness is not paper
   completeness.
