# Cluster jobs

**No cluster run has occurred.** No cluster was available for this work, and
nothing in the manuscript depends on one. These files are prepared so that a
human with access can launch them without re-deriving the setup, and so that the
claim "this needs a cluster" is backed by an actual job script rather than an
assertion.

All three scripts are deterministic: fixed commit, fixed dependency lock, fixed seed
per shard, immutable output directory, checkpoint interval, memory cap and time
cap. A shard that is re-run produces the same output or refuses to overwrite.

## Order of operations

1. Run `spinor_degree10_smoke.slurm` on one node. It must finish and reproduce
   the degree 4/6/8/10 ranks 1, 2, 7, 14 before anything larger is launched.
2. Then `spinor_degree10_nostop_production.slurm`. This is the computation the
   project attempted and did not finish: the degree-10 scan with the Hilbert
   stopping rule disabled, run to a terminal status. On a laptop it was killed
   after about four hours at rank 12 of 14, still inside its first batch.
3. Only then `spinor_degree12_production.slurm`.

Step 2 is the one that would change a claim. Degrees 4, 6 and 8 reach
`candidate_exhaustion`, so their ranks are derived by searching the ansatz out
rather than by stopping at a target supplied from outside. Degree 10 has no such
status, and the manuscript says "saturation" there instead of "exhaustion"
because of it. Either terminal outcome is a result: exhaustion at rank 14 shows
the port-graph ansatz suffices at degree 10, exhaustion below 14 shows it does
not and locates what is missing, as happened at degree 8. A time-limited stop is
not a result, and the recorder marks it as such rather than reporting the rank
reached.

## What the smoke test proves

That the environment on the cluster reproduces results already established on a
laptop. If it does not, the discrepancy is in the environment, and the production
job would have inherited it silently.
