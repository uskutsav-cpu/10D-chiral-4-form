# Exact resume commands

Every script writes incrementally and skips work already recorded, so an
interrupted run is resumed by re-issuing the identical command.

```bash
# smoke test (one node)
REPO_URL=<url> REPO_COMMIT=<sha> ARCHIVE_DIR=<path> \
  sbatch cluster/spinor_degree10_smoke.slurm

# production array; shards that already wrote COMPLETE exit immediately
REPO_URL=<url> REPO_COMMIT=<sha> ARCHIVE_DIR=<path> \
  sbatch cluster/spinor_degree12_production.slurm

# resume a single failed shard
REPO_URL=<url> REPO_COMMIT=<sha> ARCHIVE_DIR=<path> \
  sbatch --array=17 cluster/spinor_degree12_production.slurm
```

Local resume, no cluster:

```bash
python spinor_trace_bridge/scripts/run_comparison.py          # skips finished degrees
python spinor_trace_bridge/scripts/run_exact_jacobian.py      # skips finished (prime, seed)
python spinor_trace_bridge/scripts/run_archive_jacobian_exact.py --archive PATH
```
