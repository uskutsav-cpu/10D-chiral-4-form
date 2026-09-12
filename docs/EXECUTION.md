# Execution guide

## Ordinary laptop workflow

```bash
source .venv/bin/activate
make verify
make plan
make smoke
make degree8
make completion8
make map-proof
```

Use `bash scripts/run_research.sh core` for the same workflow plus rational
lifting and the symbolic reports. The script stops on the first failure.
A repeated command reuses only matching, hash-valid checkpoints. No repeated
code submissions are needed while a foreground command is running.

Set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, and
`VECLIB_MAXIMUM_THREADS=1` to avoid oversubscription. The supplied script does
this by default. Computation uses CPU tensor contractions; it does not need
CUDA, MATLAB, a GPU, or a new large-RAM server for the demonstrated degree-eight
experiments. Higher-degree memory requirements are determined by `plan`.

## Source freeze and larger jobs

```bash
python -m chiral4form source --fetch --checkout .cache/upstream --output data/imported
python -m chiral4form plan --source .cache/upstream --config configs/degree12.json
```

Review every `within_budget` result before launching. The default degree-12
configuration is a conservative starting policy, not a claim it fits every
12-vertex graph. Raise budgets only after checking actual available memory.
The runner refuses an infeasible plan before expensive tensor work. There
is no automatic inference of completion from imported rank tables.

```bash
python -m chiral4form run --source .cache/upstream \
  --config configs/degree12.json --output runs/degree12
```

Alternatively use the same source JSON exporter with your own execution
backend. The polynomial analysis accepts explicit, schema-checked models
without tying them to a particular tensor implementation.

## Cluster workers

Each prime is independently checkpointed:

```bash
python -m chiral4form worker --prime-index 0 --config configs/degree8.json \
  --output runs/cluster-degree8
```

The optional `scripts/cluster.slurm` runs a three-prime array. It does not
invent an account, partition, or filesystem for your institution. Set those
through your cluster's normal `sbatch` options. An index outside the config's
prime list fails. Sharded models can be lifted by passing their three model
JSON paths explicitly to `chiral4form lift`.

## Interruption and recovery

Ctrl-C records an interrupted stage and removes its lock. Rerun the same
command to resume. A hard kill or power loss can leave a `.lock` file. Read
its PID and verify that the process is dead **on the machine that created
it** before removing that one lock. Do not delete an entire run directory or
its successful checkpoints to recover one interrupted stage.

Use `python -m chiral4form status --output runs/degree8` to inspect state.
The bounded-map proof has `summary.json` on completion and content-addressed
sample receipts during execution. `verify-map8 --reevaluate-tensors` bypasses
stored sample values for a full numerical recheck.

## Independent verification and publication

A successful `make verify` validates infrastructure and packaged certificates.
It is not equivalent to rerunning a degree-twelve tensor experiment. The
release ledger states which experiments were actually run. Keep source
commits, exact configurations, formula/basis fingerprints, input tensors,
primes, fit/holdout separation, and failure logs with any paper artifact.
