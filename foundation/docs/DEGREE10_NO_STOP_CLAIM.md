# What the degree-10 no-stop run would certify — and what depends on it

Required before running anything: state the claim, then check whether the
manuscript makes it. This document does that, and the answer changes what should
be done.

## The claim such a run would support

Running the spinor enumeration with the Hilbert stopping rule disabled makes the
terminal status informative. Every archived scan stopped at the Hilbert target,
a number supplied from outside the computation, so none of them says anything
about the candidate family. With the rule disabled the terminal status
distinguishes:

| terminal status | what it certifies |
|---|---|
| `candidate_exhaustion` at rank 14 | the port-graph ansatz **suffices** at degree ten: the family was searched out and reached the full graded dimension |
| `candidate_exhaustion` below rank 14 | the ansatz **does not** suffice, and the shortfall is located, as happened at degree eight |
| `time_limit` / `killed` | nothing, in either direction |

The precise property is **exhaustion of the declared candidate grammar at degree
ten**, not enumeration completeness in any absolute sense. Even a successful run
would certify only that the port-graph generator produced no further distinct
candidate, within the configured bounds.

## What the manuscript actually claims

Searched both manuscripts for every no-stop-dependent statement. There is exactly
one, in the long form:

> Running the spinor enumeration with its stopping rule disabled, so that it
> terminates on candidate exhaustion rather than on an externally supplied
> target, reaches rank seven at degree eight --- but only when a separate family
> of structured tensor-word candidates is included alongside the port graphs.

That is a **degree-eight** claim, and the run supporting it **completed**:

    degree  4   rank 1   stopped_by candidate_exhaustion
    degree  6   rank 2   stopped_by candidate_exhaustion
    degree  8   rank 7   stopped_by candidate_exhaustion   <- the claim
    degree 10   ---      stopped_by killed

`verification/spinor_degree10_no_stop.json`.

**No claim in either manuscript depends on the degree-ten no-stop run.** The
degree-ten agreement between the two descriptions is established by the
common-sample span comparison, which is exact over `F_p` and does not use the
enumerator's terminal status at all. The manuscript already says "saturation" at
degree ten, never "exhaustion", precisely because that run did not finish.

## Decision

The run is **not performed**, and this is a decision rather than an omission.

- It certifies a property the manuscript does not assert.
- The previous attempt consumed roughly four hours of CPU and was killed at rank
  12 of 14, still inside its first batch of 64 candidates, so a local completion
  is not close.
- No authorised cluster exists for this work. Buying cluster time is out of scope
  and using an unauthorised account is not an option.

Running it would be doing work because a checklist asks for it, which is the
thing this phase explicitly forbids.

## What is prepared instead

`cluster/spinor_degree10_nostop_production.slurm` is a complete, deterministic
job — pinned commit, verified archive, immutable output directory that refuses to
overwrite, `--no-hilbert-stop`, greedy contraction ordering (the `auto` strategy
switches to exact dynamic programming around twenty operands and a degree-ten
amputated contraction has nineteen, which is what silently stalled two earlier
runs), and the dense-`I` plan that a 64 GB node can afford and a laptop cannot.
It hands the run to the recorder, which distinguishes `candidate_exhaustion` from
a resource stop rather than reporting a rank as if it were terminal.

**No cluster run has occurred.** Anyone with a node can settle this in one
submission; nothing in the paper waits on it.

## Wording the manuscript is held to

Permitted, and used:

> degrees four, six and eight terminate on candidate exhaustion
> at degree ten the enumeration saturates

Forbidden, and gated:

> the degree-ten enumeration is complete
> the ansatz is exhaustive at degree ten

The distinction between *enumeration complete by proof*, *complete within a
declared grammar*, *exhausted under configured bounds*, and *merely failed to
find another direction* is maintained: degrees 4/6/8 are the second of those,
degree ten is none of them.
