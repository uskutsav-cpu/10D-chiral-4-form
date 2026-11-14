# Degree-10 no-stop: terminal status

Artifact: `results/degree10/no_stop_terminal.json`.
Claim definition and the decision not to run it:
`docs/DEGREE10_NO_STOP_CLAIM.md`.

## Terminal status by degree

| degree | final rank | new + product | unique graphs tried | stopped by |
|---:|---:|---|---:|---|
| 4 | 1 | 1 + 0 | 2 | **candidate_exhaustion** |
| 6 | 2 | 2 + 0 | 5 | **candidate_exhaustion** |
| 8 | 7 | 6 + 1 | 32 | **candidate_exhaustion** |
| 10 | — | — | — | **killed** |

Degrees 4, 6 and 8 terminated by searching the candidate grammar out, not by
stopping at the Hilbert target. That is what makes their ranks informative about
the ansatz rather than about a number supplied from outside.

Degree ten has **no terminal status**. The process was killed after roughly four
hours of CPU time, at rank 12 of the target 14, without leaving its first batch
of 64 candidates. That is not a result in either direction and is not reported
as one.

## Completeness taxonomy

The four things that could be meant, and which applies where:

| meaning | degrees |
|---|---|
| enumeration complete **by proof** | none |
| complete **within the declared candidate grammar** | 4, 6, 8 |
| search **exhausted under configured bounds** | none |
| search **merely failed** to find another direction | none |
| **no terminal status** | 10 |

Degrees 4/6/8 are the second row and nothing stronger: the port-graph generator
produced no further distinct candidate, which is a statement about that grammar,
not about the space of all invariants.

## What depends on this

Nothing in either manuscript depends on the degree-ten run. The single
no-stop-dependent statement is the degree-eight one, and that run completed. The
degree-ten agreement between the two descriptions comes from the common-sample
span comparison, which is exact over `F_p` and never consults the enumerator's
terminal status. The manuscript says **saturation** at degree ten, never
**exhaustion**.

## What is prepared

`cluster/spinor_degree10_nostop_production.slurm`, deterministic and pinned, with
its configuration hash recorded in the artifact. **No cluster run has occurred.**
One submission on any suitable node would settle it, and nothing in the paper
waits on the answer.
