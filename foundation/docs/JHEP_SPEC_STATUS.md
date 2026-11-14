# Status against the JHEP execution specification

Specification: `JHEP_Full_Completion_Execution_Prompt.pdf`.
Last updated 2026-08-01, revised after the exact-Jacobian work landed.

## Workstream status

| WS | subject | status |
|---|---|---|
| **A** | exact trace-side degree-10 result | **COMPLETE** |
| **B** | product / primitive / reachable decomposition with incidence certificates | **COMPLETE** |
| **C** | trace reproducibility and clean clones | **COMPLETE** for clean-clone reproduction; additional primes and rational reconstruction NOT done |
| **D** | spinor baseline audit | **COMPLETE**, but re-done: the commit named in the specification did not exist and the archive was re-imported |
| **E** | rank-81 numerical certification | **COMPLETE by a different and stronger route**; the seed x scale x step float64 matrix was NOT run |
| **F** | degree-10 spinor no-Hilbert-stop | **PARTIAL** — an independent exact modular enumeration reaches rank 14 at degree 10 by saturation; the archive's own `--no-hilbert-stop` run did not complete locally |
| **G** | signature / complexification argument | **COMPLETE**, and it overturned the recorded blocker |
| **H** | gamma bridge | **COMPLETE**, exact over `F_p`, 86 tests |
| **I** | spinor–trace common-sample comparison | **COMPLETE** |
| **J** | physics interpretation | **DRAFTED**; requires coauthor confirmation |
| **K** | manuscript, figures, tables, bibliography, package, mock review | **COMPLETE** |

## G — the recorded blocker was wrong

The earlier record held that the oscillator frame's real form is Euclidean
`SO(10)`, where `*^2 = -1` admits no real self-dual five-form, and treated this
as blocking H, I and the manuscript.

Computing the metric from the archive's own operators gives real signature
`(5,5)` — **split**. There `*^2 = +1` on five-forms exactly as in Lorentzian, so
real self-dual five-forms exist and no such obstruction arises. A null frame
cannot be Euclidean at all, since Euclidean signature has no isotropic vectors.

What survives is precise and weaker: `(5,5)` and `(1,9)` are inequivalent real
forms, so the transition is complex over `R`; but both metrics have discriminant
`-1` up to squares, so over `C` and over `F_p` they are congruent. The bridge
constructs that congruence explicitly and checks it.

## H — the bridge, at both primes

    ker(forward) = the anti-self-dual 126     (span equality, not dimension)
    im(forward)  = the gamma-traceless 126    (span equality, not dimension)
    inverse . forward = the self-dual projector, exactly
    equivariant under GL(5) with character det(A)
    equivariant under Clifford reflections, which generate the full group
    no floating point and no tolerance anywhere in the package

## I — common-sample comparison, exact over `F_p`

    degree  4:  trace 1  spinor 1   spans equal, holdout validated
    degree  6:  trace 2  spinor 2   spans equal, holdout validated
    degree  8:  trace 7  spinor 6   spans NOT equal, containment strict
    degree 10:  trace 14 spinor 14  spans equal, holdout validated

The degree-8 shortfall is a property of the port-graph candidate family, not of
the bridge, and it independently reproduces the archive's own finding that
structured tensor-word candidates are needed there.

## E — what was done instead of the specified matrix

The specification asks for a float64 matrix over five seeds, three scales and
three step sizes. One configuration costs over ten minutes even with parallel
evaluation, so ninety of them was not feasible here. Two things were done
instead, and the second is strictly stronger than what was asked for:

1. Both archived Jacobians were re-analysed from the stored matrices under an
   explicit noise-floor rule. They separate cleanly: one with no rows at the
   floor, rank 81 stable across six orders of magnitude of tolerance and a gap of
   `2.0e7`; one with 48 of 83 rows at the floor, no gap, honest rank 35 — and 83
   if the normalisation rule is broken.
2. The finite-difference Jacobian was **replaced** by an exact analytic one,
   computed by amputation over `F_p`, validated by Euler's identity. Because the
   gamma-traceless basis is integral, the modular rank is a *rigorous* lower
   bound on the characteristic-zero rank rather than a probabilistic one.

An intermediate revision of this file recorded exact modular rank **59** for the
port-graph subset only, with 13 structured tensor-word candidates not
re-implemented. That is superseded. All 83 candidates are now implemented in the
same exact arithmetic and the certificate covers the complete selection:

    scheduled 83, evaluated 83, evaluation errors 0, zero rows 0
    Jacobian 83 x 126, exact modular rank 81
    cumulative rank by degree  4:1  6:3  8:9  10:21  12:81
    explicit 81 x 81 minor, determinant nonzero, two independent routines agree
    => rank over Q >= 81, unconditionally

The matching upper bound `126 - 45 = 81` remains analytic and from the
literature. The computation supplies the lower half only, at finitely many
points, and is never described as proving the generic rank.

## F — what is and is not established

Established: an independent exact modular enumeration reaches evaluation rank 14
at degree 10 by rank saturation, and its span equals the tensor-side span with
holdout validation.

Not established: a terminal `candidate_exhaustion` status from the archive's own
`--no-hilbert-stop` scan at degree 10. The run did not complete on this machine.
Cluster job files are prepared in `cluster/` but **no cluster run has occurred**.

## Remaining non-human work

| item | why not done |
|---|---|
| characteristic-zero check of `dim(B10 cap P10)` | **`dim D10` is now SETTLED** — `dim_Q D10 = 11`, `dim_Q Q10 = 3`, by exact rational closure with a held-out-prime-validated lift and an explicit non-vanishing minor. `B10 cap P10` remains a bound: three of the twelve published candidates have coordinates of too large a height to lift at the primes available. They are `P10_10`, `P10_11`, `P10_12` — which are exactly the three AMB-02 candidates. The coincidence is worth recording but is not established as causal: AMB-02 is a bracket-association ambiguity, while the lift failure is a height problem, and the plausible common cause is that the nested five-`N1050` structures are the most deeply nested expressions and so carry the largest-height rational coordinates |
| certified rational reconstruction | no claim depends on it; every modular result is used as a lower bound only, and the wording gate blocks "exact over `Q`" |
| degree-10 no-stop terminal status | run did not complete locally; `cluster/spinor_degree10_nostop_production.slurm` is prepared, no cluster run has occurred |
| incidence certificate at the four further primes | the generator exists and reproduces both stored primes exactly; ~15 min per prime, not completed here. Guards against a one-off slip; by the PO-09 analysis it discharges nothing |
| PO-03, PO-05, PO-07, PO-09 | open mathematical debts recorded in `docs/PROOF_OBLIGATIONS.md`; no claim is strengthened while its obligation is open. PO-08's cardinality half is now **discharged analytically** |

Done since the previous revision: the structured tensor-word candidates were
implemented in exact arithmetic, closing the 83-candidate gap; a systematic
primary-literature sweep was completed and is recorded in
`audit/RELATED_WORK_COMPLETE.md`. Novelty rows remain `PROVISIONAL` because
confirmation is a coauthor decision, not because the search is outstanding.

Also since then: the rank-81 certificate gained a second sample point; the
cardinality half of PO-08 was proved; PO-09 was narrowed to the two quantities
actually exposed; `SPINOR_JACOBIAN_STABILITY.md` was written; the degree-10
no-stop cluster job was written; and the incidence certificates gained a
versioned generator, having been produced inline.

## A defect class this revision kept finding

Six separate artifacts asserted a number that a certificate beside them
contradicted: two status documents on the Jacobian rank, the dimension
dictionary on the cumulative ranks, the manuscript on the same, the reproduction
records on test and gate counts, and the arXiv comments field on page and figure
counts. None was in the manuscript's *scientific* claims, which the existing
gates covered; all were in prose the gates did not reach.

The pattern in every case was a number maintained by hand next to an artifact
that already knew it. The response was the same each time: generate it, or gate
it. `check_manuscript.py` grew from 32 checks to 63, the reproduction record and
the arXiv comments field became generated, and the test counts are now collected
rather than typed.

Two of the fixes were themselves wrong on the first attempt and were caught by
negative-testing the gate — reintroducing the stale value and confirming the
build fails. Every gate added here was checked that way.
