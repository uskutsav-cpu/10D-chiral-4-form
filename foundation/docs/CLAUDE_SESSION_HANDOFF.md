# Session handoff

Last updated 2026-08-01. Branch `research/maximal-chiral-four-form-program`.

## Where things stand

A complete private scientific and technical submission candidate exists. The
manuscript compiles from the official journal class with zero errors and zero
undefined citations or references, every scientific number in it is generated
from a JSON certificate, and a fresh clone reproduces all of it.

## The findings a successor should know

1. **The signature blocker was wrong.** The oscillator frame is split `(5,5)`,
   not Euclidean. Real self-dual five-forms exist there, `⋆² = +1`, and the
   bridge was buildable all along. See `spinor_trace_bridge/docs/REAL_FORM_DICTIONARY.md`.

2. **`12 = 12` was a coincidence.** The published degree-ten span is not a
   complement to the products; it contains one product direction, so its
   non-product content is eleven. The structural decomposition is
   `A10 = G10 ⊕ P10`. See `docs/DIMENSION_DICTIONARY.md`.

3. **Degree eight disagrees, and the drop test says why.** The spinor
   port-graph stream reaches 6 where the tensor side reaches 7, with strict
   containment. Run with the *full* spinor family the two sides give
   `trace = spinor = union = 7` at four primes, and dropping each family in turn
   shows only the tensor words are indispensable. So the missing direction is not
   reachable by any single-graph contraction. See
   `verification/DEGREE8_SPAN_EQUALITY.md`. It is a positive result, not a bug.

4. **Both characteristic-zero questions are closed.** `dim_Q D10 = 11` and
   `dim_Q(B10 cap P10) = 1`, exactly, not as bounds. The route in both cases was
   the same: CRT lift + rational reconstruction + validation at a held-out prime,
   then exact rational linear algebra. The intersection has an explicit integer
   generator, verified at a further prime. See
   `docs/D10_Q10_CHARACTERISTIC_ZERO_STATUS.md` and
   `docs/B10_P10_INTERSECTION_STATUS.md`.

5. **The exact Jacobian now covers all 83 candidates and reaches rank 81.** An
   earlier revision of this file said 59 for the port-graph subset; that is
   superseded. See `docs/RANK81_EXACT_CERTIFICATE.md`. The bound is
   `rank_Q >= 81` unconditionally, witnessed by an explicit 81x81 minor. The
   matching upper bound is analytic and stays that way.

## What is deliberately incomplete

| item | state |
|---|---|
| float64 seed/scale/step Jacobian matrix | not run; one configuration costs >10 min. Replaced by an exact analytic Jacobian, which is stronger |
| degree-10 no-stop terminal status | see `verification/SPINOR_DEGREE10_NO_STOP.md`; cluster script prepared, no cluster run has occurred |
| certified rational reconstruction | not done; no claim depends on it, and every modular result is used as a lower bound only |
| PO-09 | **CERTIFIED** for every degree-ten claim: `dim_Q D10 = 11`, `dim_Q Q10 = 3`, `dim_Q(B10 cap P10) = 1`, each by an exact rational lift validated at a held-out prime |
| PO-03 / PO-05 / PO-07 | `NOT APPLICABLE` to this manuscript — checked, not assumed; see `audit/PROOF_OBLIGATIONS_FINAL.md` |
| PO-08 | now split: the cardinality bound is **discharged analytically**; removal-minimality under general `GL` is still open |
| novelty rows | all PROVISIONAL — the literature sweep is *done* (`audit/RELATED_WORK_COMPLETE.md`); what is missing is coauthor confirmation, which is a human gate |

## Long runs: check for orphans before blaming the machine

Several long runs here died silently, with the process simply gone and no
traceback, and the first diagnosis — 8 GB of RAM, too many heavy jobs — was
**wrong**. `ps` eventually showed the real cause: two separate rank-81 processes
running the same command against the same row cache, plus three abandoned
`pytest` processes from runs that had been reported as dead but were not. They
were competing for CPU, and each new run made it worse.

A background job that stops producing output has not necessarily stopped, and a
job **reported as completed** has not necessarily taken its children with it —
two `pytest` processes were found alive twenty minutes after the runs that
spawned them were reported finished with exit code 0. Check

    ps aux | grep -E "pytest|spinor_trace_bridge"

before concluding anything about resources, and kill what you find. This had to
be done twice in one session; after each cleanup the same commands ran fine.

`pytest` writing to a redirected file also buffers, so a log that has not grown
for several minutes is not evidence of a hang. Check the process's CPU
percentage instead.

The machine does have 8 GB and the degree-12 contractions are not small, so it is
still worth not stacking heavy jobs. But that was not what killed these.

Everything long resumes rather than restarts: the rank-81 runner keys its row
cache by prime, seed, candidate id and formula hash, and the incidence generator
skips primes already recorded. Re-issuing the identical command is the correct
response to a killed run.

## Where the human gates are

- `submission_candidate/AUTHORSHIP_DECISION_REQUIRED.md` — 22 unticked items
- `spinor_trace_bridge/docs/MENTOR_REVIEW_ITEMS.md` — G-1 to G-9
- `audit/NOVELTY_MATRIX.md` — every row PROVISIONAL

## How to pick this up

```bash
python -m pytest                                  # 254
cd spinor_trace_bridge && python -m pytest        # 86
python manuscript/scripts/make_numbers.py         # regenerate from artifacts
python manuscript/scripts/make_tables.py
python manuscript/scripts/make_figures.py
python scripts/build_submission_package.py        # isolated build + archives
python manuscript/scripts/check_manuscript.py     # wording, claim and doc gates
python scripts/reproduce_all.py                   # all of the above, recorded
```

Long runs write incrementally and skip completed work, so re-issuing the same
command resumes rather than restarts.

## Warnings

- The scratchpad is wiped between sessions. The working virtualenv lives at
  `../.venv`, outside the repository.
- The third-party spinor archive at `../spinor-work` is **not** redistributable
  and is excluded from the repository and the release candidate.
- Run the bridge suite from inside `spinor_trace_bridge/`. Invoking it from the
  repository root with `-c` collects the wrong tests.
