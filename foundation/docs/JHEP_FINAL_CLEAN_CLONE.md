# Final clean-clone reproduction

Machine-readable: `results/jhep/final_clean_clone.json`. Actual counts, not
expected ones.

Fresh `git clone` into `/tmp` — outside iCloud — with a new virtualenv built from
`requirements-lock.txt` plus the documented extras (`opt_einsum`, `networkx`,
`joblib`, `python-igraph`, `pynauty`, `matplotlib`). No bytecode carried over: the
clone contained zero `__pycache__` entries.

## Results

| step | commit | result |
|---|---|---|
| tensor test suite | `ba311a0` | **224 passed**, 0 failures |
| bridge test suite | `7e4ea08` | **72 passed**, 0 failures |
| manuscript gates | `7e4ea08` | **72 checks**, all passed |
| regenerate macros | `7e4ea08` | 122 lines, **0 artifacts missing** |
| regenerate tables | `7e4ea08` | 5 |
| regenerate identity | `7e4ea08` | from the generator artifact |
| regenerate figures | `7e4ea08` | 7 |
| release policy check | `7e4ea08` | facts match the repository |
| isolated manuscript build | `7e4ea08` | **24 pages, 0 errors, 0 undefined citations, 0 undefined references, 0 overfull boxes** |
| release candidate | `7e4ea08` | 173 files, **secret and path scan clean** |
| archive hashes vs working tree | `7e4ea08` | **identical**, both archives |

## Why two commits appear

The tensor suite is the expensive step and ran at `ba311a0`. Four commits landed
while it was running: `99b5aa5`, `fda76a8`, `b77ef4d`, `7e4ea08`. Every one was
inspected: none touches a file under `tests/`, `src/`,
`spinor_trace_bridge/src/`, `results/` or `verification/`, so no test can observe
them. They change figure determinism, archive packaging, audit prose and
generated artifacts.

Every other step was re-run in a second clone at `7e4ea08` itself. Reporting one
commit for everything would have been simpler and false.

## The archives are byte-reproducible

Both source archives built in an independent clone hash **identically** to the
working tree's. That was not true when this phase began: a first clean clone
produced different archive hashes, and comparing them file by file showed the
same 25-file set, all eighteen non-figure files byte-identical, and all seven
figures differing solely in matplotlib's embedded `/CreationDate`.

The content was fine and the comparison was useless — the worst combination,
because a hash check that fails on every rebuild teaches you to ignore it. Fixed
by pinning `SOURCE_DATE_EPOCH` for the figures and normalising tar/zip
timestamps, ownership and modes. A hash mismatch from here means the content
changed.

## What a clean clone does not regenerate

- `results/rank81/certificate.json`
- `results/rank81/minor81_certificate.json`
- `verification/SPINOR_JACOBIAN_RUNS.json`

These need a local copy of the third-party spinor archive, which is not
redistributed. The certificates themselves are committed, so every number in the
manuscript is present and both manuscripts build without it. Only their
*regeneration* needs the archive; `scripts/verify_archive.py` checks a reader's
own copy against a manifest of per-file hashes.

## Requirements, each actually checked

zero failed tests · zero LaTeX errors · zero undefined citations · zero undefined
references · zero overfull boxes · zero absolute private paths · zero
mentor-archive files · zero secret-scan hits · zero unlisted tracked files above
the release-policy threshold.
