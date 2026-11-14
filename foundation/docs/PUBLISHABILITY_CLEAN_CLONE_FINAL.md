# Final clean-clone reproduction

Machine-readable: `results/jhep/publishability_clean_clone_final.json`.
Actual counts, not expected ones.

Fresh clone into `/tmp` — outside iCloud — with a virtualenv built from
`requirements-lock.txt` plus the documented extras.

| step | result |
|---|---|
| tensor test suite | **254 passed**, 0 failures, exit 0 |
| bridge test suite | **86 passed**, 0 failures, exit 0 |
| aggregator self-test | **11 rejections fired**; cells unmodified; deterministic and order-independent |
| G-10 and D10 certificate tests | 34 passed |
| manuscript gates | **72 checks**, all passed |
| isolated manuscript build | **25 pages, 0 errors, 0 undefined citations, 0 undefined references, 0 overfull boxes** |
| release candidate | 199 files, **secret and path scan clean** |
| release policy | facts match the repository |
| archives rebuilt twice in the clone | **byte-identical** |
| archives vs working tree | **byte-identical** |

## Byte reproducibility

Both source archives rebuild to identical bytes inside the clone, and match the
working tree's exactly. This was not true when the goal began: a rebuild changed
every archive hash, because matplotlib stamps a creation time into each figure
and tar/zip record mtimes. Fixed by pinning `SOURCE_DATE_EPOCH` and normalising
the archive metadata.

The point is not tidiness. A hash check that fails on every rebuild teaches you
to ignore it; one that only fails on a content change is worth reading.

## Which commit

The suites ran at `69ee992`. Later commits change audit prose, the review
packages and submission metadata only — no file under `tests/`, `src/`,
`spinor_trace_bridge/src/`, `results/` or `verification/`. That was checked
commit by commit, not assumed.

## The one requirement not met, and why

`zero_unlicensed_third_party_source` is recorded as **false**, and the reason is
not the third-party files. `jheppub.sty` and `JHEP.bst` are SISSA Medialab's,
under LPPL 1.3 or later, and redistributing them in a submission archive is their
intended use.

The flag is false because **this repository has no licence at all**. Until one is
chosen, nothing here is reusable by anyone and the release cannot honestly be
called open-source software. See `audit/LICENSING_AND_PROVENANCE_FINAL.md`.

## Not regenerated without the third-party archive

`results/rank81/certificate.json`, `results/rank81/minor81_certificate.json` and
`verification/SPINOR_JACOBIAN_RUNS.json`. All are committed, so every number in
the manuscript is present and both manuscripts build without the archive; only
their regeneration needs it.
