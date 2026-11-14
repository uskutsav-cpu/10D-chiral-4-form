# Public release policy

## Authorised for public push

- User-authored source, tests, derivations, documentation, claim ledgers.
- Manuscript source, figures, tables, bibliography.
- Compact certificates and result matrices needed for reproduction.
- Environment locks, reproduction scripts, release notes.
- Project-authored audit summaries describing the external archive **without**
  redistributing its code or substantial contents.
- Annotated tags and source releases tied to verified commits.

## Prohibited

- The mentor archive in any form: ZIP, source, work directories, substantial
  extracts, raw run outputs.
- Secrets, tokens, keys, credentials.
- Personal absolute paths, caches, bytecode, temporary checkpoints.
- Regenerable large arrays in ordinary Git history.
- Any artifact whose authorship or licence cannot be established.

## Enforcement

| control | mechanism |
|---|---|
| secret scan | pattern match over every text file; build **fails**, not warns |
| absolute path scan | `/Users/<name>/` pattern; build fails |
| archive exclusion | archive lives outside the repository tree; allowlist copy |
| large files | `scripts/check_release_policy.py` fails if a tracked file exceeds 1 MB without being named, with its reason, in that script's allowed list |
| history rewriting | **never force-push**; published scientific history is immutable |

## Current state

Largest tracked file is `results/stress_flow/interacting_flow_equations.json` at 1223 KB.

1 tracked file(s) exceed 1 MB, each deliberately:

- `results/stress_flow/interacting_flow_equations.json` — the interacting flow equations through degree 12; regenerating them costs hours and every downstream certificate is keyed to them

These figures are checked by `scripts/check_release_policy.py`, which fails if a tracked file grows past the limit without being listed. They were maintained by hand once and were wrong.

## Licence

**Not yet chosen.** Without one the code is not reusable regardless of where it
is hosted. This is a human decision recorded in
`submission_candidate/AUTHORSHIP_DECISION_REQUIRED.md` item 6, and the release
carries `LICENCE-DECISION-REQUIRED.md` rather than a guessed licence.
