# Data Contract with the Invariant-Atlas Repository

This repository should consume verified artifacts from `selfdual-5form-invariants` without mutating the atlas engine.

## Required imported objects

At minimum, a frozen source snapshot should identify:
- degree-4,6,8,10,12 invariant bases;
- basis ordering/fingerprint;
- exact evaluation interface on a self-dual five-form sample;
- interacting stress-flow coefficient artifact;
- static stress rows;
- intrinsic sextic change-of-basis artifact;
- source conventions.

## Recommended local layout

When both repositories are checked out side by side:

```text
workspace/
  selfdual-5form-invariants/
  10D-chiral-4-form/
```

Run:

```bash
python scripts/freeze_source_snapshot.py ../selfdual-5form-invariants
```

The command writes only hashes/metadata into this repository by default; it does not silently copy large generated catalogs.

## Imported claim rule

An imported numerical table is a **provenance pointer**, not an independent reproduction. `data/baseline/paper2_baseline.json` therefore carries `status = imported` until the relevant source pipeline has been rerun and a new certificate is generated here.

## Future API

The desired stable upstream interface is a small export containing:
- ordered basis labels and graph/tensor records;
- exact evaluation callable or serialized evaluation recipe;
- degree metadata;
- convention ID;
- semantic fingerprint.

Paper 2 should not depend on private paths or ephemeral caches in the predecessor repository.
