"""Verify a local copy of the third-party spinor archive against the manifest.

The archive is not redistributed, so a reader supplies their own copy.  This
checks that the copy is the same one the results were computed from, before any
archive-dependent script is run against it.

    python scripts/verify_archive.py --archive PATH \
        --manifest release_candidate/spinor-archive/MANIFEST.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True)
    ap.add_argument("--manifest", required=True)
    args = ap.parse_args()

    archive = Path(args.archive).expanduser().resolve()
    manifest = json.loads(Path(args.manifest).read_text())

    missing, mismatched, ok = [], [], 0
    for entry in manifest["files"]:
        p = archive / entry["path"]
        if not p.is_file():
            missing.append(entry["path"])
            continue
        if sha256(p) != entry["sha256"]:
            mismatched.append(entry["path"])
        else:
            ok += 1

    combined = hashlib.sha256(
        "".join(e["sha256"] for e in manifest["files"]).encode()).hexdigest()

    print(f"archive: {archive}")
    print(f"  matched   : {ok} / {manifest['n_files']}")
    print(f"  missing   : {len(missing)}")
    print(f"  mismatched: {len(mismatched)}")
    for p in (missing + mismatched)[:10]:
        print(f"    - {p}")

    if missing or mismatched:
        print("VERIFICATION FAILED: this is not the archive the results came from.")
        print("Archive-dependent results will not reproduce. All other results will.")
        return 1
    print(f"verified; combined sha256 {combined[:32]}...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
