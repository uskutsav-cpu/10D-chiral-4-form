#!/usr/bin/env python3
"""Place the arXiv and journal source archives at their frozen release paths.

The archives themselves are built and isolation-tested by
`scripts/build_submission_package.py`, which refuses to emit them unless a build
from the staged files alone is clean. This step only freezes them: copies them
under `release/`, writes a `.sha256` beside each, and records what the packages
must satisfy so the checks are auditable rather than remembered.

Nothing is uploaded.

    python scripts/freeze_submission_packages.py
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import tarfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "submission_candidate"
ARXIV = ROOT / "release" / "arxiv"
JHEP = ROOT / "release" / "jhep"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main() -> int:
    ARXIV.mkdir(parents=True, exist_ok=True)
    JHEP.mkdir(parents=True, exist_ok=True)

    tar = SRC / "arxiv_source.tar.gz"
    zip_ = SRC / "jhep_source.zip"
    for p in (tar, zip_):
        if not p.exists():
            print(f"missing {p.name}; run scripts/build_submission_package.py")
            return 1

    names_tar = sorted(tarfile.open(tar).getnames())
    names_zip = sorted(zipfile.ZipFile(zip_).namelist())

    problems = []
    if names_tar != names_zip:
        problems.append("the two archives do not contain the same file set")
    if "main.tex" not in names_tar:
        problems.append("master main.tex is not at the archive root")
    if "main.bbl" not in names_tar:
        problems.append("main.bbl absent; the journal requires it with BibTeX")
    for bad in ("cover_letter.md", "cover_letter.pdf"):
        if any(bad in n for n in names_tar):
            problems.append(f"{bad} must not be inside the source archive")
    if not any(n.startswith("figures/") for n in names_tar):
        problems.append("no figures in the archive")

    dest_tar = ARXIV / "arxiv_source.tar.gz"
    dest_zip = JHEP / "jhep_submission_source.tar.gz"
    shutil.copy2(tar, dest_tar)
    # The journal package is frozen as a tarball at the requested path; the zip
    # built alongside it is kept too, since submission systems differ.
    shutil.copy2(tar, dest_zip)
    shutil.copy2(zip_, JHEP / "jhep_source.zip")

    for p in (dest_tar, dest_zip, JHEP / "jhep_source.zip"):
        (p.parent / (p.name + ".sha256")).write_text(
            f"{sha256(p)}  {p.name}\n")

    manifest = json.loads((SRC / "package_manifest.json").read_text())
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    record = {
        "schema": 1,
        "generated_by": "scripts/freeze_submission_packages.py",
        "commit": head,
        "files_in_archives": names_tar,
        "identical_file_sets": names_tar == names_zip,
        "isolated_build": manifest.get("build"),
        "comments_field": manifest.get("comments_field"),
        "checks": {
            "master_tex_at_root": "main.tex" in names_tar,
            "bbl_included": "main.bbl" in names_tar,
            "no_cover_letter_inside": not any("cover_letter" in n for n in names_tar),
            "figures_present": any(n.startswith("figures/") for n in names_tar),
        },
        "problems": problems,
        "uploaded": False,
        "submitted": False,
    }
    (ROOT / "release" / "SUBMISSION_FREEZE.json").write_text(
        json.dumps(record, indent=1) + "\n")

    print(f"frozen {len(names_tar)} files; identical sets: {names_tar == names_zip}")
    print(f"isolated build: {manifest.get('build')}")
    if problems:
        print("PROBLEMS:")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("arXiv and journal packages frozen; nothing uploaded")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
