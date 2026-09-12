#!/usr/bin/env python3
"""Freeze hashes of critical files from a local predecessor checkout.

This records provenance only. It does not copy or reinterpret scientific output.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from chiral4form.provenance import sha256_file

DEFAULT_FILES = [
    "results/10d_order12.json",
    "results/degree12_benchmarks.json",
    "docs/mentor_review_package.md",
    "docs/assumptions_limitations_and_open_questions.md",
]


def git_head(repo: Path) -> str:
    return subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_repo", type=Path)
    parser.add_argument("--out", type=Path, default=Path("data/baseline/source_manifest.json"))
    args = parser.parse_args()
    repo = args.source_repo.resolve()
    if not (repo / ".git").exists():
        raise SystemExit(f"not a git checkout: {repo}")

    files = []
    for rel in DEFAULT_FILES:
        path = repo / rel
        files.append(
            {
                "path": rel,
                "exists": path.exists(),
                "sha256": sha256_file(path) if path.exists() else None,
            }
        )
    payload = {
        "source_repo": str(repo),
        "git_head": git_head(repo),
        "files": files,
    }
    out = args.out
    if not out.is_absolute():
        out = Path.cwd() / out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(out)


if __name__ == "__main__":
    main()
