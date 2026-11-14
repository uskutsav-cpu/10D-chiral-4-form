#!/usr/bin/env python3
"""Check the claims in docs/PUBLIC_RELEASE_POLICY.md against the repository.

The policy asserted "nothing tracked above 1 MB" and named the compiled
manuscript as the largest tracked file at 378 KB. Both were false: a 1.2 MB
result file is tracked and the manuscript had grown to 434 KB. A policy document
that states a fact about the repository should be checked against it, not
maintained by hand.

    python scripts/check_release_policy.py            # check
    python scripts/check_release_policy.py --write    # update the stated facts
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "docs" / "PUBLIC_RELEASE_POLICY.md"
LIMIT = 1024 * 1024

#: Files allowed to exceed the limit, each with the reason it is not regenerable
#: cheaply. Adding to this list is a deliberate act; growing past the limit by
#: accident is what the check is for.
ALLOWED_LARGE = {
    "results/stress_flow/interacting_flow_equations.json":
        "the interacting flow equations through degree 12; regenerating them "
        "costs hours and every downstream certificate is keyed to them",
}


def tracked_sizes() -> list[tuple[int, str]]:
    out = subprocess.run(["git", "ls-files"], cwd=ROOT,
                         capture_output=True, text=True).stdout.split("\n")
    sizes = []
    for rel in out:
        if not rel.strip():
            continue
        p = ROOT / rel
        if p.is_file():
            sizes.append((p.stat().st_size, rel))
    return sorted(sizes, reverse=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()

    sizes = tracked_sizes()
    if not sizes:
        print("no tracked files found; is this a git repository?")
        return 1

    biggest_size, biggest = sizes[0]
    over = [(s, r) for s, r in sizes if s > LIMIT]
    unexpected = [(s, r) for s, r in over if r not in ALLOWED_LARGE]

    problems = []
    if unexpected:
        for s, r in unexpected:
            problems.append(f"{r} is {s/1024:.0f} KB, over the 1 MB limit and "
                            f"not in the allowed list")

    text = POLICY.read_text() if POLICY.exists() else ""
    stated = re.search(r"Largest tracked file is[^.]*?(\d[\d,]*)\s*KB", text)
    if stated:
        claimed = int(stated.group(1).replace(",", ""))
        actual = round(biggest_size / 1024)
        if claimed != actual:
            problems.append(f"policy states the largest tracked file is "
                            f"{claimed} KB; it is {biggest} at {actual} KB")
    if "nothing tracked above 1 MB" in text and over:
        problems.append("policy states nothing is tracked above 1 MB, but "
                        + ", ".join(r for _, r in over) + " is")

    for s, r in sizes[:8]:
        print(f"  {s/1024:8.1f} KB  {r}")

    if args.write and POLICY.exists():
        listed = "\n".join(
            f"- `{r}` — {ALLOWED_LARGE[r]}" for _, r in over if r in ALLOWED_LARGE)
        block = (
            f"Largest tracked file is `{biggest}` at "
            f"{biggest_size/1024:.0f} KB.\n\n"
            + (f"{len(over)} tracked file(s) exceed 1 MB, each deliberately:\n\n"
               f"{listed}\n\n" if over else "Nothing tracked exceeds 1 MB.\n\n")
            + "These figures are checked by `scripts/check_release_policy.py`, "
              "which fails if a tracked file grows past the limit without being "
              "listed. They were maintained by hand once and were wrong.\n")
        new = re.sub(r"Largest tracked file is.*?(?=\n## |\Z)", block, text,
                     flags=re.S)
        POLICY.write_text(new)
        print("policy facts updated")
        return 0

    if problems:
        print("\nRELEASE POLICY CHECK FAILED:")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("\nrelease policy facts match the repository")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
