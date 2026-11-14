#!/usr/bin/env python3
"""Emit the JHEP Stage-0 live-state audit.

Writes three files, all regenerable and all derived from the working tree
rather than from prose:

    audit/JHEP_LIVE_REPOSITORY_STATE.md
    audit/JHEP_RESULT_INVENTORY.json
    audit/JHEP_CLAIM_INPUTS.json

The inventory records, for every artifact a JHEP claim reads, its path, size,
SHA-256, the command that regenerates it, and whether it is tracked by git.
Nothing here computes science; it only records what the science left behind.

Usage:
    python scripts/emit_jhep_live_state.py [--repo .]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

# Artifacts a JHEP claim reads, with the command that regenerates each.
# (path, regeneration command, what it certifies)
RESULT_ARTIFACTS: list[tuple[str, str, str]] = [
    (
        "results/rank81/certificate.json",
        "python spinor_trace_bridge/scripts/run_rank81_certificate.py "
        "--archive <archive> --fitting-primes 32749,32719,32717 "
        "--holdout-primes 32713,32707 --seeds 11,22,33",
        "exact modular Jacobian, 83x126, rank 81, per-degree blocks, "
        "Euler homogeneity, over the sample x prime matrix",
    ),
    (
        "results/rank81/minor81_certificate.json",
        "python spinor_trace_bridge/scripts/run_minor81_certificate.py",
        "explicit 81x81 minor: pivot rows, pivot columns, nonzero determinant, "
        "two independent determinant routines",
    ),
    (
        "spinor_trace_bridge/results/bridge_validation.json",
        "cd spinor_trace_bridge && python -m pytest",
        "forward rank 126, kernel = anti-self-dual 126, image = gamma-traceless "
        "126, left inverse, equivariance",
    ),
    (
        "verification/COMMON_SAMPLE_REGISTRY.json",
        "python spinor_trace_bridge/scripts/run_comparison.py",
        "degree-resolved tensor/spinor span comparison on a common sample",
    ),
    (
        "results/10d_order8.json",
        "python scripts/run_10d.py --order 8",
        "degree-8 tensor invariant enumeration",
    ),
    (
        "results/10d_order10.json",
        "python scripts/run_10d.py --order 10",
        "degree-10 tensor invariant enumeration (atlas input)",
    ),
    (
        "results/10d_order12.json",
        "python scripts/run_10d.py --order 12",
        "degree-12 tensor invariant enumeration",
    ),
    (
        "results/10d_graph_catalog.json",
        "python scripts/generate_graph_catalog.py",
        "canonical graph catalogue underlying the tensor basis",
    ),
    (
        "results/stress_flow_exact_low_degree.json",
        "python scripts/stress_flow_pipeline.py --exact",
        "exact low-degree stress-flow map (degree-10 application input)",
    ),
    (
        "results/live_state.json",
        "python scripts/emit_jhep_live_state.py",
        "recorded repository state at generation time",
    ),
]

# Claims the JHEP manuscript intends to make, each bound to the artifacts that
# support it.  A claim with a missing artifact cannot enter the manuscript.
CLAIM_INPUTS: list[dict] = [
    {
        "claim_id": "JHEP-C1",
        "statement": "The tensor module Lambda^5_+ V and the gamma-traceless "
        "symmetric chiral bispinors Sym^2_{gamma-tr} S_+ are related by an "
        "explicit equivariant map, exact over the fields in which it is "
        "constructed.",
        "kind": "certified",
        "artifacts": ["spinor_trace_bridge/results/bridge_validation.json"],
        "tests": ["spinor_trace_bridge/tests/test_bridge.py"],
    },
    {
        "claim_id": "JHEP-C2",
        "statement": "The forward map has rank 126 on the self-dual channel, "
        "kernel exactly the anti-self-dual 126, and an exact left inverse whose "
        "composition with the forward map is the self-dual projector.",
        "kind": "certified",
        "artifacts": ["spinor_trace_bridge/results/bridge_validation.json"],
        "tests": ["spinor_trace_bridge/tests/test_bridge.py"],
    },
    {
        "claim_id": "JHEP-C3",
        "statement": "Independently constructed tensor and spinor invariant "
        "families span the same space at degrees 4, 6 and 10; at degree 8 the "
        "port-graph family alone is strictly contained and the structured "
        "tensor-word family supplies the missing direction.",
        "kind": "certified",
        "artifacts": ["verification/COMMON_SAMPLE_REGISTRY.json"],
        "tests": ["spinor_trace_bridge/tests/test_bridge.py"],
    },
    {
        "claim_id": "JHEP-C4",
        "statement": "The 83 selected invariant functions have exact modular "
        "Jacobian rank 81, which is an unconditional lower bound on the rank "
        "over Q because the gamma-traceless basis is integral.",
        "kind": "certified",
        "artifacts": [
            "results/rank81/certificate.json",
            "results/rank81/minor81_certificate.json",
        ],
        "tests": ["spinor_trace_bridge/tests/test_adversarial.py"],
    },
    {
        "claim_id": "JHEP-C5",
        "statement": "The generic functional rank is exactly 81: the lower "
        "bound is certified, the matching upper bound 126 - 45 = 81 is analytic.",
        "kind": "proved (upper bound analytic) + certified (lower bound)",
        "artifacts": [
            "results/rank81/certificate.json",
            "docs/RANK81_CHARACTERISTIC_ZERO_PROOF.md",
        ],
        "tests": ["spinor_trace_bridge/tests/test_adversarial.py"],
    },
    {
        "claim_id": "JHEP-C6",
        "statement": "Application: dim A10 = 14, dim D10 = 11, dim Q10 = 3, "
        "with compact representatives and a basis-independent cardinality bound.",
        "kind": "certified + proved (cardinality bound)",
        "artifacts": [
            "results/10d_order10.json",
            "results/stress_flow_exact_low_degree.json",
        ],
        "tests": [
            "tests/test_published_degree10_invariants.py",
            "tests/test_quotient_cardinality_bound.py",
        ],
    },
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git(repo: Path, *args: str) -> str:
    out = subprocess.run(
        ["git", *args], cwd=repo, capture_output=True, text=True, check=False
    )
    return out.stdout.strip()


def tracked_files(repo: Path) -> set[str]:
    return set(git(repo, "ls-files").splitlines())


def collect_repo_state(repo: Path) -> dict:
    tracked = tracked_files(repo)
    porcelain = git(repo, "status", "--porcelain")
    dirty = [line for line in porcelain.splitlines() if line.strip()]
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "head": git(repo, "rev-parse", "HEAD"),
        "head_short": git(repo, "rev-parse", "--short", "HEAD"),
        "head_subject": git(repo, "log", "-1", "--format=%s"),
        "head_committer_date": git(repo, "log", "-1", "--format=%cI"),
        "branch": git(repo, "rev-parse", "--abbrev-ref", "HEAD"),
        "remote": git(repo, "remote", "get-url", "origin"),
        "local_branches": [
            line.strip()
            for line in git(repo, "branch", "--format=%(refname:short) %(objectname)").splitlines()
        ],
        "remote_branches": [
            line.strip()
            for line in git(
                repo, "branch", "-r", "--format=%(refname:short) %(objectname)"
            ).splitlines()
            if "->" not in line
        ],
        "tags": [
            line.strip()
            for line in git(repo, "tag", "--list", "--format=%(refname:short) %(objectname)").splitlines()
        ],
        "stashes": git(repo, "stash", "list").splitlines(),
        "dirty_paths": dirty,
        "tracked_file_count": len(tracked),
        "describe": git(repo, "describe", "--tags", "--always"),
    }


def collect_inventory(repo: Path) -> list[dict]:
    tracked = tracked_files(repo)
    rows = []
    for rel, command, certifies in RESULT_ARTIFACTS:
        path = repo / rel
        row = {
            "path": rel,
            "certifies": certifies,
            "regeneration_command": command,
            "tracked": rel in tracked,
            "exists": path.exists(),
        }
        if path.exists():
            row["bytes"] = path.stat().st_size
            row["sha256"] = sha256(path)
            row["mtime_utc"] = (
                datetime.fromtimestamp(path.stat().st_mtime, timezone.utc)
                .isoformat(timespec="seconds")
            )
            row["status"] = "PRESENT"
        else:
            row["status"] = "MISSING"
        rows.append(row)
    return rows


def collect_test_suites(repo: Path) -> list[dict]:
    return [
        {
            "suite": "tensor side",
            "root": ".",
            "command": "python -m pytest tests/ -q",
            "test_files": sorted(p.name for p in (repo / "tests").glob("test_*.py")),
        },
        {
            "suite": "bridge side",
            "root": "spinor_trace_bridge",
            "command": "cd spinor_trace_bridge && python -m pytest",
            "test_files": sorted(
                p.name for p in (repo / "spinor_trace_bridge" / "tests").glob("test_*.py")
            ),
        },
    ]


def render_markdown(state: dict, inventory: list[dict], suites: list[dict]) -> str:
    lines: list[str] = []
    A = lines.append
    A("# JHEP Stage 0 --- live repository state")
    A("")
    A(f"Generated {state['generated_utc']} by `scripts/emit_jhep_live_state.py`.")
    A("Every value below is read from the working tree, not from a prior report.")
    A("")
    A("## Head")
    A("")
    A("| field | value |")
    A("|---|---|")
    A(f"| branch | `{state['branch']}` |")
    A(f"| local HEAD | `{state['head']}` |")
    A(f"| subject | {state['head_subject']} |")
    A(f"| committed | {state['head_committer_date']} |")
    A(f"| remote | {state['remote']} |")
    A(f"| tracked files | {state['tracked_file_count']} |")
    A(f"| stashes | {len(state['stashes']) or 'none'} |")
    A("")
    A("## Branches")
    A("")
    A("### Local")
    A("")
    for b in state["local_branches"]:
        A(f"- `{b}`")
    A("")
    A("### Remote")
    A("")
    for b in state["remote_branches"]:
        A(f"- `{b}`")
    A("")
    A("## Tags")
    A("")
    for t in state["tags"]:
        A(f"- `{t}`")
    A("")
    A("## Uncommitted paths")
    A("")
    if state["dirty_paths"]:
        for d in state["dirty_paths"]:
            A(f"- `{d}`")
    else:
        A("None.")
    A("")
    A("## Result artifacts a JHEP claim reads")
    A("")
    A("| artifact | status | bytes | sha256 (first 16) | tracked |")
    A("|---|---|---|---|---|")
    for row in inventory:
        digest = row.get("sha256", "")[:16]
        size = row.get("bytes", "")
        A(
            f"| `{row['path']}` | {row['status']} | {size} | `{digest}` | "
            f"{'yes' if row['tracked'] else 'no'} |"
        )
    A("")
    A("Regeneration commands and what each artifact certifies are in")
    A("`audit/JHEP_RESULT_INVENTORY.json`.")
    A("")
    A("## Test suites")
    A("")
    for s in suites:
        A(f"### {s['suite']}")
        A("")
        A(f"    {s['command']}")
        A("")
        A(f"{len(s['test_files'])} test files.")
        A("")
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".", type=Path)
    args = ap.parse_args()
    repo = args.repo.resolve()

    state = collect_repo_state(repo)
    inventory = collect_inventory(repo)
    suites = collect_test_suites(repo)

    audit = repo / "audit"
    audit.mkdir(exist_ok=True)

    (audit / "JHEP_LIVE_REPOSITORY_STATE.md").write_text(
        render_markdown(state, inventory, suites), encoding="utf-8"
    )
    (audit / "JHEP_RESULT_INVENTORY.json").write_text(
        json.dumps(
            {"repository_state": state, "artifacts": inventory, "test_suites": suites},
            indent=1,
        )
        + "\n",
        encoding="utf-8",
    )

    missing = [r["path"] for r in inventory if r["status"] == "MISSING"]
    claim_rows = []
    present = {r["path"] for r in inventory if r["status"] == "PRESENT"}
    for claim in CLAIM_INPUTS:
        unmet = [a for a in claim["artifacts"] if a not in present and not (repo / a).exists()]
        claim_rows.append({**claim, "missing_artifacts": unmet,
                           "input_status": "READY" if not unmet else "BLOCKED"})
    (audit / "JHEP_CLAIM_INPUTS.json").write_text(
        json.dumps(
            {
                "generated_utc": state["generated_utc"],
                "head": state["head"],
                "claims": claim_rows,
                "missing_artifacts": missing,
            },
            indent=1,
        )
        + "\n",
        encoding="utf-8",
    )

    blocked = [c["claim_id"] for c in claim_rows if c["input_status"] == "BLOCKED"]
    print(f"head {state['head_short']} on {state['branch']}")
    print(f"artifacts present {len(present)}/{len(inventory)}")
    if missing:
        print("missing: " + ", ".join(missing))
    print(f"claims ready {len(claim_rows) - len(blocked)}/{len(claim_rows)}")
    if blocked:
        print("blocked: " + ", ".join(blocked))
    return 0


if __name__ == "__main__":
    sys.exit(main())
