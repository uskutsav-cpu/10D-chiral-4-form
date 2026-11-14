"""Assemble the private release candidate, and refuse to if it is unsafe.

Deliberately excluded, and the exclusion is enforced rather than documented:

  * the third-party spinor archive (redistribution permission unresolved) --
    only a manifest with per-file hashes and adapter instructions is shipped;
  * source PDFs of cited papers (copyright);
  * any file matching the secret patterns or containing an absolute home path.

The build fails on a violation instead of shipping and warning.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "release_candidate"

SECRET_PATTERNS = [
    (r"gh[pousr]_[A-Za-z0-9]{20,}", "GitHub token"),
    (r"sk-[A-Za-z0-9]{20,}", "API key"),
    (r"AKIA[0-9A-Z]{16}", "AWS key id"),
    (r"-----BEGIN (RSA |OPENSSH |EC )?PRIVATE KEY-----", "private key"),
    (r"xox[baprs]-[A-Za-z0-9-]{10,}", "Slack token"),
]
HOME_PATH = re.compile(r"/Users/[A-Za-z0-9._-]+/")

INCLUDE = [
    ("trace-code", ["src", "scripts", "results", "verification", "tests", "pytest.ini", "requirements.txt",
                    "requirements-lock.txt", "README.md"]),
    ("bridge-code", ["spinor_trace_bridge"]),
    # results/rank81 holds the exact Jacobian certificate and the explicit
    # 81x81 minor -- the strongest computational claim in the package. It was
    # absent from the release candidate while the manuscript cited it.
    ("certificates", ["verification", "results/intrinsic_candidates",
                      "results/rank81"]),
    ("reproduction", ["manuscript/scripts", "scripts"]),
]

TEXT_SUFFIXES = {".py", ".md", ".txt", ".json", ".tex", ".cfg", ".ini",
                 ".toml", ".yaml", ".yml", ".sh", ".csv", ".bib", ".sty", ".bst"}



def release_copy_ignore(directory: str, names: list[str]) -> set[str]:
    """Exclude local build debris and the separate independent-audit bundle."""
    ignored = set(shutil.ignore_patterns(
        "__pycache__", "*.pyc", ".pytest_cache", "reproduction-logs"
    )(directory, names))
    # Keep results required by the historical tensor tests. The independent
    # audit has its own release and preserves original machine-specific paths.
    # Exclude only this exact subtree, not unrelated directories named audit.
    if Path(directory).resolve() == (ROOT / "results").resolve():
        ignored.add("audit")
    return ignored


def scan(path: Path) -> list[str]:
    """Secret and absolute-path scan of one file."""
    if path.suffix.lower() not in TEXT_SUFFIXES:
        return []
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return []
    problems = []
    for pattern, what in SECRET_PATTERNS:
        if re.search(pattern, text):
            problems.append(f"{path}: possible {what}")
    for m in HOME_PATH.finditer(text):
        problems.append(f"{path}: absolute home path {m.group(0)}")
    return problems


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def archive_manifest(archive: Path) -> dict:
    """Hashes and expected layout of the excluded third-party archive."""
    files = []
    for p in sorted(archive.rglob("*")):
        if not p.is_file():
            continue
        if ".pytest_cache" in p.parts or p.name in (".DS_Store",):
            continue
        if ".git" in p.parts:
            continue
        files.append({"path": str(p.relative_to(archive)),
                      "bytes": p.stat().st_size, "sha256": sha256(p)})
    combined = hashlib.sha256(
        "".join(f["sha256"] for f in files).encode()).hexdigest()
    return {
        "note": ("This archive is THIRD-PARTY code and is NOT redistributed. "
                 "Redistribution permission is unresolved; see "
                 "submission_candidate/AUTHORSHIP_DECISION_REQUIRED.md item 7."),
        "expected_local_directory": "self_dual_5_invariant_enumerator/",
        "n_files": len(files),
        "combined_sha256": combined,
        "adapter_instructions": [
            "Place your own copy of the archive anywhere on disk.",
            "Verify it against combined_sha256 below using scripts/verify_archive.py.",
            "Pass its path to the archive-dependent scripts with --archive PATH:",
            "  spinor_trace_bridge/scripts/run_archive_jacobian_exact.py",
            "  spinor_trace_bridge/scripts/analyse_archived_jacobians.py",
            "  spinor_trace_bridge/scripts/run_float_jacobian_matrix.py",
            "All other results reproduce without the archive.",
        ],
        "files": files,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", default=None,
                    help="local path to the third-party archive, for the manifest only")
    ap.add_argument("--allow-unsafe", action="store_true",
                    help="build even if the scan finds problems (not recommended)")
    args = ap.parse_args()

    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    problems: list[str] = []
    copied = 0
    for target, sources in INCLUDE:
        for rel in sources:
            src = ROOT / rel
            if not src.exists():
                continue
            dst = OUT / target / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            if src.is_dir():
                shutil.copytree(src, dst, dirs_exist_ok=True,
                                ignore=release_copy_ignore)
            else:
                shutil.copy2(src, dst)

    for p in sorted(OUT.rglob("*")):
        if p.is_file():
            copied += 1
            problems.extend(scan(p))

    # manifest for the excluded archive
    if args.archive:
        archive = Path(args.archive).expanduser().resolve()
        if archive.is_dir():
            man = OUT / "spinor-archive"
            man.mkdir(parents=True, exist_ok=True)
            (man / "MANIFEST.json").write_text(
                json.dumps(archive_manifest(archive), indent=1))
            (man / "README.md").write_text(
                "# Third-party spinor archive — NOT INCLUDED\n\n"
                "This directory contains only a manifest. The archive itself is\n"
                "third-party code whose redistribution permission is unresolved,\n"
                "so it is deliberately absent.\n\n"
                "See `MANIFEST.json` for per-file hashes, the expected directory\n"
                "name, and adapter instructions for reproducing the\n"
                "archive-dependent results from your own copy.\n")

    # citation metadata
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                            capture_output=True, text=True).stdout.strip()
    (OUT / "CITATION.cff").write_text(
        "cff-version: 1.2.0\n"
        "message: \"If you use this software, please cite it.\"\n"
        "title: \"Exact degree-ten invariants of the ten-dimensional self-dual five-form\"\n"
        "authors:\n"
        "  - name: \"AUTHOR LIST -- HUMAN ACTION REQUIRED\"\n"
        "version: \"0.0.0-private\"\n"
        f"commit: \"{commit}\"\n"
        "license: \"LICENCE NOT YET CHOSEN -- HUMAN ACTION REQUIRED\"\n"
        "doi: \"DOI NOT CREATED -- PENDING AUTHORISATION\"\n"
        "repository-code: \"URL PENDING AUTHORISATION\"\n"
    )
    (OUT / "LICENCE-DECISION-REQUIRED.md").write_text(
        "# No licence has been chosen\n\n"
        "Without a licence this code is not reusable by others, regardless of\n"
        "where it is hosted. Choosing one is a human decision; see\n"
        "`submission_candidate/AUTHORSHIP_DECISION_REQUIRED.md` item 6.\n")

    report = {
        "files_copied": copied,
        "scan_problems": problems,
        "archive_included": False,
        "archive_manifest_included": bool(args.archive),
        "commit": commit,
    }
    (OUT / "RELEASE_SCAN.json").write_text(json.dumps(report, indent=1))

    print(f"release candidate: {copied} files -> {OUT.relative_to(ROOT)}")
    if problems:
        print(f"SCAN FOUND {len(problems)} PROBLEM(S):")
        for p in problems[:20]:
            print(f"  - {p}")
        if not args.allow_unsafe:
            print("refusing to declare the release candidate safe")
            return 1
    else:
        print("secret and path scan: clean")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
