#!/usr/bin/env python3
"""Run the full reproduction sequence and record what actually happened.

The reproduction record was written by hand twice and went stale both times --
it claimed 199 and 49 tests when the suites were 207 and 72, and 32 gates when
there were 60. A record of a run should be produced by the run.

This script executes each step, captures the real numbers, and writes
`verification/REPRODUCTION_RECORD.json` plus a markdown table beside it. It exits
nonzero if any step fails, so it is a check and not only a report.

    python scripts/reproduce_all.py                  # full sequence
    python scripts/reproduce_all.py --skip-tests     # everything but the suites
    python scripts/reproduce_all.py --skip-build     # no LaTeX required

Archive-dependent artifacts are NOT regenerated: they need a copy of the
third-party spinor archive that is not redistributed. They are committed, so
every number in both manuscripts is present without it. The record says which
steps were skipped rather than implying a fuller run than occurred.
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VENV = ROOT.parent / ".venv" / "bin" / "python"
PY = str(VENV) if VENV.exists() else sys.executable
OUT_JSON = ROOT / "verification" / "REPRODUCTION_RECORD.json"
OUT_MD = ROOT / "verification" / "REPRODUCTION_RECORD.md"


#: Step output goes to files, not pipes, so it survives the child being killed.
#: That mattered here: several runs of this script were reported as dead while
#: their pytest children were still alive and competing with the next attempt.
#: A log on disk shows what actually happened; a pipe held by a dead parent does
#: not.
LOGDIR = ROOT / "verification" / "reproduction-logs"


class Result:
    """Enough of CompletedProcess for the callers, plus where it came from."""

    def __init__(self, returncode: int, stdout: str, resumed: bool = False):
        self.returncode, self.stdout, self.resumed = returncode, stdout, resumed


def complete(tag: str, text: str) -> bool:
    """Whether a step log records a step that ran to completion.

    A long step on this machine is sometimes killed part-way, leaving a log that
    looks like progress. Resuming from one of those would report a partial run
    as a finished one, so each tag says what its own completion marker is.
    """
    if tag.startswith("tests-"):
        return bool(re.search(r"\d+ (passed|failed|error)", text)) or \
            bool(re.search(r"\[100%\]", text))
    if tag == "build":
        return "packaged ->" in text
    if tag == "gates":
        return "manuscript gates:" in text
    return bool(text.strip())


def run(cmd: list[str], cwd: Path = ROOT, timeout: int = 3600, tag: str = "step",
        resume: bool = False):
    LOGDIR.mkdir(parents=True, exist_ok=True)
    log = LOGDIR / f"{tag}.log"
    if resume and log.exists():
        text = log.read_text(errors="replace")
        if complete(tag, text):
            return Result(0, text, resumed=True), 0.0
    t0 = time.time()
    with log.open("w") as fh:
        proc = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT,
                              text=True, timeout=timeout)
    return Result(proc.returncode, log.read_text(errors="replace")), \
        round(time.time() - t0, 1)


def count_tests(output: str) -> int | None:
    """Passed count from pytest, whether or not it printed a summary line.

    This configuration prints the progress dots and no summary line, so the
    fallback is not an edge case -- it is the normal path, and it has to be
    right. An earlier version matched only the first run of dots on each line
    and reported 1 for a 72-test suite.
    """
    m = re.search(r"(\d+) passed", output)
    if m:
        return int(m.group(1))
    # One dot per passing test. Strip the progress markers first so their
    # punctuation is not counted, then count what remains.
    body = re.sub(r"\[\s*\d+%\]", "", output)
    body = re.sub(r"^\S+\.py\b.*$", "", body, flags=re.MULTILINE)
    dots = body.count(".")
    return dots or None


def flush(record: dict) -> None:
    """Write the record after every step.

    A killed run must leave the steps it finished, not an empty file. This
    project has now been bitten three times by a summary written only at the
    end of a loop whose body wrote everything else.
    """
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(record, indent=1) + "\n")


def step_tests(record: dict, resume: bool) -> bool:
    ok = True
    for label, cwd in (("tensor test suite", ROOT),
                       ("bridge test suite", ROOT / "spinor_trace_bridge")):
        tag = "tests-" + ("bridge" if "bridge" in label else "tensor")
        proc, secs = run([PY, "-m", "pytest", "-q", "-p", "no:cacheprovider"],
                         cwd=cwd, tag=tag, resume=resume)
        n = count_tests(proc.stdout)
        passed = proc.returncode == 0
        record["steps"].append({
            "step": label, "passed": passed, "tests": n, "seconds": secs,
            # A killed process exits nonzero with no summary; say so rather than
            # reporting the dots it managed to print as a result.
            "resumed_from_log": proc.resumed,
            "note": None if passed else "nonzero exit; suite did not complete",
        })
        ok = ok and passed
    flush(record)
    return ok


def step_regenerate(record: dict, resume: bool) -> bool:
    ok = True
    for label, script in (("numbers", "manuscript/scripts/make_numbers.py"),
                          ("tables", "manuscript/scripts/make_tables.py"),
                          ("figures", "manuscript/scripts/make_figures.py")):
        proc, secs = run([PY, str(ROOT / script)], tag=f"make-{label}",
                         resume=False)
        passed = proc.returncode == 0
        detail = None
        if label == "numbers":
            m = re.search(r"\((\d+) lines\)", proc.stdout)
            macros = len(re.findall(r"\\newcommand",
                                    (ROOT / "manuscript/generated/numbers.tex").read_text()))
            missing = len(re.findall(r"ARTIFACTMISSING",
                                     (ROOT / "manuscript/generated/numbers.tex").read_text()))
            detail = {"macros": macros, "artifacts_missing": missing}
            if missing:
                passed = False
        record["steps"].append({"step": f"regenerate {label}", "passed": passed,
                                "seconds": secs, "detail": detail})
        ok = ok and passed
    flush(record)
    return ok


def step_build(record: dict, resume: bool) -> bool:
    proc, secs = run([PY, str(ROOT / "scripts/build_submission_package.py")],
                     tag="build", resume=resume)
    diag = None
    m = re.search(r"isolated build: (\{.*\})", proc.stdout)
    if m:
        diag = json.loads(m.group(1))
    passed = proc.returncode == 0 and bool(diag) and diag.get("errors") == 0 \
        and diag.get("undefined_citations") == 0 \
        and diag.get("undefined_references") == 0
    record["steps"].append({"step": "isolated manuscript build", "passed": passed,
                            "seconds": secs, "detail": diag,
                            "resumed_from_log": proc.resumed})
    flush(record)
    return passed


def step_policy(record: dict) -> bool:
    """Release-policy facts must match the repository they describe."""
    proc, secs = run([PY, str(ROOT / "scripts/check_release_policy.py")],
                     tag="policy", resume=False)
    passed = proc.returncode == 0
    record["steps"].append({
        "step": "release policy check", "passed": passed, "seconds": secs,
        "note": None if passed else proc.stdout.strip().splitlines()[-3:],
    })
    flush(record)
    return passed


def step_gates(record: dict, resume: bool) -> bool:
    proc, secs = run([PY, str(ROOT / "manuscript/scripts/check_manuscript.py")],
                     tag="gates", resume=False)
    m = re.search(r"manuscript gates: (\d+) checks", proc.stdout)
    passed = proc.returncode == 0
    record["steps"].append({
        "step": "manuscript gates", "passed": passed, "seconds": secs,
        "detail": {"checks": int(m.group(1)) if m else None},
        "note": None if passed else proc.stdout.strip().splitlines()[-1:],
    })
    flush(record)
    return passed


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-tests", action="store_true")
    ap.add_argument("--skip-build", action="store_true",
                    help="skip the LaTeX build, for a machine with no TeX")
    ap.add_argument("--resume", action="store_true",
                    help="accept a completed step log instead of re-running the "
                         "step. Only logs that record a COMPLETED step are "
                         "accepted; a partial one is re-run.")
    args = ap.parse_args()

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--short"], cwd=ROOT,
                                capture_output=True, text=True).stdout.strip())

    record = {
        "schema": 1,
        "generated_by": "scripts/reproduce_all.py",
        "commit": head,
        "working_tree_clean": not dirty,
        "when": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            # Deliberately NOT the interpreter's absolute path: this file ships
            # in the release candidate, whose path scan rejects a home
            # directory, and it caught this.
            "interpreter": "repository virtualenv" if PY != sys.executable
                           else "ambient interpreter",
        },
        "not_regenerated": {
            "results/rank81/certificate.json": "needs the third-party archive",
            "results/rank81/minor81_certificate.json": "needs the third-party archive",
            "verification/SPINOR_JACOBIAN_RUNS.json": "needs the third-party archive",
        },
        "steps": [],
    }

    ok = True
    if args.skip_tests:
        record["steps"].append({"step": "test suites", "passed": None,
                                "note": "skipped by --skip-tests"})
    else:
        ok &= step_tests(record, args.resume)
    ok &= step_regenerate(record, args.resume)
    if args.skip_build:
        record["steps"].append({"step": "isolated manuscript build",
                                "passed": None, "note": "skipped by --skip-build"})
    else:
        ok &= step_build(record, args.resume)
    ok &= step_policy(record)
    ok &= step_gates(record, args.resume)

    record["all_steps_passed"] = ok
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(record, indent=1) + "\n")

    rows = ["# Reproduction record", "",
            "Generated by `scripts/reproduce_all.py`. Do not edit by hand --",
            "this file was maintained manually twice and went stale both times.",
            "",
            f"Commit `{head}`, working tree "
            f"{'clean' if not dirty else '**dirty**'}, {record['when']}.", "",
            "| step | result | detail | seconds |", "|---|---|---|---:|"]
    for s in record["steps"]:
        mark = {True: "**pass**", False: "**FAIL**", None: "skipped"}[s.get("passed")]
        detail = s.get("detail") or s.get("note") or ""
        if s.get("tests") is not None:
            detail = f"{s['tests']} tests"
        if s.get("resumed_from_log"):
            # Without this a resumed step reads as a fresh run that took 0.0 s,
            # which is precisely the impression this file must not give.
            mark = "pass (resumed)" if s.get("passed") else mark
            detail = f"{detail} — read from a completed step log, not re-run"
        rows.append(f"| {s['step']} | {mark} | {detail} | {s.get('seconds', '')} |")
    rows += ["", "## Not regenerated here", ""]
    for path, why in record["not_regenerated"].items():
        rows.append(f"- `{path}` --- {why}")
    rows += ["",
             "These certificates are committed, so every number in both",
             "manuscripts is present and both build without the archive. Only",
             "their *regeneration* needs it. See",
             "`release_candidate/spinor-archive/MANIFEST.json` and",
             "`scripts/verify_archive.py`.", ""]
    OUT_MD.write_text("\n".join(rows))

    print(f"wrote {OUT_JSON.relative_to(ROOT)} and {OUT_MD.relative_to(ROOT)}")
    print("ALL STEPS PASSED" if ok else "SOME STEPS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
