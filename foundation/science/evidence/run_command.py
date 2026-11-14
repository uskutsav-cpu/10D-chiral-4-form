#!/usr/bin/env python3
"""Capture one independent-audit command without reusing any prior output."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import signal
import time

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True)
    parser.add_argument("--cwd", type=Path, default=ROOT / "frozen-checkout")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    directory = ROOT / "commands" / args.name
    directory.mkdir(parents=True, exist_ok=False)
    overrides = {"PYTHONNOUSERSITE": "1", "PYTHONPYCACHEPREFIX": str(ROOT / "pycache"),
                 "XDG_CACHE_HOME": str(ROOT / "cache"), "PIP_NO_CACHE_DIR": "1",
                 "PIP_CONFIG_FILE": "/dev/null", "OMP_NUM_THREADS": "1",
                 "OPENBLAS_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
    record = {"argv": command, "cwd": str(args.cwd.resolve()), "environment_overrides": overrides,
              "started_unix": time.time(), "status": "running"}
    metadata = directory / "command.json"
    metadata.write_text(json.dumps(record, indent=2) + "\n")
    started = time.monotonic()
    with (directory / "stdout.log").open("wb") as stdout, (directory / "stderr.log").open("wb") as stderr:
        process = subprocess.Popen(command, cwd=args.cwd, env={**os.environ, **overrides}, stdout=stdout, stderr=stderr)
        record["pid"] = process.pid
        metadata.write_text(json.dumps(record, indent=2) + "\n")
        def interrupted(signum, _frame):
            record["interrupted_signal"] = signum
            if process.poll() is None:
                process.terminate()
        signal.signal(signal.SIGTERM, interrupted)
        signal.signal(signal.SIGINT, interrupted)
        record["exit_code"] = process.wait()
    record.update({"seconds": time.monotonic() - started, "finished_unix": time.time(), "status": "completed"})
    record["logs"] = {name: {"sha256": hashlib.sha256((directory / name).read_bytes()).hexdigest(),
                              "bytes": (directory / name).stat().st_size} for name in ("stdout.log", "stderr.log")}
    metadata.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2), flush=True)
    return record["exit_code"]


if __name__ == "__main__":
    raise SystemExit(main())
