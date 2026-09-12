#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path

from chiral4form.predecessor_import import import_predecessor_certificates

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--existing", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument(
        "--obstruction",
        type=Path,
        default=Path("verification/degree10/obstruction_certificate.json"),
    )
    p.add_argument(
        "--degree10-model",
        type=Path,
        default=Path("verification/degree10/rational_model.json"),
    )
    args = p.parse_args()
    result = import_predecessor_certificates(
        args.source,
        args.existing,
        args.output,
        args.obstruction,
        args.degree10_model,
    )
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
