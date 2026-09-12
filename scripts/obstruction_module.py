#!/usr/bin/env python3
"""Compute an exact dual obstruction certificate from a reachable-row matrix.

Input JSON:
{
  "prime": 32749,
  "ambient_dim": 4,
  "reachable_rows": [[...], ...]
}
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from chiral4form.certificates import build_quotient_certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    cert = build_quotient_certificate(
        payload["reachable_rows"], int(payload["ambient_dim"]), int(payload["prime"])
    )
    out = {
        "prime": cert.prime,
        "ambient_dim": cert.ambient_dim,
        "reachable_rank": cert.reachable_rank,
        "quotient_dim": cert.quotient_dim,
        "annihilators": [list(v) for v in cert.annihilators],
    }
    text = json.dumps(out, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
