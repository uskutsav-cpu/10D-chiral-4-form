#!/usr/bin/env python3
"""Research driver placeholder for *genuine* generalized-flow generators.

This command intentionally refuses to equate seed augmentation with generator
augmentation.  Physics-specific forcing rows must be supplied by the upstream
stress-tensor engine before a completion claim can be produced.
"""

from __future__ import annotations

import argparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--declare-generator", action="append", default=[])
    parser.parse_args()
    raise SystemExit(
        "NOT YET IMPLEMENTED: generator augmentation requires metric/stress variation "
        "of each declared S_i. Seed augmentation is not accepted as a substitute."
    )


if __name__ == "__main__":
    main()
