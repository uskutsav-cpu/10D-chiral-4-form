#!/usr/bin/env python3
"""Require an exact, clean predecessor commit and every required source file."""
import argparse,json
from pathlib import Path
from chiral4form.source import acquire_source
from chiral4form.provenance import atomic_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_repo',type=Path)
    parser.add_argument('--out',type=Path,default=Path('data/imported/source_manifest.json'))
    args=parser.parse_args()
    _,manifest=acquire_source(args.source_repo)
    atomic_json(args.out,manifest)
    print(args.out)

if __name__=='__main__':main()
