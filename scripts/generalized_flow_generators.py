#!/usr/bin/env python3
"""Run genuine f(tau,S) evaluations using an explicit experiment config."""
from pathlib import Path
import argparse,json
from chiral4form.pipeline import run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path('configs/degree8_completion.json'))
    parser.add_argument('--output',type=Path,default=Path('runs/degree8-generalized'))
    parser.add_argument('--source',type=Path)
    args=parser.parse_args()
    config=json.loads(args.config.read_text())
    if not config.get('extras'):
        parser.error('configuration must declare actual extra scalar generators')
    run(config,Path.cwd(),args.output,args.source)

if __name__=='__main__':main()
