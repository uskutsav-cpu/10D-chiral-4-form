"""Command line entry points; errors never turn into success receipts."""
from __future__ import annotations
import argparse
import json
import sys
import subprocess
from pathlib import Path
from .provenance import atomic_json


def _json(path):
    return json.loads(Path(path).read_text())


def main(argv=None):
    parser=argparse.ArgumentParser(prog='chiral4form')
    parser.add_argument('--root',type=Path,default=Path.cwd(),help='repository root for configuration/fixture paths')
    sub=parser.add_subparsers(dest='command',required=True)
    for cmd in ('run','plan','worker'):
        p=sub.add_parser(cmd);p.add_argument('--config',type=Path,default=Path('configs/smoke.json'))
        p.add_argument('--source',type=Path);p.add_argument('--output',type=Path,default=Path('results/smoke'))
        if cmd=='worker':p.add_argument('--prime-index',type=int,required=True)
    p=sub.add_parser('source');p.add_argument('--checkout',type=Path,default=Path('.cache/upstream'))
    p.add_argument('--output',type=Path,default=Path('data/imported'));p.add_argument('--fetch',action='store_true')
    p=sub.add_parser('status');p.add_argument('--output',type=Path,default=Path('results/smoke'))
    for cmd in ('sextic','degree8','completion8','modmax'):
        p=sub.add_parser(cmd);p.add_argument('--output',type=Path)
    p=sub.add_parser('analyze');p.add_argument('model',type=Path);p.add_argument('--depth',type=int,default=3)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('lift');p.add_argument('models',nargs='+',type=Path);p.add_argument('--bound',type=int,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('complete');p.add_argument('model',type=Path);p.add_argument('--candidates',nargs='+',required=True)
    p.add_argument('--target-rank',type=int,required=True);p.add_argument('--metric',choices=['lie_rank','linear_hull'],default='lie_rank')
    p.add_argument('--depth',type=int,default=3);p.add_argument('--max-subsets',type=int,default=1024)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('prove-map8');p.add_argument('--output',type=Path,default=Path('results/degree8-map-proof'));p.add_argument('--prime-index',type=int)
    p=sub.add_parser('verify-map8');p.add_argument('certificate',type=Path);p.add_argument('--reevaluate-tensors',action='store_true')
    args=parser.parse_args(argv)
    try:
        if args.command in ('run','plan','worker'):
            from .pipeline import run,execution_plan
            config=_json(args.root/args.config)
            if args.command=='worker':
                from .worker import run_prime
                run_prime(config,args.prime_index,args.root,args.output,args.source);return 0
            if args.command=='plan':
                result=execution_plan(config,args.root,args.source)
            else:
                run(config,args.root,args.output,args.source);return 0
        elif args.command=='source':
            from .source import export_source
            result=export_source(args.checkout,args.output,args.fetch)
        elif args.command=='status':
            result=_json(args.output/'status.json')
        elif args.command=='sextic':
            from .sextic import sextic_report
            result=sextic_report()
        elif args.command=='degree8':
            from .degree8_orbit import report
            result=report()
        elif args.command=='completion8':
            from .degree8_completion import report
            result=report()
        elif args.command=='modmax':
            from .localization import modmax_symbolic_report
            result=modmax_symbolic_report()
        elif args.command=='analyze':
            from .fitting import fields_from_json
            from .reachability import linear_invariant_hull,lie_rank
            ids,fields=fields_from_json(_json(args.model))
            result={'coordinate_ids':list(ids),'linear_hull':linear_invariant_hull(list(fields.values())),
                    'lie':lie_rank(list(fields.values()),max_depth=args.depth)}
        elif args.command=='lift':
            from .lifting import lift_models
            result=lift_models([_json(p) for p in args.models],args.bound)
        elif args.command=='complete':
            from .completion import search
            result=search(_json(args.model),args.candidates,args.target_rank,args.metric,args.depth,args.max_subsets)
        elif args.command=='prove-map8':
            from .map_proof import run as prove_map8
            certificate=prove_map8(args.output,args.root,prime_index=args.prime_index)
            result={'certificate_hash':certificate['certificate_hash'],'status':certificate['status'],'output':str(args.output)}
        elif args.command=='verify-map8':
            from .map_proof import verify
            valid=verify(_json(args.certificate),reevaluate_tensors=args.reevaluate_tensors)
            if not valid:raise ValueError('degree-eight map certificate verification failed')
            result={'verified':True,'tensor_values_reevaluated':args.reevaluate_tensors,
                    'scope':'bounded arithmetic and stated external inputs; not all orders'}
        output=getattr(args,'output',None)
        if output and args.command not in ('source','status','plan','prove-map8'):
            atomic_json(output,result);print(output)
        else:
            print(json.dumps(result,indent=2))
        return 0
    except (ValueError,TypeError,RuntimeError,OSError,MemoryError,OverflowError,subprocess.SubprocessError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        return 2

if __name__=='__main__':
    raise SystemExit(main())
