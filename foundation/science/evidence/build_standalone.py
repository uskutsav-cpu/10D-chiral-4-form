from pathlib import Path
import ast,hashlib,json,random,secrets
REPO=Path(__file__).resolve().parent.parent/'independent-audit-20260907-02/repair-checkout'
base=REPO/'results/audit/math-review'
selected={'independent_conventions.py':['parity','metric_product','star_image','compact_form','slot_records'],'independent_oracle.py':['dense_form','terms','plan','modular_dot'],'efficient_oracle.py':['PointEngine'],'independent_reverse.py':['sparse_coordinate_pullback','reverse_gradient'],'bareiss_minors.py':['bareiss']}
header='''#!/usr/bin/env python3
"""Recompute a rank81 theorem from graphs, coordinates and exact points only.

No stored values/Jacobians, candidate search, Hilbert counts, literature maps,
repository imports, files of intermediate arrays, or resume logic are used.
NumPy/opt_einsum provide bounded exact contractions and path selection.
"""
import argparse
from collections import Counter, OrderedDict
import hashlib
import itertools
import json
from pathlib import Path
import time
import numpy as np
import opt_einsum as oe

TUPLES=list(itertools.combinations(range(10),5))
POSITIONS={t:k for k,t in enumerate(TUPLES)}
COORDINATES=[k for k,t in enumerate(TUPLES) if 0 in t]
DIM=10
'''
parts=[header];provenance=[]
for filename,names in selected.items():
 text=(base/filename).read_text();tree=ast.parse(text)
 for name in names:
  node=next(n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name)
  source=ast.get_source_segment(text,node)
  parts.append(source)
  provenance.append({'source':str((base/filename).relative_to(REPO)),'symbol':name,'source_sha256':hashlib.sha256(source.encode()).hexdigest()})
parts.append('PULLBACK_INDICES,PULLBACK_SIGNS=sparse_coordinate_pullback()')
parts.append((Path(__file__).parent/'standalone_main.py').read_text())
(REPO/'science/rank81/verify.py').write_text('\n\n'.join(parts)+'\n')
(REPO/'science/rank81/implementation_provenance.json').write_text(json.dumps({'reused':'Previously independently validated program bodies, copied exactly; no prior numerical arrays reused.','symbols':provenance},indent=2)+'\n')
inputpath=REPO/'science/rank81/input.json'
assert not inputpath.exists(),'Never overwrite the fresh frozen point'
seed=secrets.randbits(64);rng=random.Random(seed);p=50021
old=json.loads((REPO/'results/rank81_basis.json').read_text())
data={'schema':1,'convention':{'metric':[-1]+[1]*9,'coordinate_tuples':[list(TUPLES) for TUPLES in __import__('itertools').combinations(range(10),5) if 0 in TUPLES],'hodge':'output-first epsilon; epsilon_0123456789=+1; F=sum A_I(e_I+star e_I)'},'graphs':[{k:i[k] for k in ('id','degree','graph')} for i in old['invariants']],'points':[{'prime':p,'seed':seed,'coordinates':[rng.randrange(p) for _ in range(126)]}]}
inputpath.write_text(json.dumps(data,indent=2)+'\n');print('Fresh standalone point',p,seed)
