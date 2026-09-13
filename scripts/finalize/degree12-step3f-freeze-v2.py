#!/usr/bin/env python3
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--discovery',type=Path,default=Path('verification/degree12/generalized_discovery.json'));ap.add_argument('--completion',type=Path,default=Path('verification/degree12/generalized_completion_theorem.json'));ap.add_argument('--minimality',type=Path,default=Path('verification/degree12/generalized_removal_minimality.json'));ap.add_argument('--output',type=Path,default=Path('verification/degree12/generalized_completion_manifest.json'));a=ap.parse_args()
    d=json.loads(a.discovery.read_text());c=json.loads(a.completion.read_text());m=json.loads(a.minimality.read_text())
    if c.get('global_generalized_completion') is not True or c.get('full_polynomial_interaction_space_reachable_through_degree12') is not True: raise SystemExit('generalized completion theorem missing')
    if d['selected_extras']!=c['selected_extras']: raise SystemExit('selected sets differ')
    out={'schema':2,'status':'frozen_global_generalized_completion_through_degree12','selected_extras':c['selected_extras'],'selected_count':c['selected_count'],'characteristic_zero_global_completion':True,'full_polynomial_reachability_through_degree12':True,
         'catalogue_optimal_count_on_discovery_primes':d['catalogue_optimal_count_mod_each_prime'],
         'fresh_prime_removal_audit':{'all_selected_fresh_prime_nonredundant':m['all_selected_extras_fresh_prime_nonredundant'],'exact_count_nonredundant_count':len(m['exact_count_nonredundant_extras']),'fresh_prime_nonredundant_count':len(m['fresh_prime_nonredundant_extras'])},
         'paper_safe_claims':['The selected genuine S generators globally complete the finite polynomial generalized-flow system through degree 12.','The proof uses full-rank global constant translations at every homogeneous layer.','The selected primitive catalogue has the recorded removal-nonredundancy status on independent fresh-prime models.'],
         'do_not_claim':['global cardinality minimality over arbitrary scalar combinations','minimality for every possible nonlinear dynamic generalized completion','all-orders completion','nonanalytic completion'],
         'files_sha256':{str(p):sha(p) for p in (a.discovery,a.completion,a.minimality)}}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print('PASS: STEP 3F GENERALIZED COMPLETION FROZEN');print('selected extras:',out['selected_count']);print('WROTE:',a.output)
if __name__=='__main__':main()
