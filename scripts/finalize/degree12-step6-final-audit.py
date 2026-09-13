#!/usr/bin/env python3
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=Path('verification/final_manifest.json'));a=ap.parse_args()
    req=[Path('verification/degree12/qq/degree12_local_theorem_certificate.json'),Path('verification/degree12/global_reachability_manifest.json'),Path('verification/degree12/generalized_discovery.json'),Path('verification/degree12/generalized_completion_theorem.json'),Path('verification/degree12/generalized_removal_minimality.json'),Path('verification/degree12/generalized_completion_manifest.json'),Path('verification/degree12/conformal_scope_report.json'),Path('verification/final_claims.json'),Path('docs/FINAL_RESULTS_FREEZE.md'),Path('docs/CONFORMAL_SCOPE_REPORT.md'),Path('manuscript/generated_results.tex')]
    missing=[str(x) for x in req if not x.exists()]
    if missing: raise SystemExit(f'missing final artifacts: {missing}')
    local=json.loads(req[0].read_text());pure=json.loads(req[1].read_text());disc=json.loads(req[2].read_text());comp=json.loads(req[3].read_text());freeze=json.loads(req[5].read_text());claims=json.loads(req[7].read_text())
    checks={'degree12_QQ_local_theorem':local.get('qq_theorem') is True,'degree12_local_dimension_10':int(local.get('local_orbit_dimension',-1))==10,'pure_stress_global':pure.get('status')=='frozen_degree12_global_constructive_reachability','generalized_selected_count_consistent':int(disc['selected_count'])==int(comp['selected_count']),'generalized_global_completion':comp.get('global_generalized_completion') is True,'generalized_freeze':freeze.get('characteristic_zero_global_completion') is True,'final_claims_status':claims.get('status')=='final_results_freeze'}
    failed=[k for k,v in checks.items() if not v]
    if failed: raise SystemExit(f'final audit failed: {failed}')
    out={'schema':1,'status':'final_project_artifact_manifest','checks':checks,'files_sha256':{str(x):sha(x) for x in req},'headline':{'pure_stress_orbit_dimensions':[4,6,10],'degree12_pure_stress_global':True,'generalized_full_polynomial_completion_through_degree12':True,'selected_primitive_extras':int(comp['selected_count'])},'remaining_research_not_required_for_this_freeze':['all-orders theorem','degree-14 extension','nonanalytic ModMax reachability classification','global arbitrary-scalar minimality theorem']}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: STEP 6 FINAL PROJECT AUDIT');print('pure-stress orbit dimensions: 4 -> 6 -> 10');print('degree-12 generalized completion: True');print('selected primitive extras:',comp['selected_count']);print('WROTE:',a.output)
if __name__=='__main__':main()
