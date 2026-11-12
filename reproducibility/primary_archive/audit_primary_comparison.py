"""Read-only replay of sealed primary outputs; never regenerate timed searches.

Run from any checkout root: python tmp/audit_primary_comparison.py
Only the new audit JSON is written; original locks, scores and panels are untouched.
"""
import contextlib
import hashlib
import io
import json
import platform
from pathlib import Path
from statistics import median
from unittest.mock import patch
import unittest
import numpy as np
import openpyxl
import balanced_assay_initialization_time as study
import test_balanced_assay_initialization_time as boundary

ROOT=Path(__file__).resolve().parent.parent
OUT=study.OUT
DEST=ROOT/'output/research'
checked={}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def resolve(record):
    # Original Windows absolute keys and tmp-relative keys have different bases.
    normalized=record.replace('\\','/')
    if '/tmp/' in normalized:return ROOT/'tmp'/normalized.split('/tmp/',1)[1]
    return ROOT/'tmp'/normalized

def locks(path):
    relative=path.relative_to(ROOT).as_posix()
    if relative in checked:return
    checked[relative]=sha(path)
    for record,expected in read(path)['files'].items():
        f=resolve(record)
        assert f.is_file() and sha(f)==expected,record
        checked[f.relative_to(ROOT).as_posix()]=expected
        if f.name=='lock.json':
            checked.pop(f.relative_to(ROOT).as_posix(),None)
            locks(f)

def verify_versions_and_locks():
    locks(OUT/'lock.json')
    # Some historical lock() calls invoke a predecessor without listing its
    # manifest in their own file map. Explicitly include the whole code chain.
    for stage in ['balanced_assay_gate_control','balanced_assay_ep_adapter',
                  'balanced_assay_timecheck','balanced_assay_tradeoff','balanced_assay_pilot']:
        locks(study.ROOT/stage/'lock.json')
    frozen=read(OUT/'lock.json')
    assert frozen['config']==study.CONFIG
    assert platform.python_version()==frozen['python']
    assert np.__version__==frozen['numpy']
    original=read(study.ROOT/'balanced_assay_pilot/lock.json')
    assert openpyxl.__version__==original['openpyxl']

def main():
    verify_versions_and_locks()
    before={f.name:sha(f) for f in [OUT/'scores.json',OUT/'sealed_panels.json',OUT/'lock.json']}
    saved=read(OUT/'scores.json'); captured={}
    def capture(path,data):captured[path.name]=data
    # Replay all 360 capacities, coverage, input/split hashes, timing, scores/pairs.
    # Redirect only writes and platform-specific lock entry points, never the solver.
    with patch.object(study,'lock',verify_versions_and_locks),patch.object(study,'write',capture),contextlib.redirect_stdout(io.StringIO()):study.score()
    assert captured['scores.json']==saved,'Replay differs from sealed results'
    rows=saved['rows'];pairs=saved['paired']; panels=study.verify_seal()
    assert len(rows)==360 and len(pairs)==180
    assert len({(x['source'],x['seed'],x['multiplier'],x['method']) for x in rows})==360
    for pair in pairs:
        assert pair['delta']==pair['b']-pair['a']
        assert pair['relative']==(pair['delta']/pair['a'] if pair['a'] else None)
    for s in saved['summary']:
        group=[x for x in pairs if x['source']==s['source'] and x['multiplier']==s['multiplier']]
        assert len(group)==20
        assert median(x['delta'] for x in group)==s['delta_median']
        assert median(x['relative'] for x in group if x['relative'] is not None)==s['relative_median']
    for row in rows:
        assert row['total_seconds']>=row['initialization_seconds']+row['search_seconds']
        assert row['search_allowance']<=max(0,row['allowed_seconds']-row['initialization_seconds'])+1e-9
        if not row['initialization_fallback']:assert row['train_cross']>=row['initialization_train']
    # Independently count each evaluation edge with Python set membership.
    yeast_names,sources=study.inputs()
    human=study.external.restricted(study.external.novel(study.external.edges(study.external.INPUTS['outcome']),[study.external.edges(study.external.INPUTS[k]) for k in study.external.CONFIG['exclusion']]),study.external.cohort())
    indexed={(r['source'],r['seed'],r['multiplier'],r['method']):r for r in saved['rows']}
    for s in study.CONFIG['sources']:
        for seed in study.CONFIG['seeds']:
            if s=='SpaceI':test=human
            else:
                e=sources[s];order=np.random.default_rng(seed).permutation(len(e));test=e[order[int(.8*len(e)):]]
            for row in [x for x in panels if x['source']==s and x['seed']==seed]:
                db,ad=set(row['db']),set(row['ad'])
                n=sum((u in db and v in ad) or (u in ad and v in db) for u,v in test.tolist())
                assert n==indexed[s,seed,row['multiplier'],row['method']]['test_cross']
    tests=unittest.TextTestRunner(stream=io.StringIO()).run(unittest.defaultTestLoader.loadTestsFromModule(boundary))
    assert tests.wasSuccessful() and tests.testsRun==6
    # Calibration and old-start provenance files are verified without retiming.
    for calpath in sorted(OUT.glob('instance_*.json')):
        cal=read(calpath)
        old=study.external.OUT/'generation.json' if cal['source']=='SpaceI' else study.ROOT/f"balanced_assay_timecheck/{cal['source']}_{cal['split_seed']}.json"
        assert sha(old)==cal['original_file_sha256']
        checked[old.relative_to(ROOT).as_posix()]=sha(old)
    for path in [OUT/'sealed_panels.json',OUT/'scores.json',*OUT.glob('panel_*.json'),*OUT.glob('instance_*.json')]:checked[path.relative_to(ROOT).as_posix()]=sha(path)
    assert before=={f.name:sha(f) for f in [OUT/'scores.json',OUT/'sealed_panels.json',OUT/'lock.json']}
    result={'date':'2026-10-25','checks':{'frozen_dependency_hashes':True,'runtime_versions_match':True,'sealed_panels':360,'calibration_files':41,'independent_evaluation_counts':360,'paired_statistics':180,'summary_groups':9,'boundary_tests':6,'score_replay_identical':True,'primary_timing_valid':sum(x['timing_valid'] for x in rows if x['multiplier']==5),'extended_timing_valid':sum(x['timing_valid'] for x in rows if x['multiplier']==20),'all_timing_valid':sum(x['timing_valid'] for x in rows),'zero_relative_denominators':sum(x['a']==0 for x in pairs),'frozen_files_unchanged':True},'implementation_details':{'shared_solver':'balanced_assay_timecheck.multistart(initial,floor,remaining,seed,True)','restart_anchors':'start 0: own initial; subsequent odd: best; subsequent even: own initial','search_deadline_reserve_seconds':0.01,'floor_source':'shared scalar from calibration, formal paths regenerate panels','timing_scope':'pipeline entry through joint solver return; graph building, calibration, final output verification and I/O excluded','relative_statistic':'median of within-pair (B-A)/A; zero denominators unavailable'},'versions':{'python':platform.python_version(),'numpy':np.__version__,'openpyxl':openpyxl.__version__},'inputs_and_frozen_outputs':dict(sorted(checked.items())),'score_sha256':before['scores.json']}
    DEST.mkdir(exist_ok=True)
    (DEST/'primary_comparison_audit_2026-10-25.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result['checks'],indent=2))

if __name__=='__main__':main()
