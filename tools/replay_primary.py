"""Replay the primary results without requiring the two omitted historical locks."""
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
import argparse
import sys
from datetime import datetime, timezone
from prepare_legacy_layout import ROOT as REPO, COMPATIBILITY, prepare
import verify_release

parser = argparse.ArgumentParser(description="Replay primary panels using retained hashes and portable metadata.")
parser.add_argument("--allow-version-difference", action="store_true", help="Test numerical replay in another environment and report that difference explicitly")
args = parser.parse_args()
verify_release.main()
verify_release.verify_followups()
prepare()
sys.path.insert(0, str(COMPATIBILITY))
import balanced_assay_initialization_time as study
import test_balanced_assay_initialization_time as boundary

ROOT=REPO/'reproducibility'
METADATA=json.loads((REPO/'results/primary/replay_metadata.json').read_text(encoding='utf-8'))
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
    for rel, expected in METADATA['retained_dependency_sha256'].items():
        path=REPO/rel
        assert path.is_file() and sha(path)==expected,rel
        checked[rel]=expected
    for stage in ['balanced_assay_gate_control','balanced_assay_ep_adapter',
                  'balanced_assay_timecheck','balanced_assay_tradeoff','balanced_assay_pilot']:
        locks(study.ROOT/stage/'lock.json')
    assert METADATA['primary_config']==study.CONFIG
    actual={'python':platform.python_version(),'numpy':np.__version__,'openpyxl':openpyxl.__version__}
    if not args.allow_version_difference:
        assert actual==METADATA['recorded_versions'], 'Use the recorded environment, or explicitly pass --allow-version-difference for a portability check'

def verify_retained_seal():
    seal=read(OUT/'sealed_panels.json')
    expected=set(study.expected_files())|{f'instance_{s}_{seed}.json' for s in ['YDIP','YBIOGRID'] for seed in study.CONFIG['seeds']}|{'instance_SpaceI.json'}
    assert set(seal['files'])==expected and seal['panels']==360
    original_digest=METADATA['omitted_locks_sha256']['reproducibility/primary_archive/ppi_data_audit/balanced_assay_initialization_time_20261024/lock.json']
    assert seal['lock_sha256']==original_digest
    assert all(sha(OUT/name)==value for name,value in seal['files'].items())
    return [read(OUT/name) for name in study.expected_files()]

def main():
    verify_versions_and_locks()
    before={f.name:sha(f) for f in [OUT/'scores.json',OUT/'sealed_panels.json']}
    saved=read(OUT/'scores.json'); captured={}
    def capture(path,data):captured[path.name]=data
    # Replay all 360 capacities, coverage, input/split hashes, timing, scores/pairs.
    # Redirect only writes and platform-specific lock entry points, never the solver.
    with patch.object(study,'lock',verify_versions_and_locks),patch.object(study,'verify_seal',verify_retained_seal),patch.object(study,'write',capture),contextlib.redirect_stdout(io.StringIO()):study.score()
    assert captured['scores.json']==saved,'Replay differs from sealed results'
    rows=saved['rows'];pairs=saved['paired']; panels=verify_retained_seal()
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
    assert before=={f.name:sha(f) for f in [OUT/'scores.json',OUT/'sealed_panels.json']}
    result={'replayed_utc':datetime.now(timezone.utc).isoformat(),'checks':{'retained_dependency_hashes':True,'omitted_historical_locks':2,'runtime_versions_match':({'python':platform.python_version(),'numpy':np.__version__,'openpyxl':openpyxl.__version__}==METADATA['recorded_versions']),'sealed_panels':360,'calibration_files':41,'independent_evaluation_counts':360,'paired_statistics':180,'summary_groups':9,'boundary_tests':6,'score_replay_identical':True,'primary_timing_valid':sum(x['timing_valid'] for x in rows if x['multiplier']==5),'extended_timing_valid':sum(x['timing_valid'] for x in rows if x['multiplier']==20),'all_timing_valid':sum(x['timing_valid'] for x in rows),'zero_relative_denominators':sum(x['a']==0 for x in pairs),'frozen_files_unchanged':True},'implementation_details':{'shared_solver':'balanced_assay_timecheck.multistart(initial,floor,remaining,seed,True)','restart_anchors':'start 0: own initial; subsequent odd: best; subsequent even: own initial','search_deadline_reserve_seconds':0.01,'floor_source':'shared scalar from calibration, formal paths regenerate panels','timing_scope':'pipeline entry through joint solver return; graph building, calibration, final output verification and I/O excluded','relative_statistic':'median of within-pair (B-A)/A; zero denominators unavailable'},'versions':{'python':platform.python_version(),'numpy':np.__version__,'openpyxl':openpyxl.__version__},'inputs_and_frozen_outputs':dict(sorted(checked.items())),'score_sha256':before['scores.json']}
    DEST.mkdir(parents=True,exist_ok=True)
    (DEST/'primary_comparison_replay.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result['checks'],indent=2))

if __name__=='__main__':main()
