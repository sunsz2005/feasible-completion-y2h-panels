"""Frozen Space-I mapping-cohort retrospective evaluation. No parameter tuning."""
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from balanced_assay_pilot import State, coverage_order, matrix, partition, sha
from balanced_assay_tradeoff import feasible_role_greedy
from balanced_assay_timecheck import multistart, role_blind
from balanced_assay_ep_adapter import adapter
from balanced_assay_gate_control import post_pipeline, lock as previous_lock

ROOT = Path(__file__).resolve().parent / 'ppi_data_audit'
OUT = ROOT / 'balanced_assay_external_20261024'
INPUTS = {
    'candidates': ROOT/'balanced_assay_space1_20261024/historical_mapping_candidates.json',
    'train': ROOT/'balanced_assay_source_audit/raw/HI-II-14.tsv',
    'outcome': ROOT/'HuRI.tsv',
    'earlier': ROOT/'balanced_assay_source_audit/raw/H-I-05.tsv',
    'calibration': ROOT/'balanced_assay_source_audit/raw/Test_space_screens-19.tsv',
}
CONFIG = dict(budget=200, db=100, ad=100, guard=.98, q=.8,
              seeds=list(range(20261020,20261040)),
              deterministic=['coverage_partition','feasible_role_greedy','role_blind'],
              timed=['multistart_roles','multistart_joint','ep_adapter_repaired_swap','after_role_search'],
              primary_comparator='multistart_joint',
              exclusion=['train','earlier','calibration'],
              timing_tolerance_fraction=.05, timing_tolerance_seconds=.02)


def write(path, value):
    pending = path.with_suffix('.pending')
    pending.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8')
    pending.replace(path)


def edges(path):
    result=set()
    for line in path.read_text(encoding='utf-8-sig').splitlines():
        parts=line.split()
        if not parts or parts[0].startswith('#'):
            continue
        if len(parts)!=2 or not all(x.startswith('ENSG') for x in parts):
            raise ValueError(f'Unexpected edge row in {path}')
        u,v=(x.split('.')[0] for x in parts)
        if u!=v:
            result.add(tuple(sorted((u,v))))
    return result


def cohort():
    rows=json.loads(INPUTS['candidates'].read_text(encoding='utf-8'))
    names=sorted(r['ensembl_gene_id'] for r in rows)
    assert len(names)==6277 and len(set(names))==len(names)
    assert len(set(r['orf_id'] for r in rows))==len(names)
    assert all(r['historical_collection_membership'] and r['db_pass'] is None and r['ad_pass'] is None for r in rows)
    return names


def restricted(values,names):
    ids={name:i for i,name in enumerate(names)}
    return np.array(sorted((ids[u],ids[v]) for u,v in values if u in ids and v in ids),dtype=np.int32).reshape(-1,2)


def lock(create=False):
    previous_lock()
    payload=dict(config=CONFIG,files={str(p.resolve()):sha(p) for p in [Path(__file__),OUT/'protocol.md',*INPUTS.values()]})
    target=OUT/'lock.json'
    if create:
        assert not target.exists(), 'Protocol already frozen'
        write(target,payload)
    else:
        assert json.loads(target.read_text(encoding='utf-8'))==payload, 'Frozen input/protocol changed'
    return payload


def panel_record(method,seed,s,elapsed,meta=None):
    assert np.count_nonzero(s.roles==1)==100 and np.count_nonzero(s.roles==2)==100
    return dict(method=method,seed=seed,db=np.flatnonzero(s.roles==1).tolist(),
                ad=np.flatnonzero(s.roles==2).tolist(),train_cross=int(s.score),
                coverage=float(s.coverage()),search_seconds=elapsed,meta=meta)


def generate():
    lock()
    assert not (OUT/'sealed_panels.json').exists(), 'Panels already sealed'
    assert not (OUT/'generation.json').exists(), 'Do not overwrite an existing experiment'
    names=cohort(); train=restricted(edges(INPUTS['train']),names)
    a=matrix(len(names),train)
    t=time.perf_counter(); selected=coverage_order(a,200); base=partition(a,selected)
    common=time.perf_counter()-t; reference=base.coverage(); floor=.98*reference
    t=time.perf_counter(); proposed,moves=feasible_role_greedy(a,selected,floor)
    budget=time.perf_counter()-t
    rows=[panel_record('coverage_partition',None,base,0.),
          panel_record('feasible_role_greedy',None,proposed,budget,dict(external_swaps=moves))]
    t=time.perf_counter(); blind=role_blind(a,selected,floor)
    rows.append(panel_record('role_blind',None,blind,time.perf_counter()-t))
    setup=dict(names=names,train_edges=len(train),train_nodes=int(np.count_nonzero(a.sum(axis=1))),
               reference_coverage=reference,coverage_floor=floor,common_seconds=common,
               allowed_search_seconds=budget,rows=rows)
    write(OUT/'generation.json',setup)
    print(f'Fixed train {len(train)} edges; search allowance {budget:.4f}s',flush=True)
    for i,seed in enumerate(CONFIG['seeds']):
        rows=[]
        methods=CONFIG['timed'][i%4:]+CONFIG['timed'][:i%4]
        for method in methods:
            if method.startswith('multistart'):
                s,meta=multistart(base,floor,budget,seed,method=='multistart_joint')
                elapsed=meta['search_seconds']
            elif method.startswith('ep_'):
                s,_,_,meta=adapter(a,train,base,floor,budget,seed)
                elapsed=meta['search_seconds']
            else:
                _,_,s,meta=post_pipeline(a,selected,base,floor,budget,'role')
                elapsed=meta['total_seconds']
            assert s.coverage()+1e-7>=floor
            rows.append(panel_record(method,seed,s,elapsed,meta))
        write(OUT/f'panels_{seed}.json',dict(rows=rows))
        print(f'Panels complete: {i+1}/20 seeds',flush=True)
    paths=[OUT/'generation.json',*[OUT/f'panels_{seed}.json' for seed in CONFIG['seeds']]]
    write(OUT/'sealed_panels.json',dict(lock_sha256=sha(OUT/'lock.json'),
                                      files={p.name:sha(p) for p in paths},expected_rows=83))
    print('SEALED all 83 panels; no outcome graph parsed during generation.',flush=True)


def verify_seal():
    seal=json.loads((OUT/'sealed_panels.json').read_text(encoding='utf-8'))
    expected={'generation.json',*(f'panels_{s}.json' for s in CONFIG['seeds'])}
    assert set(seal['files'])==expected
    assert seal['lock_sha256']==sha(OUT/'lock.json')
    assert all(sha(OUT/name)==digest for name,digest in seal['files'].items())
    setup=json.loads((OUT/'generation.json').read_text(encoding='utf-8'))
    rows=setup['rows'][:]
    assert [(r['method'],r['seed']) for r in rows]==[(m,None) for m in CONFIG['deterministic']]
    for seed in CONFIG['seeds']:
        part=json.loads((OUT/f'panels_{seed}.json').read_text(encoding='utf-8'))['rows']
        assert {(r['method'],r['seed']) for r in part}=={(m,seed) for m in CONFIG['timed']}
        assert len(part)==4
        rows+=part
    assert len(rows)==seal['expected_rows']==83
    return setup,rows


def cross(roles,e):
    u,v=roles[e[:,0]],roles[e[:,1]]
    return int(np.count_nonzero((u!=0)&(v!=0)&(u!=v)))


def novel(outcome,excluded):
    return outcome-set().union(*excluded)


def score():
    lock(); setup,rows=verify_seal(); names=cohort()
    assert names==setup['names']
    train_all=edges(INPUTS['train']); outcome=edges(INPUTS['outcome'])
    excluded=[edges(INPUTS[key]) for key in CONFIG['exclusion']]
    later=novel(outcome,excluded)
    train=restricted(train_all,names); test=restricted(later,names)
    a=matrix(len(names),train)
    scored=[]
    for row in rows:
        db,ad=row['db'],row['ad']
        assert len(db)==len(set(db))==100 and len(ad)==len(set(ad))==100
        assert not set(db)&set(ad) and all(0<=v<len(names) for v in db+ad)
        roles=np.zeros(len(names),dtype=np.int8); roles[db]=1; roles[ad]=2
        counts=a[:,db+ad].sum(axis=1)
        terms=1.-np.power(.2,counts); terms[roles!=0]=1.
        coverage=float(terms.sum()); train_cross=cross(roles,train)
        assert train_cross==row['train_cross'] and abs(coverage-row['coverage'])<1e-6
        assert coverage+1e-7>=setup['coverage_floor']
        elapsed=row['search_seconds']; allowance=setup['allowed_search_seconds']
        scored.append({k:v for k,v in row.items() if k not in ['db','ad','meta']} |
                      dict(test_cross=cross(roles,test),coverage_ratio=coverage/setup['reference_coverage'],
                           timing_valid=elapsed<=allowance+max(.02,.05*allowance) if row['method'] in CONFIG['timed'] else None))
    summary={}
    proposed=next(r['test_cross'] for r in scored if r['method']=='feasible_role_greedy')
    for method in CONFIG['deterministic']+CONFIG['timed']:
        group=[r for r in scored if r['method']==method]
        values=[r['test_cross'] for r in group]; median=float(np.median(values))
        summary[method]=dict(n=len(group),test_median=median,test_min=min(values),test_max=max(values),
             train_median=float(np.median([r['train_cross'] for r in group])),
             coverage_ratio_min=min(r['coverage_ratio'] for r in group),
             coverage_ratio_max=max(r['coverage_ratio'] for r in group),
             proposed_delta_median=proposed-median,
             proposed_relative_to_median=(proposed/median-1) if median else None,
             proposed_positive=sum(proposed>v for v in values),proposed_equal=sum(proposed==v for v in values),
             proposed_negative=sum(proposed<v for v in values),
             timing_invalid=sum(r['timing_valid'] is False for r in group))
    info=dict(candidates=len(names),train_edges=len(train),train_nodes=setup['train_nodes'],
              train_isolates=len(names)-setup['train_nodes'],outcome_full=len(outcome),
              outcome_overlaps={k:len(outcome&v) for k,v in zip(CONFIG['exclusion'],excluded)},
              novel_full=len(later),novel_in_cohort=len(test),
              novel_train_visible=sum(a[u].any() and a[v].any() for u,v in test),
              search_allowance=setup['allowed_search_seconds'],common_seconds=setup['common_seconds'])
    write(OUT/'scores.json',dict(info=info,summary=summary,rows=scored,
                               seal_sha256=sha(OUT/'sealed_panels.json')))
    with (OUT/'scores.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(scored[0])); writer.writeheader(); writer.writerows(scored)
    print(json.dumps(dict(info=info,summary=summary),ensure_ascii=False,indent=2))


if __name__=='__main__':
    {'freeze':lambda:lock(True),'generate':generate,'score':score}[sys.argv[1]]()
