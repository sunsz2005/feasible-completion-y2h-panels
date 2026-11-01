"""One fixed diagnostic: identical feasible solver, two paid initializations."""
import hashlib
import json
import platform
import sys
import time
from pathlib import Path
import numpy as np
import balanced_assay_external as external
from balanced_assay_pilot import State,inputs,matrix,coverage_order,partition,sha
from balanced_assay_tradeoff import feasible_role_greedy
from balanced_assay_timecheck import multistart

ROOT=Path(__file__).resolve().parent/'ppi_data_audit'
OUT=ROOT/'balanced_assay_initialization_time_20261024'
CONFIG=dict(sources=['YDIP','YBIOGRID','SpaceI'],seeds=list(range(20261020,20261040)),
            multipliers=[1,5,20],methods=['A','B'],budget=200,guard=.98,q=.8,
            primary_multiplier=5,relative_gate=.05,positive_gate=14,
            timing_fraction=.05,timing_seconds=.02,expected_panels=360)


def scalar(value):
    if isinstance(value,np.generic): return value.item()
    raise TypeError(str(type(value)))


def write(path,data):
    pending=path.with_suffix('.pending')
    pending.write_text(json.dumps(data,ensure_ascii=False,indent=2,default=scalar),encoding='utf-8')
    pending.replace(path)


def digest(array):
    return hashlib.sha256(array.tobytes()).hexdigest()


def lock(create=False):
    external.lock(); external.verify_seal()
    paths=[Path(__file__),OUT/'protocol.md',OUT/'closure_snapshot.md',external.OUT/'lock.json']
    value=dict(config=CONFIG,files={str(p.resolve()):sha(p) for p in paths},
               python=platform.python_version(),numpy=np.__version__)
    path=OUT/'lock.json'
    if create:
        assert not path.exists()
        write(path,value)
    else:
        assert json.loads(path.read_text(encoding='utf-8'))==value,'Frozen protocol or implementation changed'


def initialize(a,method,floor):
    selected=coverage_order(a,CONFIG['budget']); reference=partition(a,selected)
    candidate=reference
    if method=='B': candidate,_=feasible_role_greedy(a,selected,floor)
    return reference,candidate


def pipeline(a,method,floor,budget,seed):
    started=time.perf_counter()
    reference,initial=initialize(a,method,floor)
    initialization_seconds=time.perf_counter()-started
    tolerance=max(CONFIG['timing_seconds'],CONFIG['timing_fraction']*budget)
    fallback=initialization_seconds>budget+tolerance
    if fallback: initial=reference
    initialization_train=int(initial.score)
    remaining=max(0.,budget-(time.perf_counter()-started))
    meta=dict(search_seconds=0.,starts=0,perturb_accepted=0)
    result=initial
    if remaining>0 and not fallback:
        # Both arms invoke precisely the same frozen feasible joint-search function.
        result,meta=multistart(initial,floor,remaining,seed,True)
    total=time.perf_counter()-started
    return result,dict(initialization_seconds=initialization_seconds,
       initialization_train=initialization_train,initialization_fallback=fallback,
       initialization_over_target=initialization_seconds>budget,
       allowed_seconds=budget,search_allowance=remaining,total_seconds=total,
       timing_valid=total<=budget+tolerance,**meta)


def calibration(a,names,train,source,seed,key):
    path=OUT/f'instance_{key}.json'
    if path.exists():
        value=json.loads(path.read_text(encoding='utf-8'))
        assert value['names']==names and value['train_sha256']==digest(train)
        return value
    started=time.perf_counter()
    selected=coverage_order(a,CONFIG['budget']); reference=partition(a,selected)
    reference_coverage=reference.coverage(); floor=.98*reference_coverage
    proposed,_=feasible_role_greedy(a,selected,floor)
    t0=time.perf_counter()-started
    original_path=(external.OUT/'generation.json' if source=='SpaceI' else ROOT/f'balanced_assay_timecheck/{source}_{seed}.json')
    original=json.loads(original_path.read_text(encoding='utf-8'))
    if source!='SpaceI':
        assert original['train_sha256']==digest(train)
        original_panels=original['panels']
        def ids(method,role): return [names.index(x) for x in original_panels[method][role]]
        original_a=State(a,ids('coverage_partition','DB'),ids('coverage_partition','AD')).roles
        original_b=State(a,ids('feasible_role_greedy','DB'),ids('feasible_role_greedy','AD')).roles
    else:
        original_a,original_b=[State(a,r['db'],r['ad']).roles for method in ['coverage_partition','feasible_role_greedy']
                              for r in original['rows'] if r['method']==method]
    assert np.array_equal(reference.roles,original_a) and np.array_equal(proposed.roles,original_b)
    value=dict(source=source,split_seed=seed if source!='SpaceI' else None,names=names,
       train_sha256=digest(train),train_edges=len(train),reference_coverage=reference_coverage,
       floor=floor,t0_seconds=t0,calibration_train_a=int(reference.score),calibration_train_b=int(proposed.score),
       original_initializations_identical=True,original_file_sha256=sha(original_path))
    write(path,value)
    print(f'CALIBRATED {key} T0={t0:.3f}s',flush=True)
    return value


def instance(source,seed,yeast_names,sources):
    if source=='SpaceI':
        names=external.cohort(); train=external.restricted(external.edges(external.INPUTS['train']),names)
        key='SpaceI'
    else:
        names=yeast_names; all_edges=sources[source]
        order=np.random.default_rng(seed).permutation(len(all_edges)); cut=int(.8*len(order))
        train=all_edges[order[:cut]]; key=f'{source}_{seed}'
    return key,names,train


def expected_files():
    return [f'panel_{source}_{seed}_{mult}_{method}.json'
            for seed in CONFIG['seeds'] for source in CONFIG['sources']
            for mult in CONFIG['multipliers'] for method in CONFIG['methods']]


def generate():
    lock(); assert not (OUT/'sealed_panels.json').exists()
    names,sources=inputs(); done=sum((OUT/f).exists() for f in expected_files())
    for i,seed in enumerate(CONFIG['seeds']):
        for j,source in enumerate(CONFIG['sources']):
            key,nodes,train=instance(source,seed,names,sources); a=matrix(len(nodes),train)
            cal=calibration(a,nodes,train,source,seed,key)
            multipliers=CONFIG['multipliers'][i%3:]+CONFIG['multipliers'][:i%3]
            for mult in multipliers:
                methods=['A','B'] if (i+j)%2==0 else ['B','A']
                for method in methods:
                    path=OUT/f'panel_{source}_{seed}_{mult}_{method}.json'
                    if path.exists(): continue
                    progress=dict(done=done,total=360,source=source,seed=seed,multiplier=mult,method=method,
                                  target_seconds=cal['t0_seconds']*mult,phase='generating')
                    write(OUT/'progress.json',progress)
                    print(f'START {done+1}/360 {source} {seed} {mult}x {method}',flush=True)
                    state,meta=pipeline(a,method,cal['floor'],cal['t0_seconds']*mult,seed)
                    row=dict(source=source,seed=seed,multiplier=mult,method=method,
                         db=np.flatnonzero(state.roles==1).tolist(),ad=np.flatnonzero(state.roles==2).tolist(),
                         train_cross=int(state.score),coverage=state.coverage(),instance=key,
                         instance_sha256=sha(OUT/f'instance_{key}.json'),**meta)
                    assert len(row['db'])==len(row['ad'])==100 and row['coverage']+1e-7>=cal['floor']
                    write(path,row); done+=1
                    print(f'DONE {done}/360 timing_valid={meta["timing_valid"]}',flush=True)
    files=expected_files()+[f'instance_{s}_{seed}.json' for s in ['YDIP','YBIOGRID'] for seed in CONFIG['seeds']]+['instance_SpaceI.json']
    assert done==360
    write(OUT/'sealed_panels.json',dict(lock_sha256=sha(OUT/'lock.json'),files={name:sha(OUT/name) for name in files},panels=360))
    write(OUT/'progress.json',dict(done=360,total=360,phase='sealed_not_scored'))
    print('SEALED 360 panels; evaluation not yet scored.',flush=True)


def verify_seal():
    seal=json.loads((OUT/'sealed_panels.json').read_text(encoding='utf-8'))
    expected=set(expected_files())|{f'instance_{s}_{seed}.json' for s in ['YDIP','YBIOGRID'] for seed in CONFIG['seeds']}|{'instance_SpaceI.json'}
    assert set(seal['files'])==expected and seal['panels']==360
    assert seal['lock_sha256']==sha(OUT/'lock.json')
    assert all(sha(OUT/name)==value for name,value in seal['files'].items())
    return [json.loads((OUT/name).read_text(encoding='utf-8')) for name in expected_files()]


def score():
    lock(); rows=verify_seal(); names,sources=inputs(); scored=[]
    # Outcome edges are only parsed/used here, after every output has been sealed.
    outcome=external.novel(external.edges(external.INPUTS['outcome']),
                         [external.edges(external.INPUTS[key]) for key in external.CONFIG['exclusion']])
    for seed in CONFIG['seeds']:
        for source in CONFIG['sources']:
            key,nodes,train=instance(source,seed,names,sources)
            cal=json.loads((OUT/f'instance_{key}.json').read_text(encoding='utf-8'))
            assert digest(train)==cal['train_sha256']
            if source=='SpaceI': test=external.restricted(outcome,nodes)
            else:
                e=sources[source]; order=np.random.default_rng(seed).permutation(len(e)); test=e[order[int(.8*len(e)):]]
                old=json.loads((ROOT/f'balanced_assay_timecheck/{source}_{seed}.json').read_text())
                assert digest(test)==old['test_sha256']
            a=matrix(len(nodes),train)
            for row in [r for r in rows if r['seed']==seed and r['source']==source]:
                db,ad=row['db'],row['ad']
                assert len(db)==len(set(db))==len(ad)==len(set(ad))==100
                assert not set(db)&set(ad) and all(0<=u<len(nodes) for u in db+ad)
                roles=np.zeros(len(nodes),dtype=np.int8); roles[db]=1; roles[ad]=2
                terms=1.-np.power(.2,a[:,db+ad].sum(axis=1)); terms[roles!=0]=1.
                coverage=float(terms.sum()); train_score=external.cross(roles,train)
                assert abs(coverage-row['coverage'])<1e-6 and train_score==row['train_cross']
                assert coverage+1e-7>=cal['floor'] and row['instance_sha256']==sha(OUT/f'instance_{key}.json')
                assert abs(row['allowed_seconds']-row['multiplier']*cal['t0_seconds'])<1e-8
                valid=row['total_seconds']<=row['allowed_seconds']+max(.02,.05*row['allowed_seconds'])
                assert bool(row['timing_valid'])==bool(valid)
                scored.append({k:v for k,v in row.items() if k not in ['db','ad']}|
                             dict(test_cross=external.cross(roles,test),coverage_ratio=coverage/cal['reference_coverage']))
    paired=[]; summaries=[]
    for source in CONFIG['sources']:
        for mult in CONFIG['multipliers']:
            pair=[]
            for seed in CONFIG['seeds']:
                x={r['method']:r for r in scored if r['source']==source and r['seed']==seed and r['multiplier']==mult}
                va,vb=x['A'],x['B']; delta=vb['test_cross']-va['test_cross']
                row=dict(source=source,seed=seed,multiplier=mult,a=va['test_cross'],b=vb['test_cross'],delta=delta,
                         relative=delta/va['test_cross'] if va['test_cross'] else None,
                         train_delta=vb['train_cross']-va['train_cross'],timing_valid=va['timing_valid'] and vb['timing_valid'])
                pair.append(row); paired.append(row)
            relative=[r['relative'] for r in pair if r['relative'] is not None]
            summary=dict(source=source,multiplier=mult,a_median=float(np.median([r['a'] for r in pair])),
                         b_median=float(np.median([r['b'] for r in pair])),
                         delta_median=float(np.median([r['delta'] for r in pair])),
                         relative_median=float(np.median(relative)) if relative else None,
                         positive=sum(r['delta']>0 for r in pair),equal=sum(r['delta']==0 for r in pair),
                         negative=sum(r['delta']<0 for r in pair),valid_pairs=sum(r['timing_valid'] for r in pair),
                         train_delta_median=float(np.median([r['train_delta'] for r in pair])))
            summaries.append(summary)
    primary={s['source']:s for s in summaries if s['multiplier']==5}
    valid=all(s['valid_pairs']==20 for s in primary.values())
    gain=all(primary[s]['relative_median'] is not None and primary[s]['relative_median']>=.05 and primary[s]['positive']>=14 for s in ['YDIP','YBIOGRID'])
    human=primary['SpaceI']['delta_median']>0
    decision='CONTINUE_METHOD_DRAFT' if valid and gain and human else 'CLOSE_AS_LIMITED_ENGINEERING' if valid else 'INVALID_PRIMARY_TIMING'
    write(OUT/'scores.json',dict(summary=summaries,paired=paired,rows=scored,decision=decision,
                 primary_checks=dict(timing_valid=valid,yeast_gate=gain,human_positive=human),
                 seal_sha256=sha(OUT/'sealed_panels.json')))
    write(OUT/'progress.json',dict(done=360,total=360,phase='scored',decision=decision))
    print(json.dumps(dict(summary=summaries,decision=decision),ensure_ascii=False,indent=2))


if __name__=='__main__':
    {'freeze':lambda:lock(True),'generate':generate,'score':score}[sys.argv[1]]()
