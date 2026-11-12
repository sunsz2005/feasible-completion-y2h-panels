"""Controlled selection-policy x coverage-enforcement-timing development study."""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from balanced_assay_pilot import ROOT,coverage_order,inputs,matrix,partition,sha
from balanced_assay_tradeoff import feasible_role_greedy
from balanced_assay_timecheck import role_blind,timed_descent
from balanced_assay_ep_adapter import repair,lock as previous_lock

OUT=ROOT/'balanced_assay_gate_control'
CONFIG=dict(seeds=list(range(20261020,20261040)),budget=200,guard=.98,
            policies=['role','degree'],repair='reference-directed 8x8 coverage repair',
            timing_tolerance_fraction=.05,timing_tolerance_seconds=.02)


def construct(a,panel,floor,policy):
    return feasible_role_greedy(a,panel,floor)[0] if policy=='role' else role_blind(a,panel,floor)


def post_pipeline(a,panel,base,floor,budget,policy):
    start=time.perf_counter(); deadline=start+max(0.,budget-.01)
    raw=construct(a,panel,0.,policy)
    construction_seconds=time.perf_counter()-start
    raw_score,raw_coverage=raw.score,raw.coverage()
    repaired,meta=repair(raw.copy(),base,floor,deadline)
    repair_seconds=time.perf_counter()-start-construction_seconds
    repaired_actual=repaired.copy()
    best_repaired=repaired.copy() if repaired.score>base.score else base.copy()
    search_start=time.perf_counter()
    moves=timed_descent(repaired,floor,deadline,True)
    best=repaired if repaired.score>best_repaired.score else best_repaired.copy()
    elapsed=time.perf_counter()-start
    return raw,best_repaired,best,dict(policy=policy,construction_seconds=construction_seconds,
                                    repair_seconds=repair_seconds,postsearch_seconds=time.perf_counter()-search_start,
                                    total_seconds=elapsed,postsearch_moves=moves,
                                    raw_score=raw_score,raw_coverage=raw_coverage,
                                    repaired_actual_score=repaired_actual.score,
                                    repaired_actual_coverage=repaired_actual.coverage(),
                                    raw_proteins_retained=len(set(np.flatnonzero(raw.roles))&set(np.flatnonzero(repaired_actual.roles))),
                                    repaired_selected_reference=best_repaired.score==base.score and np.array_equal(best_repaired.roles,base.roles),
                                    final_selected_reference=np.array_equal(best.roles,base.roles),**meta)


def lock():
    previous_lock(); OUT.mkdir(exist_ok=True)
    files=[Path(__file__),ROOT/'balanced_assay_gate_control_protocol.md',ROOT/'balanced_assay_ep_adapter/lock.json']
    payload=dict(config=CONFIG,files={str(p.relative_to(ROOT.parent)):sha(p) for p in files},
                 interpretation='Same previously inspected splits; controlled development pipeline comparison')
    path=OUT/'lock.json'
    if path.exists():
        assert json.loads(path.read_text())==payload
    else:
        path.write_text(json.dumps(payload,indent=2),encoding='utf-8')


def run():
    lock(); names,sources=inputs()
    for source,edges in sources.items():
        for seed in CONFIG['seeds']:
            path=OUT/f'{source}_{seed}.json'
            if path.exists():
                continue
            order=np.random.default_rng(seed).permutation(len(edges)); cut=int(.8*len(edges))
            train,test=edges[order[:cut]],edges[order[cut:]]; a=matrix(len(names),train)
            start=time.perf_counter(); panel=coverage_order(a,200); base=partition(a,panel)
            common=time.perf_counter()-start; reference=base.coverage(); floor=.98*reference
            start=time.perf_counter(); during_role=construct(a,panel,floor,'role')
            budget=time.perf_counter()-start
            start=time.perf_counter(); during_degree=construct(a,panel,floor,'degree')
            degree_time=time.perf_counter()-start
            states={'coverage_partition':(base,0.),'during_role':(during_role,budget),'during_degree':(during_degree,degree_time)}
            diagnostics=[]
            policies=['role','degree'] if seed%2==0 else ['degree','role']
            for policy in policies:
                raw,after,search,meta=post_pipeline(a,panel,base,floor,budget,policy)
                states[f'after_{policy}_raw']=(raw,meta['construction_seconds'])
                states[f'after_{policy}_repair']=(after,meta['construction_seconds']+meta['repair_seconds'])
                states[f'after_{policy}_search']=(search,meta['total_seconds'])
                diagnostics.append(meta)
            # Evaluate only after all policies and enforcement variants are fixed.
            rows,panels=[],{}
            for method,(s,elapsed) in states.items():
                u,v=s.roles[test[:,0]],s.roles[test[:,1]]
                rows.append(dict(source=source,seed=seed,method=method,budget=200,guard=.98,
                                 train_cross=s.score,test_cross=int(np.count_nonzero((u!=0)&(v!=0)&(u!=v))),
                                 coverage=s.coverage(),reference_coverage=reference,coverage_ratio=s.coverage()/reference,
                                 feasible=s.coverage()+1e-8>=floor,search_seconds=elapsed,common_seconds=common,
                                 allowed_search_seconds=budget,timing_valid=elapsed<=budget+max(.02,.05*budget)))
                panels[method]={role:[names[i] for i in np.flatnonzero(s.roles==code)] for role,code in [('DB',1),('AD',2)]}
            result=dict(rows=rows,panels=panels,diagnostics=diagnostics,
                        train_sha256=hashlib.sha256(train.tobytes()).hexdigest(),test_sha256=hashlib.sha256(test.tobytes()).hexdigest())
            pending=path.with_suffix('.pending'); pending.write_text(json.dumps(result,indent=2),encoding='utf-8'); pending.replace(path)
            print(f'Finished {source} {seed}',flush=True)


if __name__=='__main__':
    run()
