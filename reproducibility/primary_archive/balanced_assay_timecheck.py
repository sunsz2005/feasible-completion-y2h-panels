"""Third development round: matched-time baselines and targeted ablations."""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from balanced_assay_pilot import ROOT, State, coverage_order, descent, inputs, matrix, partition, sha
from balanced_assay_tradeoff import feasible_role_greedy, lock as tradeoff_lock

OUT = ROOT / 'balanced_assay_timecheck'
CONFIG = dict(budget=200, guard=.98, seeds=list(range(20261020,20261040)),
              replacement_shortlist=128, perturb_steps=[2,4,8],
              timing_tolerance_fraction=.05, timing_tolerance_seconds=.02)


def role_blind(a, panel, floor):
    """Same completion rule, ranking only by degree (no cross-role gain)."""
    n,k=len(a),len(panel)
    witness=np.zeros(n,dtype=bool); witness[list(panel)]=True
    counts=a[:,list(panel)].sum(axis=1).astype(np.int32)
    committed=np.zeros(n,dtype=bool)
    roles=[[],[]]
    ranked=np.lexsort((np.arange(n),-a.sum(axis=1).astype(np.int32)))
    for step in range(k):
        fallback=next(int(u) for u in ranked if witness[u] and not committed[u])
        removable=np.flatnonzero(witness & ~committed)
        losses=np.where(a[:,removable],.8*np.power(.2,np.maximum(counts[:,None]-1,0)),0.)
        losses[witness]=0.
        losses=losses.sum(axis=0)+np.power(.2,counts[removable])
        removals=removable[np.lexsort((removable,losses))[:8]]
        chosen,removed=fallback,None
        for u in [int(x) for x in ranked if not committed[x]][:16]:
            if witness[u]:
                chosen=u
                break
            new=counts[:,None]+a[:,u,None].astype(np.int32)-a[:,removals]
            terms=1.-np.power(.2,new)
            terms[witness]=1.; terms[u]=1.
            terms[removals,np.arange(len(removals))]=1.-np.power(.2,new[removals,np.arange(len(removals))])
            scores=terms.sum(axis=0)
            if float(scores.max())+1e-8>=floor:
                removed=int(removals[int(np.argmax(scores))]); chosen=u
                break
        if removed is not None:
            witness[removed]=False; witness[chosen]=True
            counts+=a[:,chosen].astype(np.int32)-a[:,removed]
        committed[chosen]=True; roles[step%2].append(chosen)
    # Re-optimize roles exactly as in the full construction. The ablation is in selection.
    s=State(a,*roles)
    assert s.coverage()+1e-7>=floor
    descent(s,0.,replace=False)
    return s


def timed_descent(s,floor,deadline,replace):
    moves=0
    while time.perf_counter()<deadline:
        d,p=np.flatnonzero(s.roles==1),np.flatnonzero(s.roles==2)
        gains=s.cd[d,None]-s.cp[d,None]+s.cp[p][None,:]-s.cd[p][None,:]+2*s.a[np.ix_(d,p)]
        i,j=np.unravel_index(int(np.argmax(gains)),gains.shape)
        best=int(gains[i,j]),int(d[i]),int(p[j])
        if replace:
            chosen,outside=np.flatnonzero(s.roles),np.flatnonzero(s.roles==0)
            if len(outside):
                opposite=np.where(s.roles[chosen,None]==1,s.cp[outside],s.cd[outside])
                current=np.where(s.roles[chosen]==1,s.cp[chosen],s.cd[chosen])
                rg=opposite-current[:,None]
                count=min(CONFIG['replacement_shortlist'],rg.size)
                top=np.argpartition(rg.ravel(),-count)[-count:]
                top=sorted(top,key=lambda x:(-int(rg.ravel()[x]),int(x)))
                for f in top:
                    if time.perf_counter()>=deadline:
                        return moves
                    ii,jj=np.unravel_index(f,rg.shape)
                    gain=int(rg[ii,jj])
                    if gain<=max(0,best[0]):
                        break
                    u,v=int(chosen[ii]),int(outside[jj])
                    trial=s.copy(); trial.swap(u,v)
                    if trial.coverage()+1e-8>=floor:
                        best=gain,u,v
                        break
        if time.perf_counter()>=deadline or best[0]<=0:
            break
        s.swap(best[1],best[2]); moves+=1
    return moves


def multistart(base,floor,seconds,seed,replace):
    started=time.perf_counter()
    deadline=started+max(0.,seconds-.01)
    rng=np.random.default_rng(seed)
    panel=np.flatnonzero(base.roles)
    degree_order=np.argsort(-base.a.sum(axis=1).astype(np.int32),kind='stable')
    best=base.copy(); starts=0; perturb_accepted=0
    # First start is the original deterministic baseline.
    while time.perf_counter()<deadline:
        if starts==0:
            s=base.copy()
        elif not replace:
            perm=rng.permutation(panel)
            s=State(base.a,perm[:len(panel)//2],perm[len(panel)//2:])
        else:
            s=(best if starts%2 else base).copy()
            steps=CONFIG['perturb_steps'][(starts-1)%3]
            for _ in range(steps):
                if time.perf_counter()>=deadline:
                    break
                u=int(rng.choice(np.flatnonzero(s.roles)))
                if rng.random()<.5:
                    candidates=degree_order[:min(1000,len(base.a))]
                    candidates=candidates[s.roles[candidates]!=s.roles[u]]
                else:
                    candidates=np.flatnonzero(s.roles!=s.roles[u])
                if not len(candidates):
                    continue
                v=int(rng.choice(candidates))
                trial=s.copy(); trial.swap(u,v)
                if trial.coverage()+1e-8>=floor:
                    s=trial; perturb_accepted+=1
                    if s.score>best.score:
                        best=s.copy()
        timed_descent(s,floor,deadline,replace)
        starts+=1
        if s.score>best.score:
            best=s.copy()
    elapsed=time.perf_counter()-started
    return best,dict(search_seconds=elapsed,starts=starts,perturb_accepted=perturb_accepted)


def lock():
    tradeoff_lock()
    OUT.mkdir(exist_ok=True)
    paths=[Path(__file__),ROOT/'balanced_assay_timecheck_protocol.md',
           ROOT/'balanced_assay_tradeoff/lock.json',
           Path(__file__).with_name('balanced_assay_tradeoff.py'),
           Path(__file__).with_name('balanced_assay_pilot.py')]
    value=dict(config=CONFIG,files={str(p.relative_to(ROOT.parent)):sha(p) for p in paths},
               interpretation='third development round; same previously used edge splits')
    path=OUT/'lock.json'
    if path.exists():
        assert json.loads(path.read_text())==value
    else:
        path.write_text(json.dumps(value,indent=2),encoding='utf-8')


def run():
    lock()
    names,sources=inputs()
    for source,edges in sources.items():
        for seed in CONFIG['seeds']:
            path=OUT/f'{source}_{seed}.json'
            if path.exists():
                continue
            order=np.random.default_rng(seed).permutation(len(edges))
            cut=int(.8*len(edges)); train,test=edges[order[:cut]],edges[order[cut:]]
            a=matrix(len(names),train)
            started=time.perf_counter()
            panel=coverage_order(a,CONFIG['budget'])
            base=partition(a,panel)
            common=time.perf_counter()-started
            reference=base.coverage(); floor=CONFIG['guard']*reference
            started=time.perf_counter()
            proposed,_=feasible_role_greedy(a,panel,floor)
            budget=time.perf_counter()-started
            states={'coverage_partition':(base,0.,{}),'feasible_role_greedy':(proposed,budget,{})}
            # Alternate the two baseline orders by seed to reduce systematic order bias.
            modes=[False,True] if seed%2==0 else [True,False]
            for replace in modes:
                s,meta=multistart(base,floor,budget,seed+(100000 if replace else 0),replace)
                states['multistart_joint' if replace else 'multistart_roles']=(s,meta['search_seconds'],meta)
            started=time.perf_counter()
            blind=role_blind(a,panel,floor)
            states['ablation_no_role_gain']=(blind,time.perf_counter()-started,{})
            started=time.perf_counter()
            unguarded,_=feasible_role_greedy(a,panel,0.)
            states['ablation_no_coverage_guard']=(unguarded,time.perf_counter()-started,{})
            # All panels are fixed before the internal test edges are evaluated.
            rows,panels=[],{}
            for method,(s,elapsed,meta) in states.items():
                u,v=s.roles[test[:,0]],s.roles[test[:,1]]
                idx=np.flatnonzero(s.roles)
                internal=int(a[np.ix_(idx,idx)].sum())//2
                external_nodes=len(set(idx)-set(panel))
                roles_changed=int(sum(s.roles[i]!=base.roles[i] for i in set(idx)&set(panel)))
                matched=method in ('multistart_roles','multistart_joint')
                rows.append(dict(source=source,seed=seed,budget=CONFIG['budget'],method=method,
                                 train_cross=s.score,test_cross=int(np.count_nonzero((u!=0)&(v!=0)&(u!=v))),
                                 reference_coverage=reference,coverage=s.coverage(),coverage_ratio=s.coverage()/reference,
                                 feasible=s.coverage()+1e-8>=floor,common_seconds=common,search_seconds=elapsed,
                                 total_seconds=common+elapsed,allowed_search_seconds=budget,
                                 timing_valid=(elapsed<=budget+max(.02,.05*budget)) if matched else True,
                                 starts=meta.get('starts',0),perturb_accepted=meta.get('perturb_accepted',0),
                                 panel_internal_train_edges=internal,replaced_proteins=external_nodes,
                                 shared_proteins_roles_changed=roles_changed,
                                 cut_fraction=s.score/internal if internal else 0.))
                panels[method]={role:[names[i] for i in np.flatnonzero(s.roles==code)]
                                for role,code in [('DB',1),('AD',2)]}
            value=dict(rows=rows,panels=panels,train_sha256=hashlib.sha256(train.tobytes()).hexdigest(),
                       test_sha256=hashlib.sha256(test.tobytes()).hexdigest())
            temp=path.with_suffix('.pending')
            temp.write_text(json.dumps(value,indent=2),encoding='utf-8'); temp.replace(path)
            print(f'Finished {source} {seed} candidate_search={budget:.3f}s',flush=True)


if __name__=='__main__':
    run()
