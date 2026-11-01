"""Independent EP-Prox-inspired adapter; not the authors' original implementation."""
import itertools
import json
import time
from pathlib import Path

import numpy as np
from balanced_assay_pilot import ROOT,State,coverage_order,inputs,matrix,partition,sha
from balanced_assay_tradeoff import feasible_role_greedy
from balanced_assay_timecheck import timed_descent,lock as previous_lock

OUT=ROOT/'balanced_assay_ep_adapter'
CONFIG=dict(seeds=list(range(20261020,20261040)),budget=200,guard=.98,
            iterations=100,lambda_initial=1e-10,lambda_multiplier=10,
            beta_fraction=.9,conflict_multiplier=.1,spectral_iterations=30,jitter=.001,
            repair_shortlist=8,initializations=['uniform','jitter'],timing_tolerance=.05)


def multiply(n,edges,v):
    if not len(edges):
        return np.zeros(n)
    u,w=edges[:,0],edges[:,1]
    return np.bincount(u,weights=v[w],minlength=n)+np.bincount(w,weights=v[u],minlength=n)


def prox(z,k,mu):
    # Proposition 2 applied independently to the two role segments.
    order=np.lexsort((np.arange(len(z)),-z))
    result=np.clip(z-mu,0.,1.)
    top=order[:k]; result[top]=np.clip(z[top]+mu,0.,1.)
    return result


def smooth(x,y,edges,beta,rho):
    ay=multiply(len(x),edges,y)
    return float(-2*x@ay+beta*(x@x+y@y)+rho*(x@y))


def gradient(x,y,edges,beta,rho):
    return (-2*multiply(len(x),edges,y)+2*beta*x+rho*y,
            -2*multiply(len(x),edges,x)+2*beta*y+rho*x)


def ep_scores(n,edges,k,degree_bound,seed,jitter=False,deadline=float('inf')):
    rng=np.random.default_rng(seed)
    x=np.full(n,1./(2*k)); y=x.copy()
    if jitter:
        x=np.clip(x+CONFIG['jitter']*rng.uniform(-1,1,n),0,1)
        y=np.clip(y+CONFIG['jitter']*rng.uniform(-1,1,n),0,1)
    power=np.full(n,1./np.sqrt(n))
    spectral=0.
    for _ in range(CONFIG['spectral_iterations']):
        value=multiply(n,edges,power); norm=float(np.linalg.norm(value))
        if norm<=1e-15:
            break
        power=value/norm; spectral=norm
    beta=CONFIG['beta_fraction']*spectral
    rho=CONFIG['conflict_multiplier']*spectral
    lipschitz=max(1.,2*(degree_bound+beta)+rho)
    step=1./lipschitz
    oldx,oldy=x.copy(),y.copy(); momentum=1.; penalty=CONFIG['lambda_initial']; since=0
    done=0
    for _ in range(CONFIG['iterations']):
        if time.perf_counter()>=deadline:
            break
        next_momentum=(1.+np.sqrt(1.+4.*momentum*momentum))/2.
        factor=(momentum-1.)/next_momentum
        zx=x+factor*(x-oldx); zy=y+factor*(y-oldy)
        gx,gy=gradient(zx,zy,edges,beta,rho)
        nx=prox(zx-step*gx,k,step*penalty); ny=prox(zy-step*gy,k,step*penalty)
        difference=float(np.sum((nx-x)**2)+np.sum((ny-y)**2))
        scale=max(float(np.sqrt(nx@nx+ny@ny)),1e-15)
        oldx,oldy=x,y; x,y=nx,ny; momentum=next_momentum
        done+=1; since+=1
        if np.sqrt(difference)/scale<.5 or since>=10:
            penalty*=CONFIG['lambda_multiplier']; since=0
        if difference<=1e-15:
            break
    return x,y,dict(iterations=done,penalty=penalty,spectral_estimate=spectral)


def round_roles(a,x,y,k):
    # Explicit adapter: descending score, unique proteins and fixed role capacities.
    n=len(a)
    values=np.concatenate((x,y)); ids=np.tile(np.arange(n),2); role=np.repeat([1,2],n)
    order=np.lexsort((role,ids,-values))
    chosen=set(); groups=[[],[]]
    for j in order:
        u=int(ids[j]); side=int(role[j])-1
        if u not in chosen and len(groups[side])<k:
            chosen.add(u); groups[side].append(u)
        if len(chosen)==2*k:
            break
    assert all(len(g)==k for g in groups)
    return State(a,*groups)


def repair(s,reference,floor,deadline=float('inf')):
    reference_set=set(np.flatnonzero(reference.roles))
    steps=0
    while s.coverage()+1e-8<floor:
        if time.perf_counter()>=deadline:
            return reference.copy(),dict(repair_steps=steps,fallback=True)
        selected=set(np.flatnonzero(s.roles))
        remove=np.array(sorted(selected-reference_set),dtype=np.int32)
        add=np.array(sorted(reference_set-selected),dtype=np.int32)
        assert len(remove)==len(add)>0
        counts=s.cd+s.cp; unselected=s.roles==0
        # Exact single removal loss and addition gain rank the fixed shortlists.
        losses=np.where(s.a[:,remove],.8*np.power(.2,np.maximum(counts[:,None]-1,0)),0.)
        losses[~unselected]=0.
        losses=losses.sum(axis=0)+np.power(.2,counts[remove])
        gains=.8*(np.power(.2,counts)*unselected)@s.a[:,add]
        gains+=np.power(.2,counts[add])
        remove=remove[np.lexsort((remove,losses))[:CONFIG['repair_shortlist']]]
        add=add[np.lexsort((add,-gains))[:CONFIG['repair_shortlist']]]
        best=None
        for u in remove:
            for v in add:
                if time.perf_counter()>=deadline:
                    return reference.copy(),dict(repair_steps=steps,fallback=True)
                t=s.copy(); t.swap(int(u),int(v))
                key=(t.coverage(),t.score,-int(u),-int(v))
                if best is None or key>best[0]:
                    best=key,t
        s=best[1]; steps+=1
    return s,dict(repair_steps=steps,fallback=False)


def adapter(a,edges,reference,floor,seconds,seed):
    started=time.perf_counter(); deadline=started+max(0.,seconds-.01)
    best=reference.copy(); best_before=reference.copy(); logs=[]; raw=None
    degree_bound=float(a.sum(axis=1).max())
    for i,initial in enumerate(CONFIG['initializations']):
        if time.perf_counter()>=deadline:
            break
        x,y,meta=ep_scores(len(a),edges,np.count_nonzero(reference.roles==1),degree_bound,
                            seed+i*100000,jitter=initial=='jitter',deadline=deadline)
        if time.perf_counter()>=deadline:
            logs.append(dict(initialization=initial,stage='deadline_after_ep',**meta)); break
        independent_top_x=set(np.argsort(-x,kind='stable')[:np.count_nonzero(reference.roles==1)])
        independent_top_y=set(np.argsort(-y,kind='stable')[:np.count_nonzero(reference.roles==2)])
        s=round_roles(a,x,y,np.count_nonzero(reference.roles==1))
        if raw is None or s.score>raw.score:
            raw=s.copy()
        raw_coverage=s.coverage(); raw_cut=s.score
        s,repair_meta=repair(s,reference,floor,deadline)
        if s.score>best_before.score:
            best_before=s.copy()
        if s.score>best.score:
            best=s.copy()
        moves=timed_descent(s,floor,deadline,True)
        if s.score>best.score:
            best=s.copy()
        logs.append(dict(initialization=initial,stage='complete',raw_coverage=raw_coverage,raw_cut=raw_cut,
                         independent_rounding_role_conflicts=len(independent_top_x&independent_top_y),
                         postsearch_moves=moves,**repair_meta,**meta))
    return best,best_before,raw,dict(search_seconds=time.perf_counter()-started,starts=logs)


def lock():
    previous_lock(); OUT.mkdir(exist_ok=True)
    files=[Path(__file__),ROOT/'balanced_assay_ep_adapter_protocol.md',ROOT/'balanced_assay_timecheck/lock.json']
    payload=dict(config=CONFIG,files={str(p.relative_to(ROOT.parent)):sha(p) for p in files},
                 interpretation='EP-Prox-inspired independent adapter; not an original-author replication or new external confirmation')
    target=OUT/'lock.json'
    if target.exists():
        assert json.loads(target.read_text())==payload
    else:
        target.write_text(json.dumps(payload,indent=2),encoding='utf-8')


def run():
    lock(); names,sources=inputs()
    for source,edges in sources.items():
        for seed in CONFIG['seeds']:
            path=OUT/f'{source}_{seed}.json'
            if path.exists():
                continue
            order=np.random.default_rng(seed).permutation(len(edges)); cut=int(.8*len(edges))
            train,test=edges[order[:cut]],edges[order[cut:]]; a=matrix(len(names),train)
            t=time.perf_counter(); panel=coverage_order(a,CONFIG['budget']); base=partition(a,panel)
            common=time.perf_counter()-t; ref=base.coverage(); floor=.98*ref
            t=time.perf_counter(); proposed,_=feasible_role_greedy(a,panel,floor)
            budget=time.perf_counter()-t
            adapted,repaired,raw,meta=adapter(a,train,base,floor,budget,seed)
            states={'coverage_partition':base,'feasible_role_greedy':proposed,
                    'ep_adapter_repaired':repaired,'ep_adapter_repaired_swap':adapted}
            if raw is not None:
                states['ep_adapter_unguarded_reference']=raw
            rows,panels=[],{}
            for method,s in states.items():
                u,v=s.roles[test[:,0]],s.roles[test[:,1]]
                rows.append(dict(source=source,seed=seed,method=method,budget=200,
                                 train_cross=s.score,test_cross=int(np.count_nonzero((u!=0)&(v!=0)&(u!=v))),
                                 coverage=s.coverage(),reference_coverage=ref,coverage_ratio=s.coverage()/ref,
                                 feasible=s.coverage()+1e-8>=floor,
                                 common_seconds=common,proposed_search_seconds=budget,
                                 adapter_search_seconds=meta['search_seconds'],
                                 timing_valid=meta['search_seconds']<=budget+max(.02,.05*budget)))
                panels[method]={role:[names[j] for j in np.flatnonzero(s.roles==code)]
                                for role,code in [('DB',1),('AD',2)]}
            result=dict(rows=rows,panels=panels,adapter_log=meta,
                        train_sha256=__import__('hashlib').sha256(train.tobytes()).hexdigest(),
                        test_sha256=__import__('hashlib').sha256(test.tobytes()).hexdigest())
            pending=path.with_suffix('.pending'); pending.write_text(json.dumps(result,indent=2),encoding='utf-8'); pending.replace(path)
            print(f'Finished {source} {seed}',flush=True)


if __name__=='__main__':
    run()
