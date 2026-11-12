"""Second-round development: feasible construction and coverage/cut tradeoffs."""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from balanced_assay_pilot import (ROOT, State, coverage_order, descent, inputs,
                                 lock as original_lock, matrix, partition, sha)

OUT = ROOT / 'balanced_assay_tradeoff'
CONFIG = dict(seeds=list(range(20261020, 20261040)), budgets=[50, 100, 200],
              guards=[.90, .95, .98, .99, 1.], candidate_shortlist=16,
              removal_shortlist=8, train_fraction=.8)


def feasible_role_greedy(a, panel, floor):
    """Role-gain greedy with a feasible completion witness, using training only.

    The witness always contains every committed protein and k proteins total.
    Outside candidates replace an uncommitted witness protein only if the full
    witness stays coverage-feasible. Hence the finished balanced panel is feasible.
    """
    n, k = len(a), len(panel)
    witness = np.zeros(n, dtype=bool)
    witness[list(panel)] = True
    counts = a[:, list(panel)].sum(axis=1).astype(np.int32)
    committed = np.zeros(n, dtype=bool)
    roles = [[], []]
    degree = a.sum(axis=1)
    external = 0
    for step in range(k):
        side = step % 2
        opposite = a[:, roles[1-side]].sum(axis=1)
        ranked = np.lexsort((np.arange(n), -degree.astype(np.int32), -opposite.astype(np.int32)))
        # The best available witness member is always a feasible fallback.
        fallback = next(int(u) for u in ranked if witness[u] and not committed[u])
        # Rank uncommitted removals by exact coverage loss, then ID.
        removable = np.flatnonzero(witness & ~committed)
        base = 1. - np.power(.2, counts)
        losses = np.where(a[:, removable], .8 * np.power(.2, np.maximum(counts[:, None]-1, 0)), 0.)
        losses[witness] = 0.
        losses = losses.sum(axis=0) + np.power(.2, counts[removable])
        removals = removable[np.lexsort((removable, losses))[:CONFIG['removal_shortlist']]]
        chosen, removed = fallback, None
        for u in [int(x) for x in ranked if not committed[x]][:CONFIG['candidate_shortlist']]:
            if witness[u]:
                chosen = u
                break
            new_counts = counts[:, None] + a[:, u, None].astype(np.int32) - a[:, removals]
            coverage_terms = 1. - np.power(.2, new_counts)
            coverage_terms[witness] = 1.
            coverage_terms[u] = 1.
            coverage_terms[removals, np.arange(len(removals))] = 1. - np.power(.2, new_counts[removals, np.arange(len(removals))])
            scores = coverage_terms.sum(axis=0)
            if float(scores.max()) + 1e-8 >= floor:
                # Preserve as much coverage as possible at this role gain.
                removed = int(removals[int(np.argmax(scores))])
                chosen = u
                break
        if removed is not None:
            witness[removed] = False
            witness[chosen] = True
            counts += a[:, chosen].astype(np.int32) - a[:, removed]
            external += 1
        committed[chosen] = True
        roles[side].append(chosen)
    s = State(a, *roles)
    assert s.coverage() + 1e-7 >= floor
    assert np.array_equal(witness, s.roles != 0)
    descent(s, 0., replace=False)
    return s, external


def lock():
    original_lock()
    OUT.mkdir(exist_ok=True)
    paths = [Path(__file__), Path(__file__).with_name('balanced_assay_pilot.py'),
             ROOT / 'balanced_assay_tradeoff_protocol.md',
             ROOT / 'balanced_assay_pilot/lock.json',
             ROOT / 'YeRI_supplementary_data_15.xlsx',
             ROOT / 'YDIP/interactions.txt', ROOT / 'YBIOGRID/interactions.txt']
    payload = dict(config=CONFIG, files={str(p.relative_to(ROOT.parent)): sha(p) for p in paths},
                   interpretation='development on previously inspected edge splits; not confirmation')
    path = OUT / 'lock.json'
    if path.exists():
        assert json.loads(path.read_text()) == payload
    else:
        path.write_text(json.dumps(payload, indent=2), encoding='utf-8')


def run():
    lock()
    names, source_edges = inputs()
    for source, edges in source_edges.items():
        for seed in CONFIG['seeds']:
            path = OUT / f'{source}_{seed}.json'
            if path.exists():
                continue
            order = np.random.default_rng(seed).permutation(len(edges))
            cut = int(.8 * len(edges))
            train, test = edges[order[:cut]], edges[order[cut:]]
            a = matrix(len(names), train)
            t = time.perf_counter()
            cov_order = coverage_order(a, 200)
            cov_time = time.perf_counter() - t
            records, panels, selected = [], {}, {}
            for k in CONFIG['budgets']:
                t = time.perf_counter()
                base = partition(a, cov_order[:k])
                base_time = time.perf_counter() - t + cov_time
                ref = base.coverage()
                states = {}
                for guard in CONFIG['guards']:
                    floor = guard * ref
                    states[guard, 'coverage_partition'] = base, base_time, 0
                    t = time.perf_counter()
                    greedy, changed = feasible_role_greedy(a, cov_order[:k], floor)
                    greedy_time = time.perf_counter() - t + cov_time
                    states[guard, 'feasible_role_greedy'] = greedy, greedy_time, changed
                    t = time.perf_counter()
                    joint = base.copy()
                    moves = descent(joint, floor)
                    states[guard, 'joint_from_coverage'] = joint, base_time + time.perf_counter() - t, moves
                    t = time.perf_counter()
                    strong = greedy.copy()
                    moves = descent(strong, floor)
                    states[guard, 'feasible_role_greedy_swap'] = strong, greedy_time + time.perf_counter() - t, moves
                    # Deployable portfolio: selection ONLY by training cut, not test scores.
                    winner = min(('joint_from_coverage', 'feasible_role_greedy_swap'),
                                 key=lambda m: (-states[guard, m][0].score, m))
                    chosen = states[guard, winner][0]
                    elapsed = states[guard, 'joint_from_coverage'][1] + states[guard, 'feasible_role_greedy_swap'][1]
                    states[guard, 'train_selected_portfolio'] = chosen, elapsed, 0
                    selected[f'{k}/{guard}'] = winner
                # Panels are all frozen before evaluating ANY test edge for this budget.
                for (guard, method), (s, elapsed, moves) in states.items():
                    key = f'{k}/{guard}/{method}'
                    u, v = s.roles[test[:, 0]], s.roles[test[:, 1]]
                    hit = int(np.count_nonzero((u != 0) & (v != 0) & (u != v)))
                    records.append(dict(source=source, seed=seed, budget=k, guard=guard, method=method,
                                        coverage=s.coverage(), reference_coverage=ref,
                                        coverage_ratio=s.coverage()/ref, train_cross=s.score,
                                        test_cross=hit, test_edges=len(test), test_fraction=hit/len(test),
                                        feasible=s.coverage()+1e-8 >= guard*ref, seconds=elapsed, moves=moves))
                    panels[key] = {role: [names[i] for i in np.flatnonzero(s.roles == code)]
                                   for role, code in [('DB', 1), ('AD', 2)]}
            path.write_text(json.dumps(dict(rows=records, panels=panels, train_selected=selected,
                                           train_sha256=hashlib.sha256(train.tobytes()).hexdigest(),
                                           test_sha256=hashlib.sha256(test.tobytes()).hexdigest()), indent=2), encoding='utf-8')
            print(f'Finished {source} {seed}', flush=True)


if __name__ == '__main__':
    run()
