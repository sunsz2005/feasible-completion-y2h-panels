"""Frozen internal development pilot; no external interaction outcomes are loaded."""
import argparse
import csv
import hashlib
import heapq
import json
import platform
import time
from pathlib import Path

import numpy as np
import openpyxl

ROOT = Path(__file__).resolve().parent / 'ppi_data_audit'
OUT = ROOT / 'balanced_assay_pilot'
CONFIG = dict(seeds=list(range(20261020, 20261040)), budgets=[50, 100, 200],
              guards=[0.98, 0.95, 1.0], q=0.8, train_fraction=0.8,
              single_rounds=80, replacement_shortlist=128,
              escape_rounds=3, proposals_per_length=64, lengths=[2, 3, 4])


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_edges(name):
    nodes, edges = set(), set()
    for line in (ROOT / name / 'interactions.txt').read_text().splitlines():
        u, v = line.split()
        nodes.update((u, v))
        if u != v:
            edges.add(tuple(sorted((u, v))))
    return nodes, edges


def inputs():
    raw = {name: read_edges(name) for name in ('YDIP', 'YBIOGRID')}
    common = raw['YDIP'][0] & raw['YBIOGRID'][0]
    book = openpyxl.load_workbook(ROOT / 'YeRI_supplementary_data_15.xlsx', read_only=True, data_only=True)
    rows = book.active.values
    assert next(rows) == ('orf_name', 'GenBank_accession', 'ORFeome_collection', 'screened_as_DB', 'screened_as_AD')
    eligible, seen = set(), set()
    for node, _, _, db, ad in rows:
        assert node not in seen
        seen.add(node)
        if db is True and ad is True and node in common:
            eligible.add(node)
    book.close()
    names = sorted(eligible)
    assert len(common) == 4329 and len(names) == 3408
    idx = {node: i for i, node in enumerate(names)}
    edges = {name: np.array(sorted((idx[u], idx[v]) for u, v in value[1]
                                  if u in idx and v in idx), dtype=np.int32)
             for name, value in raw.items()}
    return names, edges


def matrix(n, edges):
    a = np.zeros((n, n), dtype=bool)
    if len(edges):
        a[edges[:, 0], edges[:, 1]] = True
        a[edges[:, 1], edges[:, 0]] = True
    return a


def coverage_order(a, k):
    neighbors = [np.flatnonzero(row) for row in a]
    left = np.ones(len(a))
    heap = [(-float(1 + .8 * len(nb)), i) for i, nb in enumerate(neighbors)]
    heapq.heapify(heap)
    selected, order = set(), []
    while len(order) < k:
        old, u = heapq.heappop(heap)
        if u in selected:
            continue
        gain = left[u] + .8 * left[neighbors[u]].sum()
        if abs(gain + old) > 1e-9:
            heapq.heappush(heap, (-float(gain), u))
            continue
        selected.add(u)
        order.append(u)
        left[u] = 0
        left[neighbors[u]] *= .2
    return order


class State:
    def __init__(self, a, d, p):
        self.a = a
        self.roles = np.zeros(len(a), dtype=np.int8)
        self.roles[list(d)] = 1
        self.roles[list(p)] = 2
        self.cd = a[:, list(d)].sum(axis=1).astype(np.int32)
        self.cp = a[:, list(p)].sum(axis=1).astype(np.int32)
        self.score = int(self.cp[list(d)].sum())

    def copy(self):
        s = object.__new__(State)
        s.a, s.roles, s.cd, s.cp, s.score = self.a, self.roles.copy(), self.cd.copy(), self.cp.copy(), self.score
        return s

    def coverage(self):
        return float(np.where(self.roles != 0, 1., 1. - np.power(.2, self.cd + self.cp)).sum())

    def delta(self, u, v):
        ru, rv = self.roles[u], self.roles[v]
        assert ru and ru != rv
        own, other = (self.cd, self.cp) if ru == 1 else (self.cp, self.cd)
        if rv == 0:
            return int(other[v] - other[u])
        return int(own[u] - other[u] + other[v] - own[v] + 2 * int(self.a[u, v]))

    def swap(self, u, v):
        gain = self.delta(u, v)
        ru, rv = self.roles[u], self.roles[v]
        change = self.a[:, v].astype(np.int32) - self.a[:, u]
        if ru == 1:
            self.cd += change
        else:
            self.cp += change
        if rv:
            if rv == 1:
                self.cd -= change
            else:
                self.cp -= change
        self.roles[u], self.roles[v] = rv, ru
        self.score += gain


def partition(a, panel):
    degree = a.sum(axis=1)
    roles = [[], []]
    remaining = set(panel)
    while remaining:
        side = 0 if len(roles[0]) <= len(roles[1]) else 1
        opposite = roles[1 - side]
        u = min(remaining, key=lambda x: (-int(a[x, opposite].sum()), -int(degree[x]), x))
        roles[side].append(u)
        remaining.remove(u)
    s = State(a, *roles)
    # All fixed-panel baselines receive identical, exhaustive role-swap search.
    descent(s, 0., replace=False)
    return s


def role_greedy(a, k):
    degree = a.sum(axis=1)
    roles = [[], []]
    selected = set()
    for _ in range(k):
        side = 0 if len(roles[0]) <= len(roles[1]) else 1
        counts = a[:, roles[1 - side]].sum(axis=1)
        u = min((i for i in range(len(a)) if i not in selected),
                key=lambda x: (-int(counts[x]), -int(degree[x]), x))
        roles[side].append(u)
        selected.add(u)
    s = State(a, *roles)
    descent(s, 0., replace=False)
    return s


def descent(s, floor, replace=True):
    moves = 0
    for _ in range(CONFIG['single_rounds']):
        d, p = np.flatnonzero(s.roles == 1), np.flatnonzero(s.roles == 2)
        gains = s.cd[d, None] - s.cp[d, None] + s.cp[p][None, :] - s.cd[p][None, :] + 2 * s.a[np.ix_(d, p)]
        flat = int(np.argmax(gains))
        i, j = np.unravel_index(flat, gains.shape)
        best = (int(gains[i, j]), int(d[i]), int(p[j]))
        if replace:
            chosen, outside = np.flatnonzero(s.roles), np.flatnonzero(s.roles == 0)
            if len(outside):
                opposite = np.where(s.roles[chosen, None] == 1, s.cp[outside], s.cd[outside])
                current = np.where(s.roles[chosen] == 1, s.cp[chosen], s.cd[chosen])
                rg = opposite - current[:, None]
                count = min(CONFIG['replacement_shortlist'], rg.size)
                top = np.argpartition(rg.ravel(), -count)[-count:]
                top = sorted(top, key=lambda x: (-int(rg.ravel()[x]), int(x)))
                for f in top:
                    ii, jj = np.unravel_index(f, rg.shape)
                    gain = int(rg[ii, jj])
                    if gain <= max(0, best[0]):
                        break
                    u, v = int(chosen[ii]), int(outside[jj])
                    trial = s.copy()
                    trial.swap(u, v)
                    if trial.coverage() + 1e-8 >= floor:
                        best = gain, u, v
                        break
        if best[0] <= 0:
            break
        s.swap(best[1], best[2])
        moves += 1
    return moves


def dynamic(s, floor, seed):
    rng = np.random.default_rng(seed)
    accepted = 0
    for _ in range(CONFIG['escape_rounds']):
        found = False
        for length in CONFIG['lengths']:
            best = None
            for _ in range(CONFIG['proposals_per_length']):
                trial = s.copy()
                for _ in range(length):
                    chosen = np.flatnonzero(trial.roles)
                    u = int(rng.choice(chosen))
                    v = int(rng.choice(np.flatnonzero(trial.roles != trial.roles[u])))
                    trial.swap(u, v)
                if trial.score > s.score and trial.coverage() + 1e-8 >= floor:
                    if best is None or trial.score > best.score:
                        best = trial
            if best is not None:
                s = best
                descent(s, floor)
                accepted += 1
                found = True
                break
        if not found:
            break
    return s, accepted


def lock():
    OUT.mkdir(exist_ok=True)
    paths = [Path(__file__), ROOT / 'balanced_assay_design.md',
             ROOT / 'balanced_assay_implementation_lock.md',
             ROOT / 'YeRI_supplementary_data_15.xlsx',
             ROOT / 'YDIP/interactions.txt', ROOT / 'YBIOGRID/interactions.txt']
    payload = dict(config=CONFIG, files={str(p.relative_to(ROOT.parent)): sha(p) for p in paths},
                   python=platform.python_version(), numpy=np.__version__, openpyxl=openpyxl.__version__)
    target = OUT / 'lock.json'
    if target.exists():
        assert json.loads(target.read_text()) == payload, 'Frozen inputs/code changed; do not reuse results.'
    else:
        target.write_text(json.dumps(payload, indent=2), encoding='utf-8')


def run():
    lock()
    names, source_edges = inputs()
    for source, edges in source_edges.items():
        for seed in CONFIG['seeds']:
            target = OUT / f'{source}_{seed}.json'
            if target.exists():
                continue
            permutation = np.random.default_rng(seed).permutation(len(edges))
            cut = int(.8 * len(edges))
            train, test = edges[permutation[:cut]], edges[permutation[cut:]]
            a = matrix(len(names), train)
            assert not np.any(a[test[:, 0], test[:, 1]])
            started = time.perf_counter()
            order = coverage_order(a, 200)
            cov_time = time.perf_counter() - started
            degree = np.argsort(-a.sum(axis=1).astype(np.int32), kind='stable')
            rows, panels = [], {}
            for k in CONFIG['budgets']:
                baseline = {}
                for method in ('coverage_partition', 'degree_partition', 'role_greedy'):
                    t = time.perf_counter()
                    s = role_greedy(a, k) if method == 'role_greedy' else partition(a, order[:k] if method == 'coverage_partition' else degree[:k])
                    elapsed = time.perf_counter() - t + (cov_time if method == 'coverage_partition' else 0.)
                    baseline[method] = s, elapsed
                reference = baseline['coverage_partition'][0].coverage()
                for guard in CONFIG['guards']:
                    floor = guard * reference
                    states = dict(baseline)
                    t = time.perf_counter()
                    single = baseline['coverage_partition'][0].copy()
                    n_single = descent(single, floor)
                    single_time = baseline['coverage_partition'][1] + time.perf_counter() - t
                    states['single_swap'] = single, single_time
                    t = time.perf_counter()
                    candidate, escapes = dynamic(single.copy(), floor, seed + k)
                    states['dynamic_step'] = candidate, single_time + time.perf_counter() - t
                    for method, (s, elapsed) in states.items():
                        key = f'{k}/{guard}/{method}'
                        count = int(np.count_nonzero(((s.roles[test[:, 0]] == 1) & (s.roles[test[:, 1]] == 2)) | ((s.roles[test[:, 0]] == 2) & (s.roles[test[:, 1]] == 1))))
                        assert np.count_nonzero(s.roles == 1) == np.count_nonzero(s.roles == 2) == k // 2
                        rows.append(dict(source=source, seed=seed, budget=k, guard=guard, method=method,
                                         train_edges=len(train), test_edges=len(test), train_cross=s.score,
                                         test_cross=count, test_fraction=count / len(test), coverage=s.coverage(),
                                         reference_coverage=reference, feasible=s.coverage() + 1e-8 >= floor,
                                         seconds=elapsed, single_moves=n_single if method in ('single_swap', 'dynamic_step') else 0,
                                         escapes=escapes if method == 'dynamic_step' else 0))
                        panels[key] = {role: [names[i] for i in np.flatnonzero(s.roles == code)] for role, code in [('DB', 1), ('AD', 2)]}
            target.write_text(json.dumps(dict(rows=rows, panels=panels, train_sha256=hashlib.sha256(train.tobytes()).hexdigest(), test_sha256=hashlib.sha256(test.tobytes()).hexdigest()), indent=2), encoding='utf-8')
            print(f'Finished {source} seed={seed}', flush=True)
    # Full-graph timing is separate from internal test-edge scores.
    timing = []
    for source, edges in source_edges.items():
        a = matrix(len(names), edges)
        for k in CONFIG['budgets']:
            for repeat in range(3):
                t = time.perf_counter()
                role_greedy(a, k)
                baseline_time = time.perf_counter() - t
                t = time.perf_counter()
                s = partition(a, coverage_order(a, k))
                floor = .98 * s.coverage()
                descent(s, floor)
                dynamic(s, floor, CONFIG['seeds'][repeat] + k)
                timing.append(dict(source=source, budget=k, repeat=repeat, baseline_seconds=baseline_time, candidate_seconds=time.perf_counter() - t))
    (OUT / 'full_graph_timing.json').write_text(json.dumps(timing, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--lock-only', action='store_true')
    args = parser.parse_args()
    lock() if args.lock_only else run()
