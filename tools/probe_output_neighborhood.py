"""Post hoc local replacement-feasibility probe at sealed output panels."""

import json
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np

import sys

repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo / "src"))
import coverage_sensitivity_benchmark as sensitivity
import independent_sa_benchmark as earlier


SAMPLES = 1000
BLOCK = 64
out = repo / "results" / "coverage_sensitivity"
sensitivity.lock()
sensitivity.verify_seals()


def feasible_fraction(state, floor, seed):
    a = state.a
    selected = np.flatnonzero(state.roles)
    outside = np.flatnonzero(state.roles == 0)
    rng = np.random.default_rng(seed)
    u = selected[rng.integers(len(selected), size=SAMPLES)]
    v = outside[rng.integers(len(outside), size=SAMPLES)]
    counts = state.cd + state.cp
    is_selected = state.roles != 0
    feasible = 0
    for start in range(0, SAMPLES, BLOCK):
        us, vs = u[start:start + BLOCK], v[start:start + BLOCK]
        changed = counts[:, None] + a[:, vs].astype(np.int16) - a[:, us].astype(np.int16)
        power = np.power(0.2, changed)
        coverage = np.where(is_selected[:, None], 1.0, 1.0 - power).sum(axis=0)
        col = np.arange(len(us))
        coverage += -power[us, col] + power[vs, col]
        feasible += int(np.count_nonzero(coverage + 1e-8 >= floor))
        if start == 0:
            for i in range(min(3, len(us))):
                trial = state.copy()
                trial.swap(int(us[i]), int(vs[i]))
                assert abs(float(coverage[i]) - trial.coverage()) < 1e-6
    return feasible / SAMPLES


yeast_names, yeast_edges = earlier.yeast_inputs()
records = []
for source in earlier.SOURCES:
    for seed in earlier.SEEDS:
        key, names, train, cal = earlier.instance(source, seed, yeast_names, yeast_edges)
        a = earlier.pilot.matrix(len(names), train)
        for guard in sensitivity.GUARDS_ALL:
            floor = guard * cal["reference_coverage"]
            for method in sensitivity.METHODS:
                path = ((earlier.OUTPUT / f"panel_{source}_{seed}_{method}.json") if guard == 0.98
                        else (out / sensitivity.panel_name(source, seed, guard, method)))
                panel = json.loads(path.read_text(encoding="utf-8"))
                state = earlier.pilot.State(a, panel["db"], panel["ad"])
                earlier.validate_panel(state, floor)
                fraction = feasible_fraction(state, floor, seed + 2000000 + earlier.SOURCES.index(source))
                records.append({
                    "source": source, "seed": seed, "tau": guard, "method": method,
                    "feasible_replacement_fraction": fraction,
                })
        print(f"PROBED {source} {seed}", flush=True)

summary = defaultdict(dict)
for source in earlier.SOURCES:
    for guard in sensitivity.GUARDS_ALL:
        tag = sensitivity.tag(guard)
        summary[source][tag] = {}
        for method in sensitivity.METHODS:
            values = [r["feasible_replacement_fraction"] for r in records
                      if r["source"] == source and r["tau"] == guard and r["method"] == method]
            assert len(values) == 20
            summary[source][tag][method] = {
                "median": statistics.median(values),
                "range": [min(values), max(values)],
            }

payload = {
    "samples_per_panel": SAMPLES,
    "interpretation": "Post hoc fixed-output local replacement probe; output panels differ across thresholds and methods; not a timed search-space count",
    "new_seal_sha256": earlier.sha(out / "seal.json"),
    "prior_seal_sha256": earlier.sha(earlier.OUTPUT / "seal.json"),
    "records": records,
    "summary": summary,
}
(out / "output_neighborhood.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))
