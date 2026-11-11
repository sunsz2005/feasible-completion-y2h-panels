"""Matched-budget coverage-threshold sensitivity for pipelines A and B.

The 0.98 arm is reused from the separately sealed independent-SA rerun, whose
software, machine, seeds, budget, and A/B implementation match this follow-up.
This script generates and seals the 0.90 and 0.95 arms before evaluation scoring.
"""

import json
import os
import sys
from pathlib import Path

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np

import balanced_assay_initialization_time as frozen
import balanced_assay_external as external
import independent_sa_benchmark as earlier


REPO = Path(__file__).resolve().parents[1]
OUTPUT = REPO / "results" / "coverage_sensitivity"
GUARDS_NEW = (0.90, 0.95)
GUARDS_ALL = (0.90, 0.95, 0.98)
METHODS = ("A", "B")
PROBE_PAIRS = 1000


def tag(guard):
    return f"{guard:.2f}"


def panel_name(source, seed, guard, method):
    return f"panel_{source}_{seed}_tau{tag(guard)}_{method}.json"


def expected_names():
    return [panel_name(source, seed, guard, method)
            for source in earlier.SOURCES for seed in earlier.SEEDS
            for guard in GUARDS_NEW for method in METHODS]


def lock(create=False):
    earlier.lock()
    files = [Path(__file__), OUTPUT / "protocol.md", REPO / "src" / "independent_sa_benchmark.py"]
    value = {
        "config": {
            "sources": earlier.SOURCES, "seeds": earlier.SEEDS,
            "guards_new": GUARDS_NEW, "guards_reused": (0.98,),
            "methods": METHODS, "q": 0.8, "k": 200,
            "multiplier": earlier.MULTIPLIER, "probe_pairs": PROBE_PAIRS,
        },
        "environment": earlier.environment(),
        "files_sha256": {str(path.relative_to(REPO)): earlier.sha(path) for path in files},
        "prior_seal_sha256": earlier.sha(earlier.OUTPUT / "seal.json"),
    }
    path = OUTPUT / "lock.json"
    if create:
        assert not path.exists(), "Lock already exists"
        earlier.write_json(path, value)
    else:
        assert json.loads(path.read_text(encoding="utf-8")) == json.loads(json.dumps(value)), "Protocol or environment changed"
    return value


def method_order(source_index, seed_index, guard_index):
    return METHODS if (source_index + seed_index + guard_index) % 2 == 0 else METHODS[::-1]


def generate():
    lock()
    assert not (OUTPUT / "seal.json").exists(), "Outputs already sealed"
    yeast_names, yeast_edges = earlier.yeast_inputs()
    done = sum((OUTPUT / name).exists() for name in expected_names())
    for seed_index, seed in enumerate(earlier.SEEDS):
        for source_index, source in enumerate(earlier.SOURCES):
            key, names, train, cal = earlier.instance(source, seed, yeast_names, yeast_edges)
            a = earlier.pilot.matrix(len(names), train)
            budget = earlier.MULTIPLIER * cal["t0_seconds"]
            for guard_index, guard in enumerate(GUARDS_NEW):
                floor = guard * cal["reference_coverage"]
                for method in method_order(source_index, seed_index, guard_index):
                    name = panel_name(source, seed, guard, method)
                    path = OUTPUT / name
                    if path.exists():
                        continue
                    print(f"START {done+1}/240 {source} {seed} tau={tag(guard)} {method}", flush=True)
                    state, timing = frozen.pipeline(a, method, floor, budget, seed)
                    db, ad, coverage = earlier.validate_panel(state, floor)
                    earlier.write_json(path, {
                        "source": source, "seed": seed, "method": method,
                        "tau": guard, "instance": key, "train_sha256": cal["train_sha256"],
                        "reference_coverage": cal["reference_coverage"],
                        "db": db, "ad": ad, "train_cross": int(state.score),
                        "coverage": coverage, **timing,
                    })
                    done += 1
                    print(f"DONE {done}/240 valid={timing['timing_valid']} train={state.score}", flush=True)
    assert done == len(expected_names())
    earlier.write_json(OUTPUT / "seal.json", {
        "lock_sha256": earlier.sha(OUTPUT / "lock.json"),
        "panels": done,
        "files_sha256": {name: earlier.sha(OUTPUT / name) for name in expected_names()},
    })
    print("SEALED: 240 panels; evaluation edges not yet read", flush=True)


def verify_seals():
    seal = json.loads((OUTPUT / "seal.json").read_text(encoding="utf-8"))
    assert seal["lock_sha256"] == earlier.sha(OUTPUT / "lock.json")
    assert set(seal["files_sha256"]) == set(expected_names())
    assert all(earlier.sha(OUTPUT / name) == sha for name, sha in seal["files_sha256"].items())
    prior = json.loads((earlier.OUTPUT / "seal.json").read_text(encoding="utf-8"))
    assert prior["lock_sha256"] == earlier.sha(earlier.OUTPUT / "lock.json")
    for seed in earlier.SEEDS:
        for source in earlier.SOURCES:
            for method in METHODS:
                name = f"panel_{source}_{seed}_{method}.json"
                assert earlier.sha(earlier.OUTPUT / name) == prior["files_sha256"][name]


def sample_replacement_feasibility(a, reference, seed, source):
    selected = np.flatnonzero(reference.roles)
    outside = np.flatnonzero(reference.roles == 0)
    rng = np.random.default_rng(seed + 1000000 + earlier.SOURCES.index(source))
    pairs = [(int(selected[rng.integers(len(selected))]), int(outside[rng.integers(len(outside))]))
             for _ in range(PROBE_PAIRS)]
    reference_coverage = reference.coverage()
    counts = {tag(guard): 0 for guard in GUARDS_ALL}
    for u, v in pairs:
        trial = reference.copy()
        trial.swap(u, v)
        ratio = trial.coverage() / reference_coverage
        for guard in GUARDS_ALL:
            counts[tag(guard)] += ratio + 1e-10 >= guard
    return {key: count / PROBE_PAIRS for key, count in counts.items()}


def score():
    lock()
    verify_seals()
    yeast_names, yeast_edges = earlier.yeast_inputs()
    outcomes = external.novel(
        external.edges(external.INPUTS["outcome"]),
        [external.edges(external.INPUTS[key]) for key in external.CONFIG["exclusion"]],
    )
    rows = []
    probes = []
    for seed in earlier.SEEDS:
        for source in earlier.SOURCES:
            key, names, train, cal = earlier.instance(source, seed, yeast_names, yeast_edges)
            if source == "SpaceI":
                test = external.restricted(outcomes, names)
            else:
                edges = yeast_edges[source]
                order = np.random.default_rng(seed).permutation(len(edges))
                test = edges[order[int(0.8 * len(edges)):]]
            a = earlier.pilot.matrix(len(names), train)
            selected = earlier.pilot.coverage_order(a, 200)
            reference = earlier.pilot.partition(a, selected)
            assert abs(reference.coverage() - cal["reference_coverage"]) < 1e-6
            probes.append({"source": source, "seed": seed,
                           "feasible_replacement_fraction": sample_replacement_feasibility(a, reference, seed, source)})
            for guard in GUARDS_ALL:
                floor = guard * cal["reference_coverage"]
                for method in METHODS:
                    path = ((earlier.OUTPUT / f"panel_{source}_{seed}_{method}.json") if guard == 0.98
                            else (OUTPUT / panel_name(source, seed, guard, method)))
                    panel = json.loads(path.read_text(encoding="utf-8"))
                    state = earlier.pilot.State(a, panel["db"], panel["ad"])
                    assert state.score == panel["train_cross"]
                    _, _, coverage = earlier.validate_panel(state, floor)
                    assert abs(coverage - panel["coverage"]) < 1e-6
                    rows.append({
                        "source": source, "seed": seed, "tau": guard,
                        "method": method, "train_cross": state.score,
                        "test_cross": external.cross(state.roles, test),
                        "coverage_ratio": coverage / cal["reference_coverage"],
                        "total_seconds": panel["total_seconds"],
                        "allowed_seconds": panel["allowed_seconds"],
                        "timing_valid": panel["timing_valid"],
                    })
    earlier.write_json(OUTPUT / "scores.json", {
        "rows": rows, "probes": probes,
        "new_seal_sha256": earlier.sha(OUTPUT / "seal.json"),
        "prior_seal_sha256": earlier.sha(earlier.OUTPUT / "seal.json"),
    })
    print("SCORED", len(rows), "sealed panels", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("lock", "generate", "score"):
        raise SystemExit("Usage: python coverage_sensitivity_benchmark.py lock|generate|score")
    {"lock": lambda: lock(True), "generate": generate, "score": score}[sys.argv[1]]()
