"""Independent simulated-annealing comparison at the archived primary budget.

This is a new, separate experiment. It does not modify the archived algorithms,
their locks, or the original sealed results. Generation never opens evaluation
edges; scoring is a separate command that requires a complete panel seal.
"""

import hashlib
import json
import math
import os
import platform
import sys
import time
from pathlib import Path

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
import openpyxl

import balanced_assay_pilot as pilot
import balanced_assay_external as external
import balanced_assay_initialization_time as frozen


REPO = Path(__file__).resolve().parents[1]
INPUTS = REPO / "data" / "inputs"
ARCHIVE = REPO / "reproducibility" / "tmp" / "ppi_data_audit"
CALIBRATION = ARCHIVE / "balanced_assay_initialization_time_20261024"
OUTPUT = REPO / "results" / "independent_sa"
SOURCES = ("YDIP", "YBIOGRID", "SpaceI")
SEEDS = tuple(range(20261020, 20261040))
METHODS = ("A", "B", "SA")
MULTIPLIER = 5
TEMPERATURE_HIGH = 4.0
TEMPERATURE_LOW = 0.2
ANNEAL_CYCLES = 4
DEGREE_POOL = 1000
TIMER_RESERVE = 0.01
TOLERANCE_FRACTION = 0.05
TOLERANCE_SECONDS = 0.02

pilot.ROOT = INPUTS
external.INPUTS = {
    "candidates": INPUTS / "balanced_assay_space1_20261024" / "historical_mapping_candidates.json",
    "train": INPUTS / "balanced_assay_source_audit" / "raw" / "HI-II-14.tsv",
    "outcome": INPUTS / "HuRI.tsv",
    "earlier": INPUTS / "balanced_assay_source_audit" / "raw" / "H-I-05.tsv",
    "calibration": INPUTS / "balanced_assay_source_audit" / "raw" / "Test_space_screens-19.tsv",
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".pending")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def environment():
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "openpyxl": openpyxl.__version__,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "threads": {key: os.environ[key] for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
    }


def config():
    return {
        "sources": SOURCES,
        "seeds": SEEDS,
        "methods": METHODS,
        "multiplier": MULTIPLIER,
        "temperature_high": TEMPERATURE_HIGH,
        "temperature_low": TEMPERATURE_LOW,
        "anneal_cycles": ANNEAL_CYCLES,
        "degree_pool": DEGREE_POOL,
        "timer_reserve": TIMER_RESERVE,
        "tolerance_fraction": TOLERANCE_FRACTION,
        "tolerance_seconds": TOLERANCE_SECONDS,
    }


def lock(create=False):
    files = [
        Path(__file__),
        OUTPUT / "protocol.md",
        *(REPO / "src" / name for name in (
            "balanced_assay_pilot.py",
            "balanced_assay_tradeoff.py",
            "balanced_assay_timecheck.py",
            "balanced_assay_initialization_time.py",
            "balanced_assay_external.py",
        )),
        INPUTS / "YeRI_supplementary_data_15.xlsx",
        INPUTS / "YDIP" / "interactions.txt",
        INPUTS / "YBIOGRID" / "interactions.txt",
        *external.INPUTS.values(),
    ]
    value = {
        "config": config(),
        "environment": environment(),
        "files_sha256": {str(path.relative_to(REPO)): sha(path) for path in files},
    }
    path = OUTPUT / "lock.json"
    if create:
        assert not path.exists(), "Lock already exists"
        write_json(path, value)
    else:
        assert json.loads(path.read_text(encoding="utf-8")) == json.loads(json.dumps(value)), "Protocol, inputs, code, or environment changed"
    return value


def yeast_inputs():
    return pilot.inputs()


def instance(source, seed, yeast_names, yeast_edges):
    if source == "SpaceI":
        names = external.cohort()
        train = external.restricted(external.edges(external.INPUTS["train"]), names)
        key = source
    else:
        names = yeast_names
        edges = yeast_edges[source]
        order = np.random.default_rng(seed).permutation(len(edges))
        train = edges[order[:int(0.8 * len(edges))]]
        key = f"{source}_{seed}"
    calibration = json.loads((CALIBRATION / f"instance_{key}.json").read_text(encoding="utf-8"))
    assert names == calibration["names"]
    assert hashlib.sha256(train.tobytes()).hexdigest() == calibration["train_sha256"]
    return key, names, train, calibration


def validate_panel(state, floor):
    db = np.flatnonzero(state.roles == 1)
    ad = np.flatnonzero(state.roles == 2)
    assert len(db) == len(ad) == 100
    assert len(np.intersect1d(db, ad)) == 0
    coverage = state.coverage()
    assert coverage + 1e-7 >= floor
    assert state.score == int(state.a[np.ix_(db, ad)].sum())
    return db.tolist(), ad.tolist(), float(coverage)


def random_external(rng, roles, degree_order):
    if rng.random() < 0.5:
        pool = degree_order[:min(DEGREE_POOL, len(degree_order))]
        for _ in range(16):
            candidate = int(pool[rng.integers(len(pool))])
            if roles[candidate] == 0:
                return candidate
    while True:
        candidate = int(rng.integers(len(roles)))
        if roles[candidate] == 0:
            return candidate


def anneal(a, floor, budget, seed):
    started = time.perf_counter()
    deadline = started + max(0.0, budget - TIMER_RESERVE)
    selected = pilot.coverage_order(a, 200)
    state = pilot.partition(a, selected)
    best = state.copy()
    rng = np.random.default_rng(seed)
    degree_order = np.argsort(-a.sum(axis=1).astype(np.int32), kind="stable")
    role_ids = [np.flatnonzero(state.roles == role) for role in (1, 2)]
    proposals = accepted = feasible = improvements = 0
    initialization_seconds = time.perf_counter() - started
    while time.perf_counter() < deadline:
        fraction = min(0.999999, (time.perf_counter() - started) / budget)
        cycle_position = (fraction * ANNEAL_CYCLES) % 1.0
        temperature = TEMPERATURE_HIGH * (TEMPERATURE_LOW / TEMPERATURE_HIGH) ** cycle_position
        proposals += 1
        if rng.random() < 0.5:
            left = int(rng.integers(len(role_ids[0])))
            right = int(rng.integers(len(role_ids[1])))
            u, v = int(role_ids[0][left]), int(role_ids[1][right])
            replacement = False
        else:
            side = int(rng.integers(2))
            index = int(rng.integers(len(role_ids[side])))
            u = int(role_ids[side][index])
            v = random_external(rng, state.roles, degree_order)
            replacement = True
        gain = state.delta(u, v)
        if gain < 0 and rng.random() >= math.exp(gain / temperature):
            continue
        if replacement:
            trial = state.copy()
            trial.swap(u, v)
            if trial.coverage() + 1e-7 < floor:
                continue
            feasible += 1
            if time.perf_counter() >= deadline:
                break
            state = trial
            role_ids[side][index] = v
        else:
            feasible += 1
            if time.perf_counter() >= deadline:
                break
            state.swap(u, v)
            role_ids[0][left], role_ids[1][right] = v, u
        accepted += 1
        if state.score > best.score:
            best = state.copy()
            improvements += 1
    elapsed = time.perf_counter() - started
    return best, {
        "initialization_seconds": initialization_seconds,
        "total_seconds": elapsed,
        "allowed_seconds": budget,
        "timing_valid": elapsed <= budget + max(TOLERANCE_SECONDS, TOLERANCE_FRACTION * budget),
        "proposals": proposals,
        "feasible_proposals": feasible,
        "accepted": accepted,
        "best_improvements": improvements,
    }


def method_order(source_index, seed_index):
    shift = (source_index + seed_index) % len(METHODS)
    return METHODS[shift:] + METHODS[:shift]


def expected_names():
    return [f"panel_{source}_{seed}_{method}.json" for seed in SEEDS for source in SOURCES for method in METHODS]


def generate():
    lock()
    assert not (OUTPUT / "seal.json").exists(), "Outputs already sealed"
    yeast_names, yeast_edges = yeast_inputs()
    done = sum((OUTPUT / name).exists() for name in expected_names())
    for seed_index, seed in enumerate(SEEDS):
        for source_index, source in enumerate(SOURCES):
            key, names, train, cal = instance(source, seed, yeast_names, yeast_edges)
            a = pilot.matrix(len(names), train)
            budget = MULTIPLIER * cal["t0_seconds"]
            for method in method_order(source_index, seed_index):
                path = OUTPUT / f"panel_{source}_{seed}_{method}.json"
                if path.exists():
                    continue
                print(f"START {done+1}/180 {source} {seed} {method} {budget:.3f}s", flush=True)
                if method == "SA":
                    state, timing = anneal(a, cal["floor"], budget, seed)
                else:
                    state, timing = frozen.pipeline(a, method, cal["floor"], budget, seed)
                db, ad, coverage = validate_panel(state, cal["floor"])
                row = {
                    "source": source, "seed": seed, "method": method, "instance": key,
                    "train_sha256": cal["train_sha256"], "db": db, "ad": ad,
                    "train_cross": int(state.score), "coverage": coverage,
                    **timing,
                }
                write_json(path, row)
                done += 1
                print(f"DONE {done}/180 valid={timing['timing_valid']} train={state.score}", flush=True)
    assert done == len(expected_names())
    value = {
        "lock_sha256": sha(OUTPUT / "lock.json"),
        "panels": done,
        "files_sha256": {name: sha(OUTPUT / name) for name in expected_names()},
    }
    write_json(OUTPUT / "seal.json", value)
    print("SEALED: all panels written before any evaluation edges were read", flush=True)


def score():
    lock()
    seal = json.loads((OUTPUT / "seal.json").read_text(encoding="utf-8"))
    assert seal["lock_sha256"] == sha(OUTPUT / "lock.json")
    assert set(seal["files_sha256"]) == set(expected_names())
    assert all(sha(OUTPUT / name) == value for name, value in seal["files_sha256"].items())
    yeast_names, yeast_edges = yeast_inputs()
    outcomes = external.novel(
        external.edges(external.INPUTS["outcome"]),
        [external.edges(external.INPUTS[key]) for key in external.CONFIG["exclusion"]],
    )
    rows = []
    for seed in SEEDS:
        for source in SOURCES:
            key, names, train, cal = instance(source, seed, yeast_names, yeast_edges)
            if source == "SpaceI":
                test = external.restricted(outcomes, names)
            else:
                edges = yeast_edges[source]
                order = np.random.default_rng(seed).permutation(len(edges))
                test = edges[order[int(0.8 * len(edges)):]]
                historical = json.loads((ARCHIVE / "balanced_assay_timecheck" / f"{source}_{seed}.json").read_text(encoding="utf-8"))
                assert hashlib.sha256(test.tobytes()).hexdigest() == historical["test_sha256"]
            a = pilot.matrix(len(names), train)
            for method in METHODS:
                panel = json.loads((OUTPUT / f"panel_{source}_{seed}_{method}.json").read_text(encoding="utf-8"))
                state = pilot.State(a, panel["db"], panel["ad"])
                assert state.score == panel["train_cross"]
                db, ad, coverage = validate_panel(state, cal["floor"])
                assert abs(coverage - panel["coverage"]) < 1e-6
                roles = state.roles
                test_cross = external.cross(roles, test)
                rows.append({
                    "source": source, "seed": seed, "method": method,
                    "train_cross": panel["train_cross"], "test_cross": test_cross,
                    "coverage_ratio": coverage / cal["reference_coverage"],
                    "total_seconds": panel["total_seconds"],
                    "allowed_seconds": panel["allowed_seconds"],
                    "timing_valid": panel["timing_valid"],
                })
    write_json(OUTPUT / "scores.json", {"rows": rows, "seal_sha256": sha(OUTPUT / "seal.json")})
    print("SCORED", len(rows), "sealed panels", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("lock", "generate", "score"):
        raise SystemExit("Usage: python independent_sa_benchmark.py lock|generate|score")
    {"lock": lambda: lock(True), "generate": generate, "score": score}[sys.argv[1]]()
