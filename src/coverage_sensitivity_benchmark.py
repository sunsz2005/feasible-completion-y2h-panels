"""Replay sealed A/B sensitivity results or generate a new two-threshold batch."""

import json
import os
import sys
from pathlib import Path

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np

import balanced_assay_initialization_time as frozen
import balanced_assay_external as external
import coverage_sensitivity_inputs as shared


REPO = Path(__file__).resolve().parents[1]
ARCHIVED = REPO / "results" / "coverage_sensitivity"
OUTPUT = ARCHIVED
REPORT = REPO / "reproducibility" / "output" / "coverage_sensitivity"
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
            for source in shared.SOURCES for seed in shared.SEEDS
            for guard in GUARDS_NEW for method in METHODS]


def lock(create=False):
    shared.lock()
    if OUTPUT == ARCHIVED:
        assert not create, "Original archived results cannot be relocked"
        return shared.provenance()
    path = OUTPUT / "lock.json"
    value = {
        "purpose": "New timed A/B sensitivity rerun, not the original experiment",
        "config": shared.provenance()["config"],
        "environment": shared.environment(),
        "files_sha256": {p.name:shared.sha(p) for p in sorted((REPO/"src").glob("*.py"))},
        "original_provenance_sha256": shared.sha(ARCHIVED / "provenance.json"),
    }
    if create:
        assert not path.exists(), "New lock already exists"
        shared.write_json(path, value)
    else:
        assert json.loads(path.read_text(encoding="utf-8")) == value, "New run configuration changed"
    return value


def seal_hashes():
    originals = shared.provenance()["original_seals"]
    return (originals["new"] if OUTPUT == ARCHIVED else shared.sha(OUTPUT / "seal.json"), originals["reference"])


def method_order(source_index, seed_index, guard_index):
    return METHODS if (source_index + seed_index + guard_index) % 2 == 0 else METHODS[::-1]


def generate():
    assert OUTPUT != ARCHIVED, "Use --output for a fresh timed run"
    lock()
    assert not (OUTPUT / "seal.json").exists(), "Outputs already sealed"
    yeast_names, yeast_edges = shared.yeast_inputs()
    done = sum((OUTPUT / name).exists() for name in expected_names())
    for seed_index, seed in enumerate(shared.SEEDS):
        for source_index, source in enumerate(shared.SOURCES):
            key, names, train, cal = shared.instance(source, seed, yeast_names, yeast_edges)
            a = shared.pilot.matrix(len(names), train)
            budget = shared.MULTIPLIER * cal["t0_seconds"]
            for guard_index, guard in enumerate(GUARDS_NEW):
                floor = guard * cal["reference_coverage"]
                for method in method_order(source_index, seed_index, guard_index):
                    name = panel_name(source, seed, guard, method)
                    path = OUTPUT / name
                    if path.exists():
                        continue
                    print(f"START {done+1}/240 {source} {seed} tau={tag(guard)} {method}", flush=True)
                    state, timing = frozen.pipeline(a, method, floor, budget, seed)
                    db, ad, coverage = shared.validate_panel(state, floor)
                    shared.write_json(path, {
                        "source": source, "seed": seed, "method": method,
                        "tau": guard, "instance": key, "train_sha256": cal["train_sha256"],
                        "reference_coverage": cal["reference_coverage"],
                        "db": db, "ad": ad, "train_cross": int(state.score),
                        "coverage": coverage, **timing,
                    })
                    done += 1
                    print(f"DONE {done}/240 valid={timing['timing_valid']} train={state.score}", flush=True)
    assert done == len(expected_names())
    shared.write_json(OUTPUT / "seal.json", {
        "lock_sha256": shared.sha(OUTPUT / "lock.json"),
        "panels": done,
        "files_sha256": {name: shared.sha(OUTPUT / name) for name in expected_names()},
    })
    print("SEALED: 240 panels; evaluation edges not yet read", flush=True)


def verify_seals():
    hashes = shared.provenance()["original_panel_hashes"]
    assert len(hashes) == 360
    for rel, digest in hashes.items():
        assert shared.sha(ARCHIVED / rel) == digest, rel
    if OUTPUT != ARCHIVED:
        seal = json.loads((OUTPUT / "seal.json").read_text(encoding="utf-8"))
        assert seal["lock_sha256"] == shared.sha(OUTPUT / "lock.json")
        assert set(seal["files_sha256"]) == set(expected_names())
        for rel, digest in seal["files_sha256"].items():
            assert shared.sha(OUTPUT / rel) == digest, rel


def sample_replacement_feasibility(a, reference, seed, source):
    selected = np.flatnonzero(reference.roles)
    outside = np.flatnonzero(reference.roles == 0)
    rng = np.random.default_rng(seed + 1000000 + shared.SOURCES.index(source))
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
    yeast_names, yeast_edges = shared.yeast_inputs()
    outcomes = external.novel(
        external.edges(external.INPUTS["outcome"]),
        [external.edges(external.INPUTS[key]) for key in external.CONFIG["exclusion"]],
    )
    rows = []
    probes = []
    for seed in shared.SEEDS:
        for source in shared.SOURCES:
            key, names, train, cal = shared.instance(source, seed, yeast_names, yeast_edges)
            if source == "SpaceI":
                test = external.restricted(outcomes, names)
            else:
                edges = yeast_edges[source]
                order = np.random.default_rng(seed).permutation(len(edges))
                test = edges[order[int(0.8 * len(edges)):]]
            a = shared.pilot.matrix(len(names), train)
            selected = shared.pilot.coverage_order(a, 200)
            reference = shared.pilot.partition(a, selected)
            assert abs(reference.coverage() - cal["reference_coverage"]) < 1e-6
            probes.append({"source": source, "seed": seed,
                           "feasible_replacement_fraction": sample_replacement_feasibility(a, reference, seed, source)})
            for guard in GUARDS_ALL:
                floor = guard * cal["reference_coverage"]
                for method in METHODS:
                    path = ((shared.OUTPUT / f"panel_{source}_{seed}_{method}.json") if guard == 0.98
                            else (OUTPUT / panel_name(source, seed, guard, method)))
                    panel = json.loads(path.read_text(encoding="utf-8"))
                    state = shared.pilot.State(a, panel["db"], panel["ad"])
                    assert state.score == panel["train_cross"]
                    _, _, coverage = shared.validate_panel(state, floor)
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
    shared.write_json(REPORT / "scores.json", {
        "rows": rows, "probes": probes,
        "new_seal_sha256": seal_hashes()[0],
        "prior_seal_sha256": seal_hashes()[1],
    })
    print("SCORED", len(rows), "sealed panels", flush=True)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("lock", "generate", "score"))
    parser.add_argument("--output", type=Path, help="New timed-run directory; never the archived results")
    args = parser.parse_args()
    if args.output:
        OUTPUT = args.output.resolve()
        forbidden = [REPO/"src", REPO/"data", REPO/"results", REPO/"reproducibility"]
        assert OUTPUT != REPO and not any(OUTPUT == p or p in OUTPUT.parents for p in forbidden), "Choose a new directory outside preserved inputs/results"
        OUTPUT.mkdir(parents=True, exist_ok=True)
        REPORT = OUTPUT / "replay"
    {"lock": lambda: lock(True), "generate": generate, "score": score}[args.command]()
