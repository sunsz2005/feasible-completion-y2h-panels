"""Inputs and calibration shared by the manuscript's A/B sensitivity paths.

The numerical input and validation functions are extracted from the historical
workflow. Archived panel bytes and original hash entries are preserved separately.
"""
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
import openpyxl
import balanced_assay_pilot as pilot
import balanced_assay_external as external

REPO = Path(__file__).resolve().parents[1]
INPUTS = REPO / "data" / "inputs"
CALIBRATION = REPO / "reproducibility" / "primary_archive" / "ppi_data_audit" / "balanced_assay_initialization_time_20261024"
OUTPUT = REPO / "results" / "coverage_sensitivity" / "reference_batch"
SOURCES = ("YDIP", "YBIOGRID", "SpaceI")
SEEDS = tuple(range(20261020, 20261040))
MULTIPLIER = 5
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


def provenance():
    return json.loads((REPO / "results/coverage_sensitivity/provenance.json").read_text(encoding="utf-8"))

def lock():
    """Verify retained input/solver hashes; do not impersonate the original lock."""
    for rel, expected in provenance()["input_and_solver_sha256"].items():
        path = REPO / rel
        if not path.is_file():
            raise FileNotFoundError(f"{rel}: restore excluded inputs with tools/restore_yeast_inputs.py")
        assert sha(path) == expected, rel

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


