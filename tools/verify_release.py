"""Verify the archived run and its browsable code, data, and result copies."""

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
ARCHIVE = ROOT / "reproducibility"
ORIGINAL = ARCHIVE / "primary_archive"
SCORES = ORIGINAL / "ppi_data_audit/balanced_assay_initialization_time_20261024/scores.json"
EXPORTS = {
    "summary": ROOT / "data/primary_summary.csv",
    "paired": ROOT / "data/primary_paired.csv",
    "rows": ROOT / "data/primary_runs.csv",
}
EXCLUDED_INPUTS = {"YDIP/interactions.txt", "YBIOGRID/interactions.txt"}
EXCLUDED_ARCHIVE = {"tmp/ppi_data_audit/" + name for name in EXCLUDED_INPUTS}


def main():
    manifest = json.loads((ARCHIVE / "manifest_sha256.json").read_text(encoding="utf-8"))
    omitted = json.loads((ROOT / "results/primary/replay_metadata.json").read_text(encoding="utf-8"))
    missing_locks = set()
    missing_excluded = set()
    for name, expected in manifest.items():
        path = (ORIGINAL / name.removeprefix("tmp/") if name.startswith("tmp/")
                else ARCHIVE / "primary_audits" / name.removeprefix("output/") if name.startswith("output/")
                else ARCHIVE / name)
        relative = path.relative_to(ROOT).as_posix()
        if relative in omitted["omitted_locks_sha256"]:
            assert not path.exists(), f"Historical lock should be absent: {relative}"
            assert omitted["omitted_locks_sha256"][relative] == expected
            missing_locks.add(relative)
            continue
        if name in EXCLUDED_ARCHIVE and not path.exists():
            missing_excluded.add(name)
            continue
        assert path.is_file(), f"Missing archived file: {name}"
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == expected, f"Archive hash mismatch: {name}"
    assert missing_locks == set(omitted["omitted_locks_sha256"])
    for rel, expected in omitted["retained_dependency_sha256"].items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == expected, rel
    assert missing_excluded in (set(), EXCLUDED_ARCHIVE), "Only one excluded yeast input is present"

    code_names = (
        "balanced_assay_ep_adapter.py", "balanced_assay_external.py",
        "balanced_assay_gate_control.py", "balanced_assay_initialization_time.py",
        "balanced_assay_pilot.py", "balanced_assay_timecheck.py",
        "balanced_assay_tradeoff.py",
    )
    for name in code_names:
        assert (ROOT / "src" / name).read_bytes() == (ORIGINAL / name).read_bytes(), name
    for name in ("audit_primary_comparison.py", "test_balanced_assay_initialization_time.py"):
        assert (ROOT / "tests" / name).read_bytes() == (ORIGINAL / name).read_bytes(), name
    inputs = (
        "HuRI.tsv", "YBIOGRID/interactions.txt", "YDIP/interactions.txt",
        "YeRI_supplementary_data_15.xlsx",
        "balanced_assay_source_audit/raw/HI-II-14.tsv",
        "balanced_assay_source_audit/raw/H-I-05.tsv",
        "balanced_assay_source_audit/raw/Test_space_screens-19.tsv",
        "balanced_assay_space1_20261024/historical_mapping_candidates.json",
    )
    for name in inputs:
        if name in EXCLUDED_INPUTS and "tmp/ppi_data_audit/" + name in missing_excluded:
            assert not (ROOT / "data/inputs" / name).exists(), f"Unpaired excluded input: {name}"
            continue
        assert (ROOT / "data/inputs" / name).read_bytes() == (ORIGINAL / "ppi_data_audit" / name).read_bytes(), name
    for name in ("scores.json", "sealed_panels.json"):
        assert (ROOT / "results/primary" / name).read_bytes() == (SCORES.parent / name).read_bytes(), name
    assert (ROOT / "requirements.txt").read_bytes() == (ARCHIVE / "requirements.txt").read_bytes()

    scores = json.loads(SCORES.read_text(encoding="utf-8"))
    expected_counts = {"summary": 9, "paired": 180, "rows": 360}
    for key, path in EXPORTS.items():
        source = scores[key]
        assert len(source) == expected_counts[key], f"Unexpected {key} count"
        with path.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            assert reader.fieldnames == list(source[0]), f"CSV columns differ: {path}"
            exported = list(reader)
        assert len(exported) == len(source), f"CSV row count differs: {path}"
        for index, (record, row) in enumerate(zip(source, exported)):
            for field, value in record.items():
                expected = "" if value is None else str(value)
                assert row[field] == expected, f"CSV mismatch: {path}, row {index}, {field}"

    seal = json.loads((SCORES.parent / "sealed_panels.json").read_text(encoding="utf-8"))
    assert seal["panels"] == 360, "Unexpected sealed panel count"
    print(f"OK: {len(manifest) - len(missing_excluded) - len(missing_locks)} archived files and browsable copies; "
          f"{len(missing_excluded)} excluded yeast inputs; {len(missing_locks)} omitted historical locks; 360 panels; 360 runs; 180 pairs; 9 summaries")


def verify_followups():
    distribution = json.loads((ROOT / "distribution_manifest.json").read_text(encoding="utf-8"))
    for rel, expected in distribution["files_sha256"].items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == expected, rel
    sens = ROOT / "results/coverage_sensitivity"
    provenance = json.loads((sens / "provenance.json").read_text(encoding="utf-8"))
    assert len(provenance["original_panel_hashes"]) == 360
    for rel, expected in provenance["original_panel_hashes"].items():
        assert hashlib.sha256((sens / rel).read_bytes()).hexdigest() == expected, rel
    for rel, expected in provenance["input_and_solver_sha256"].items():
        path = ROOT / rel
        if not path.exists() and rel in {"data/inputs/YDIP/interactions.txt", "data/inputs/YBIOGRID/interactions.txt"}:
            continue
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, rel
    scores = json.loads((sens / "scores.json").read_text(encoding="utf-8"))
    assert len(scores["rows"]) == 360
    assert {r["method"] for r in scores["rows"]} == {"A", "B"}
    fly = ROOT / "reproducibility/flybi_validation"
    seal = json.loads((fly / "pre_evaluation_seal_ab.json").read_text(encoding="utf-8"))
    assert len(seal["files"]) == 60
    for rel, expected in seal["files"].items():
        assert hashlib.sha256((fly / rel).read_bytes()).hexdigest() == expected, rel
    scores = json.loads((fly / "scores_ab.json").read_text(encoding="utf-8"))
    assert len(scores["rows"]) == 40 and {r["method"] for r in scores["rows"]} == {"A", "B"}
    print(f"OK: distribution {len(distribution['files_sha256'])} files; 360 sensitivity paths; 40 FlyBi paths and 20 calibrations")


if __name__ == "__main__":
    main()
    verify_followups()
