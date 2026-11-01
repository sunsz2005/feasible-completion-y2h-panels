"""Verify the archived run and its browsable code, data, and result copies."""

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
ARCHIVE = ROOT / "reproducibility"
SCORES = ARCHIVE / "tmp/ppi_data_audit/balanced_assay_initialization_time_20261024/scores.json"
EXPORTS = {
    "summary": ROOT / "data/primary_summary.csv",
    "paired": ROOT / "data/primary_paired.csv",
    "rows": ROOT / "data/primary_runs.csv",
}
EXCLUDED_INPUTS = {"YDIP/interactions.txt", "YBIOGRID/interactions.txt"}
EXCLUDED_ARCHIVE = {"tmp/ppi_data_audit/" + name for name in EXCLUDED_INPUTS}


def main():
    manifest = json.loads((ARCHIVE / "manifest_sha256.json").read_text(encoding="utf-8"))
    missing_excluded = set()
    for name, expected in manifest.items():
        path = ARCHIVE / name
        if name in EXCLUDED_ARCHIVE and not path.exists():
            missing_excluded.add(name)
            continue
        assert path.is_file(), f"Missing archived file: {name}"
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == expected, f"Archive hash mismatch: {name}"
    assert missing_excluded in (set(), EXCLUDED_ARCHIVE), "Only one excluded yeast input is present"

    code_names = (
        "balanced_assay_ep_adapter.py", "balanced_assay_external.py",
        "balanced_assay_gate_control.py", "balanced_assay_initialization_time.py",
        "balanced_assay_pilot.py", "balanced_assay_timecheck.py",
        "balanced_assay_tradeoff.py",
    )
    for name in code_names:
        assert (ROOT / "src" / name).read_bytes() == (ARCHIVE / "tmp" / name).read_bytes(), name
    for name in ("audit_primary_comparison.py", "test_balanced_assay_initialization_time.py"):
        assert (ROOT / "tests" / name).read_bytes() == (ARCHIVE / "tmp" / name).read_bytes(), name
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
        assert (ROOT / "data/inputs" / name).read_bytes() == (ARCHIVE / "tmp/ppi_data_audit" / name).read_bytes(), name
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
    print(f"OK: {len(manifest) - len(missing_excluded)} archived files and browsable copies; "
          f"{len(missing_excluded)} excluded yeast inputs; 360 panels; 360 runs; 180 pairs; 9 summaries")


if __name__ == "__main__":
    main()
