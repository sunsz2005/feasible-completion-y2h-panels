"""Summarize the sealed independent-SA comparison without selecting runs."""

import csv
import json
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "results" / "independent_sa"
data = json.loads((ROOT / "scores.json").read_text(encoding="utf-8"))
rows = data["rows"]
assert len(rows) == 180
sources = ("YDIP", "YBIOGRID", "SpaceI")
methods = ("A", "B", "SA")
comparisons = (("B", "SA"), ("B", "A"), ("SA", "A"))


def median(values):
    return statistics.median(values)


summary = {"seal_sha256": data["seal_sha256"], "sources": {}}
paired_rows = []
for source in sources:
    records = [row for row in rows if row["source"] == source]
    by_seed = {}
    for row in records:
        by_seed.setdefault(row["seed"], {})[row["method"]] = row
    assert len(by_seed) == 20
    assert all(set(group) == set(methods) for group in by_seed.values())
    source_summary = {"methods": {}, "comparisons": {}}
    for method in methods:
        group = [by_seed[seed][method] for seed in sorted(by_seed)]
        source_summary["methods"][method] = {
            "n": len(group),
            "evaluation_median": median([row["test_cross"] for row in group]),
            "evaluation_range": [min(row["test_cross"] for row in group), max(row["test_cross"] for row in group)],
            "training_median": median([row["train_cross"] for row in group]),
            "coverage_ratio_min": min(row["coverage_ratio"] for row in group),
            "total_seconds_median": median([row["total_seconds"] for row in group]),
            "timing_valid": sum(row["timing_valid"] for row in group),
        }
    for left, right in comparisons:
        deltas = []
        relatives = []
        training_deltas = []
        for seed in sorted(by_seed):
            a, b = by_seed[seed][left], by_seed[seed][right]
            delta = a["test_cross"] - b["test_cross"]
            relative = delta / b["test_cross"] if b["test_cross"] else None
            deltas.append(delta)
            if relative is not None:
                relatives.append(relative)
            training_deltas.append(a["train_cross"] - b["train_cross"])
            paired_rows.append({
                "source": source, "seed": seed, "comparison": f"{left}-{right}",
                "left_evaluation": a["test_cross"], "right_evaluation": b["test_cross"],
                "delta": delta, "relative": relative,
                "training_delta": training_deltas[-1],
                "timing_valid": a["timing_valid"] and b["timing_valid"],
            })
        source_summary["comparisons"][f"{left}-{right}"] = {
            "paired_delta_median": median(deltas),
            "paired_relative_median": median(relatives) if relatives else None,
            "paired_delta_range": [min(deltas), max(deltas)],
            "positive_equal_negative": [sum(x > 0 for x in deltas), sum(x == 0 for x in deltas), sum(x < 0 for x in deltas)],
            "paired_training_delta_median": median(training_deltas),
            "valid_pairs": sum(row["timing_valid"] for row in paired_rows[-20:]),
        }
    summary["sources"][source] = source_summary

(ROOT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
with (ROOT / "paired.csv").open("w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(file, fieldnames=list(paired_rows[0]))
    writer.writeheader()
    writer.writerows(paired_rows)
with (ROOT / "runs.csv").open("w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(file, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(summary["sources"], indent=2))
