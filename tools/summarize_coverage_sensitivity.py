"""Summarize sealed matched-budget coverage-threshold sensitivity results."""

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


repo = Path(__file__).resolve().parents[1]
out = repo / "results" / "coverage_sensitivity"
payload = json.loads((out / "scores.json").read_text(encoding="utf-8"))
rows = payload["rows"]
probes = payload["probes"]
assert len(rows) == 360 and len(probes) == 60
by_pair = defaultdict(dict)
for row in rows:
    key = (row["source"], row["seed"], row["tau"])
    assert row["timing_valid"]
    by_pair[key][row["method"]] = row

paired = []
for (source, seed, tau), methods in sorted(by_pair.items()):
    assert set(methods) == {"A", "B"}
    a, b = methods["A"], methods["B"]
    difference = b["test_cross"] - a["test_cross"]
    paired.append({
        "source": source, "seed": seed, "tau": tau,
        "a_evaluation": a["test_cross"], "b_evaluation": b["test_cross"],
        "evaluation_delta": difference,
        "evaluation_relative": difference / a["test_cross"] if a["test_cross"] else "",
        "training_delta": b["train_cross"] - a["train_cross"],
        "a_coverage_ratio": a["coverage_ratio"],
        "b_coverage_ratio": b["coverage_ratio"],
    })

with (out / "paired.csv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, paired[0].keys())
    writer.writeheader()
    writer.writerows(paired)

probe_by_source = defaultdict(list)
for row in probes:
    probe_by_source[row["source"]].append(row)

summary = {}
for source in ("YDIP", "YBIOGRID", "SpaceI"):
    summary[source] = {}
    for tau in (0.90, 0.95, 0.98):
        group = [r for r in paired if r["source"] == source and r["tau"] == tau]
        original = [r for r in rows if r["source"] == source and r["tau"] == tau]
        assert len(group) == 20 and len(original) == 40
        delta = [r["evaluation_delta"] for r in group]
        relative = [r["evaluation_relative"] for r in group if r["evaluation_relative"] != ""]
        probe = [r["feasible_replacement_fraction"][f"{tau:.2f}"]
                 for r in probe_by_source[source]]
        summary[source][f"{tau:.2f}"] = {
            "n": 20,
            "evaluation_delta_median": statistics.median(delta),
            "evaluation_delta_range": [min(delta), max(delta)],
            "evaluation_relative_median": statistics.median(relative),
            "b_wins_ties_a_wins": [sum(v > 0 for v in delta), sum(v == 0 for v in delta), sum(v < 0 for v in delta)],
            "training_delta_median": statistics.median(r["training_delta"] for r in group),
            "a_evaluation_median": statistics.median(r["test_cross"] for r in original if r["method"] == "A"),
            "b_evaluation_median": statistics.median(r["test_cross"] for r in original if r["method"] == "B"),
            "a_coverage_ratio_median": statistics.median(r["coverage_ratio"] for r in original if r["method"] == "A"),
            "b_coverage_ratio_median": statistics.median(r["coverage_ratio"] for r in original if r["method"] == "B"),
            "a_total_seconds_median": statistics.median(r["total_seconds"] for r in original if r["method"] == "A"),
            "b_total_seconds_median": statistics.median(r["total_seconds"] for r in original if r["method"] == "B"),
            "reference_panel_feasible_replacement_fraction_median": statistics.median(probe),
            "reference_panel_feasible_replacement_fraction_range": [min(probe), max(probe)],
        }

(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))
