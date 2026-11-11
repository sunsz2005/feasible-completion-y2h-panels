"""Plot paired training and evaluation differences in the SA comparison."""

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


root = Path(__file__).resolve().parents[1]
result_dir = root / "results" / "independent_sa"
with (result_dir / "paired.csv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle))

names = [("YDIP", "YDIP"), ("YBIOGRID", "YBIOGRID"), ("SpaceI", "Human Space-I")]
fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.4), constrained_layout=True)
report = {}
for ax, (source, label) in zip(axes, names):
    selected = [r for r in rows if r["source"] == source and r["comparison"] == "B-SA"]
    train_delta = [int(r["training_delta"]) for r in selected]
    eval_delta = [int(r["delta"]) for r in selected]
    rho, _ = spearmanr(train_delta, eval_delta)
    report[source] = {
        "n": len(selected),
        "spearman_rho": float(rho),
        "b_eval_wins": sum(v > 0 for v in eval_delta),
        "eval_ties": sum(v == 0 for v in eval_delta),
        "sa_eval_wins": sum(v < 0 for v in eval_delta),
        "training_evaluation_sign_discordances": sum(
            a * b < 0 for a, b in zip(train_delta, eval_delta)
        ),
        "training_ties": sum(v == 0 for v in train_delta),
    }
    ax.axhline(0, color="#8f8f8f", linewidth=0.8)
    ax.axvline(0, color="#8f8f8f", linewidth=0.8)
    ax.scatter(
        train_delta, eval_delta, s=34, color="#135f8c", alpha=0.78,
        edgecolor="white", linewidth=0.4,
    )
    ax.set_title(f"{label} (n=20, rho={rho:.2f})")
    ax.set_xlabel("B - SA training edges")
    ax.set_ylabel("B - SA evaluation edges")
    ax.grid(alpha=0.15)

fig.savefig(result_dir / "train_eval_deltas.pdf", dpi=300)
(result_dir / "train_eval_correlation.json").write_text(
    json.dumps(report, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(report, indent=2))
