"""Deterministic rendering of the two compact computational check figures."""
from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/riad_ridd_mpl")
import matplotlib.pyplot as plt


def render_check_figures(output_dir: Path) -> list[Path]:
    """Render the frozen model comparison and M2 probability check plots."""
    output_dir = Path(output_dir)
    with (output_dir / "metrics_summary.csv").open(encoding="utf-8", newline="") as handle:
        metrics = list(csv.DictReader(handle))
    core = [r for r in metrics if r["model"].startswith(("M0_", "M1_", "M2_", "M3_SYSTEM"))]
    model_plot = output_dir / "model_comparison.png"
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar([r["model"].split("_")[0] for r in core], [float(r["pair_accuracy"]) for r in core])
    ax.set(ylim=(0, 1), ylabel="LOSO pair accuracy", title="Frozen M0-M3 comparison")
    fig.tight_layout()
    fig.savefig(model_plot, dpi=180)
    plt.close(fig)

    with (output_dir / "pairwise_probabilities.csv").open(encoding="utf-8", newline="") as handle:
        pred = [r for r in csv.DictReader(handle) if r["model"] == "M2_PRIOR_PLUS_PSCI"]
    probability_plot = output_dir / "m2_pair_probabilities.png"
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.bar(range(24), [float(r["probability_i_gt_j"]) for r in pred])
    ax.axhline(0.5, color="black", lw=0.8)
    ax.set(xlabel="Frozen unordered pair", ylabel="P(i>j)", title="M2 out-of-fold probabilities")
    fig.tight_layout()
    fig.savefig(probability_plot, dpi=180)
    plt.close(fig)
    return [model_plot, probability_plot]
