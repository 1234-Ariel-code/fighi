from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .search import SearchResult


def _save(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_plots(outdir: Path, result: SearchResult) -> dict[str, str]:
    plot_dir = outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    if not result.interactions:
        return paths

    ranked = sorted(result.interactions, key=lambda item: item.p_value)
    q_values = np.clip(np.asarray([item.q_value for item in ranked]), 1e-300, 1.0)
    colors = ["#0f766e" if item.significant else "#94a3b8" for item in ranked]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.scatter(np.arange(1, len(ranked) + 1), -np.log10(q_values), c=colors, s=18, alpha=0.8)
    ax.axhline(-np.log10(result.config.alpha), color="#b91c1c", linestyle="--", linewidth=1.2)
    ax.set(
        xlabel="Interaction rank",
        ylabel=r"$-\log_{10}$(adjusted p-value)",
        title="Global interaction evidence",
    )
    ax.grid(axis="y", color="#e2e8f0", linewidth=0.8)
    path = plot_dir / "interaction_evidence.png"
    _save(fig, path)
    paths["interaction_evidence"] = str(path)

    orders = sorted({item.order for item in result.interactions})
    values = [
        [item.fi_gain for item in result.interactions if item.order == order] for order in orders
    ]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.boxplot(values, tick_labels=[str(order) for order in orders], showfliers=False)
    ax.set(
        xlabel="Interaction order",
        ylabel="Fisher-information gain",
        title="Evidence distribution by interaction order",
    )
    ax.grid(axis="y", color="#e2e8f0", linewidth=0.8)
    path = plot_dir / "fi_gain_by_order.png"
    _save(fig, path)
    paths["fi_gain_by_order"] = str(path)

    top = ranked[: min(30, len(ranked))][::-1]
    fig_height = max(4.8, 0.28 * len(top) + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_height))
    labels = [" × ".join(item.features) for item in top]
    scores = [-np.log10(max(item.q_value, 1e-300)) for item in top]
    bar_colors = ["#0f766e" if item.significant else "#94a3b8" for item in top]
    ax.barh(np.arange(len(top)), scores, color=bar_colors)
    ax.set_yticks(np.arange(len(top)), labels=labels)
    ax.axvline(-np.log10(result.config.alpha), color="#b91c1c", linestyle="--", linewidth=1.2)
    ax.set(xlabel=r"$-\log_{10}$(adjusted p-value)", title="Top interaction candidates")
    ax.grid(axis="x", color="#e2e8f0", linewidth=0.8)
    path = plot_dir / "top_interactions.png"
    _save(fig, path)
    paths["top_interactions"] = str(path)
    return paths
