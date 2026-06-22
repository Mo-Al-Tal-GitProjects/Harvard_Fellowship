from __future__ import annotations

import json
import platform
import subprocess
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

from genumap.analysis import StudyResults
from genumap.config import StudyConfig
from genumap.data import GenotypeDataset

plt.switch_backend("Agg")


SUPER_POPULATION_COLORS = {
    "AFR": "#7b3294",
    "AMR": "#e66101",
    "EAS": "#018571",
    "EUR": "#0571b0",
    "SAS": "#ca0020",
}


def _scatter(ax: plt.Axes, embedding: np.ndarray, labels: np.ndarray, title: str) -> None:
    for label in np.unique(labels):
        selected = labels == label
        ax.scatter(
            embedding[selected, 0],
            embedding[selected, 1],
            s=7,
            alpha=0.72,
            linewidths=0,
            label=label,
            color=SUPER_POPULATION_COLORS.get(label),
        )
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])


def _select_baseline_keys(results: StudyResults) -> tuple[str, str, str]:
    metadata = results.embedding_metadata
    pca_key = metadata.loc[metadata["method"] == "PCA", "key"].iloc[0]
    tsne_key = metadata.loc[metadata["method"] == "t-SNE", "key"].iloc[0]
    umap = metadata[metadata["method"] == "UMAP"].copy()
    maximum_dimension = int(umap["pca_components"].max())
    umap = umap[umap["pca_components"] == maximum_dimension]
    neighbor_target = 15 if 15 in umap["neighbors"].values else int(umap["neighbors"].iloc[0])
    distance_target = 0.1 if 0.1 in umap["min_dist"].values else float(umap["min_dist"].iloc[0])
    umap = umap[(umap["neighbors"] == neighbor_target) & (umap["min_dist"] == distance_target)]
    umap_key = umap.sort_values("seed")["key"].iloc[0]
    return pca_key, umap_key, tsne_key


def make_figures(run_dir: Path, dataset: GenotypeDataset, results: StudyResults) -> None:
    sns.set_theme(style="white", context="paper")
    figure_dir = run_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    pca_key, umap_key, tsne_key = _select_baseline_keys(results)
    figure, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    for axis, key, title in zip(
        axes,
        (pca_key, umap_key, tsne_key),
        ("PCA baseline", "UMAP baseline", "t-SNE baseline"),
        strict=True,
    ):
        _scatter(axis, results.embeddings[key], dataset.super_populations, title)
    handles, labels = axes[-1].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=5, frameon=False)
    figure.suptitle("Method comparison colored by 1000 Genomes super-population")
    figure.tight_layout(rect=(0, 0.1, 1, 0.95))
    figure.savefig(figure_dir / "method_comparison.png", dpi=220, bbox_inches="tight")
    plt.close(figure)

    umap_metadata = results.embedding_metadata[
        results.embedding_metadata["method"] == "UMAP"
    ].copy()
    maximum_dimension = int(umap_metadata["pca_components"].max())
    first_seed = int(umap_metadata["seed"].min())
    grid = umap_metadata[
        (umap_metadata["pca_components"] == maximum_dimension)
        & (umap_metadata["seed"] == first_seed)
    ]
    neighbors = sorted(grid["neighbors"].astype(int).unique())
    min_distances = sorted(grid["min_dist"].astype(float).unique())
    figure, axes = plt.subplots(
        len(neighbors),
        len(min_distances),
        figsize=(3.2 * len(min_distances), 2.7 * len(neighbors)),
        squeeze=False,
    )
    for row, neighbor_count in enumerate(neighbors):
        for column, min_dist in enumerate(min_distances):
            record = grid[(grid["neighbors"] == neighbor_count) & (grid["min_dist"] == min_dist)]
            key = record["key"].iloc[0]
            _scatter(
                axes[row, column],
                results.embeddings[key],
                dataset.super_populations,
                f"neighbors={neighbor_count}, min_dist={min_dist:g}",
            )
    handles, labels = axes[-1, -1].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=5, frameon=False)
    figure.suptitle(f"UMAP parameter sensitivity ({maximum_dimension} PCs, seed {first_seed})")
    figure.tight_layout(rect=(0, 0.05, 1, 0.97))
    figure.savefig(figure_dir / "umap_sensitivity.png", dpi=220, bbox_inches="tight")
    plt.close(figure)

    umap_metrics = results.metrics[results.metrics["method"] == "UMAP"].copy()
    summary = (
        umap_metrics.groupby(["pca_components", "neighbors", "min_dist"], as_index=False)
        .agg(
            trustworthiness=("trustworthiness", "mean"),
            knn_overlap=("knn_overlap", "mean"),
            seed_disparity=("seed_disparity", "mean"),
        )
        .sort_values(["pca_components", "neighbors", "min_dist"])
    )
    figure, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    for axis, metric, title in zip(
        axes,
        ("trustworthiness", "knn_overlap", "seed_disparity"),
        ("Trustworthiness", "k-NN overlap", "Seed disparity (lower is better)"),
        strict=True,
    ):
        sns.lineplot(
            data=summary,
            x="neighbors",
            y=metric,
            hue="min_dist",
            style="pca_components",
            markers=True,
            palette="viridis",
            ax=axis,
        )
        axis.set_title(title)
        axis.set_xscale("log")
        axis.get_legend().remove()
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=6, frameon=False)
    figure.suptitle("Quantitative sensitivity diagnostics")
    figure.tight_layout(rect=(0, 0.15, 1, 0.95))
    figure.savefig(figure_dir / "sensitivity_metrics.png", dpi=220, bbox_inches="tight")
    plt.close(figure)


def _git_commit() -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _git_worktree_dirty() -> bool | None:
    result = subprocess.run(
        ["git", "status", "--porcelain"],
        check=False,
        capture_output=True,
        text=True,
    )
    return bool(result.stdout.strip()) if result.returncode == 0 else None


def write_artifacts(
    run_dir: Path,
    config: StudyConfig,
    dataset: GenotypeDataset,
    results: StudyResults,
    input_hashes: dict[str, str],
    synthetic: bool,
    quick: bool,
) -> None:
    run_dir.mkdir(parents=True, exist_ok=False)
    with (run_dir / "resolved_config.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config.model_dump(mode="json"), handle, sort_keys=False)

    results.metrics.to_csv(run_dir / "metrics.csv", index=False)
    results.embedding_metadata.to_csv(run_dir / "embedding_manifest.csv", index=False)
    np.savez_compressed(run_dir / "embeddings.npz", **results.embeddings)
    pd.DataFrame(
        {
            "population": dataset.populations,
            "super_population": dataset.super_populations,
        }
    ).value_counts().rename("samples").reset_index().to_csv(
        run_dir / "sample_counts.csv",
        index=False,
    )
    pd.DataFrame(
        {
            "component": np.arange(1, len(results.explained_variance_ratio) + 1),
            "explained_variance_ratio": results.explained_variance_ratio,
        }
    ).to_csv(run_dir / "pca_variance.csv", index=False)

    umap_metrics = results.metrics[results.metrics["method"] == "UMAP"]
    pca_baseline = results.metrics[results.metrics["method"] == "PCA"].iloc[0]
    tsne_baseline = results.metrics[results.metrics["method"] == "t-SNE"].iloc[0]
    best_fidelity = umap_metrics.sort_values(
        ["trustworthiness", "knn_overlap"], ascending=False
    ).iloc[0]
    best_stability = umap_metrics.sort_values("seed_disparity", ascending=True).iloc[0]
    status = (
        "Synthetic verification run: these numbers are pipeline diagnostics, not study results."
        if synthetic
        else "Empirical run on the public 1000 Genomes Phase 3 chromosome 22 call set."
    )
    summary = f"""# GENUMAP run summary

{status}

## Data

- Source: {dataset.source}
- Samples: {dataset.n_samples:,}
- Variants after filtering/subsampling: {dataset.n_variants:,}
- Populations: {len(np.unique(dataset.populations))}
- Super-populations: {len(np.unique(dataset.super_populations))}

## Prespecified diagnostics

- Highest observed UMAP trustworthiness: {best_fidelity['trustworthiness']:.4f}
  (`{best_fidelity['key']}`)
- Associated k-NN overlap: {best_fidelity['knn_overlap']:.4f}
- Lowest mean seed disparity: {best_stability['seed_disparity']:.6f}
  (`{best_stability['key']}`)

Across all UMAP fits:

- Trustworthiness range: {umap_metrics['trustworthiness'].min():.4f} to
  {umap_metrics['trustworthiness'].max():.4f}
- k-NN overlap range: {umap_metrics['knn_overlap'].min():.4f} to
  {umap_metrics['knn_overlap'].max():.4f}
- Pairwise-distance Spearman range: {umap_metrics['distance_spearman'].min():.4f}
  to {umap_metrics['distance_spearman'].max():.4f}
- Seed-disparity range: {umap_metrics['seed_disparity'].min():.6f} to
  {umap_metrics['seed_disparity'].max():.6f}
- Super-population silhouette range:
  {umap_metrics['super_population_silhouette'].min():.4f} to
  {umap_metrics['super_population_silhouette'].max():.4f}
- Population silhouette range: {umap_metrics['population_silhouette'].min():.4f}
  to {umap_metrics['population_silhouette'].max():.4f}

For the reference views, PCA preserved sampled global distances most strongly
(Spearman {pca_baseline['distance_spearman']:.4f}), while t-SNE had local
trustworthiness {tsne_baseline['trustworthiness']:.4f} and k-NN overlap
{tsne_baseline['knn_overlap']:.4f}.

These rankings are descriptive. A visually separated embedding is not evidence that
sampling groups are discrete biological categories. Interpret the full metric
distribution and parameter sensitivity figures, not a single preferred plot.
"""
    (run_dir / "results.md").write_text(summary, encoding="utf-8")

    manifest = {
        "created_at_utc": datetime.now(UTC).isoformat(),
        "study_name": config.study_name,
        "config_hash": config.stable_hash(),
        "git_commit": _git_commit(),
        "git_worktree_dirty": _git_worktree_dirty(),
        "synthetic": synthetic,
        "quick": quick,
        "input_hashes": input_hashes,
        "dataset": {
            "source": dataset.source,
            "samples": dataset.n_samples,
            "variants": dataset.n_variants,
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": version("numpy"),
            "pandas": version("pandas"),
            "scikit-learn": version("scikit-learn"),
            "umap-learn": version("umap-learn"),
        },
    }
    (run_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    make_figures(run_dir, dataset, results)
