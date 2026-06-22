from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.spatial import procrustes
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, trustworthiness
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from umap import UMAP

from genumap.config import AnalysisConfig
from genumap.data import GenotypeDataset


@dataclass(frozen=True)
class StudyResults:
    standardized_genotypes: np.ndarray
    pca_scores: np.ndarray
    explained_variance_ratio: np.ndarray
    embeddings: dict[str, np.ndarray]
    embedding_metadata: pd.DataFrame
    metrics: pd.DataFrame


def standardize_dosages(dosages: np.ndarray) -> np.ndarray:
    values = dosages.astype(np.float32, copy=True)
    column_means = np.nanmean(values, axis=0)
    missing_rows, missing_columns = np.where(np.isnan(values))
    values[missing_rows, missing_columns] = column_means[missing_columns]
    allele_frequency = values.mean(axis=0) / 2
    scale = np.sqrt(2 * allele_frequency * (1 - allele_frequency))
    if np.any(scale <= 0):
        raise ValueError("monomorphic variants remain after filtering")
    return ((values - (2 * allele_frequency)) / scale).astype(np.float32)


def _knn_overlap(reference: np.ndarray, embedding: np.ndarray, neighbors: int) -> float:
    count = min(neighbors, len(reference) - 1)
    reference_neighbors = NearestNeighbors(n_neighbors=count + 1).fit(reference)
    embedding_neighbors = NearestNeighbors(n_neighbors=count + 1).fit(embedding)
    reference_indices = reference_neighbors.kneighbors(reference, return_distance=False)[:, 1:]
    embedding_indices = embedding_neighbors.kneighbors(embedding, return_distance=False)[:, 1:]
    overlap = [
        len(set(reference_row).intersection(embedding_row)) / count
        for reference_row, embedding_row in zip(reference_indices, embedding_indices, strict=True)
    ]
    return float(np.mean(overlap))


def _distance_spearman(
    reference: np.ndarray,
    embedding: np.ndarray,
    pairs: int,
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    first = rng.integers(0, len(reference), size=pairs * 2)
    second = rng.integers(0, len(reference), size=pairs * 2)
    keep = first != second
    first, second = first[keep][:pairs], second[keep][:pairs]
    reference_distance = np.linalg.norm(reference[first] - reference[second], axis=1)
    embedding_distance = np.linalg.norm(embedding[first] - embedding[second], axis=1)
    return float(spearmanr(reference_distance, embedding_distance).statistic)


def _group_outliers(values: np.ndarray, labels: np.ndarray, fraction: float = 0.02) -> set[int]:
    selected: set[int] = set()
    for label in np.unique(labels):
        indices = np.flatnonzero(labels == label)
        centroid = values[indices].mean(axis=0)
        distance = np.linalg.norm(values[indices] - centroid, axis=1)
        count = max(1, int(np.ceil(len(indices) * fraction)))
        selected.update(indices[np.argsort(distance)[-count:]].tolist())
    return selected


def _outlier_jaccard(reference: np.ndarray, embedding: np.ndarray, labels: np.ndarray) -> float:
    reference_outliers = _group_outliers(reference, labels)
    embedding_outliers = _group_outliers(embedding, labels)
    union = reference_outliers | embedding_outliers
    return len(reference_outliers & embedding_outliers) / len(union) if union else 1.0


def _embedding_metrics(
    reference: np.ndarray,
    embedding: np.ndarray,
    dataset: GenotypeDataset,
    analysis: AnalysisConfig,
    seed: int,
) -> dict[str, float]:
    neighbors = min(analysis.metric_neighbors, (len(reference) - 1) // 2)
    return {
        "trustworthiness": float(trustworthiness(reference, embedding, n_neighbors=neighbors)),
        "knn_overlap": _knn_overlap(reference, embedding, neighbors),
        "distance_spearman": _distance_spearman(
            reference,
            embedding,
            analysis.distance_pairs,
            seed,
        ),
        "population_silhouette": float(silhouette_score(embedding, dataset.populations)),
        "super_population_silhouette": float(
            silhouette_score(embedding, dataset.super_populations)
        ),
        "outlier_jaccard": _outlier_jaccard(reference, embedding, dataset.populations),
    }


def _embedding_key(
    method: str,
    pca_components: int,
    seed: int,
    neighbors: int | None = None,
    min_dist: float | None = None,
) -> str:
    key = f"{method}_d{pca_components}_s{seed}"
    if neighbors is not None:
        key += f"_n{neighbors}"
    if min_dist is not None:
        key += f"_m{min_dist:g}"
    return key


def run_analysis(
    dataset: GenotypeDataset,
    analysis: AnalysisConfig,
    random_seed: int,
) -> StudyResults:
    standardized = standardize_dosages(dataset.dosages)
    maximum_components = min(
        max(analysis.pca_components),
        standardized.shape[0] - 1,
        standardized.shape[1],
    )
    if maximum_components < 2:
        raise ValueError("at least two PCA components are required")
    pca = PCA(n_components=maximum_components, svd_solver="randomized", random_state=random_seed)
    pca_scores = pca.fit_transform(standardized).astype(np.float32)

    embeddings: dict[str, np.ndarray] = {}
    records: list[dict[str, float | int | str | None]] = []

    pca_key = _embedding_key("pca", maximum_components, random_seed)
    embeddings[pca_key] = pca_scores[:, :2]
    pca_metrics = _embedding_metrics(
        pca_scores,
        embeddings[pca_key],
        dataset,
        analysis,
        random_seed,
    )
    records.append(
        {
            "key": pca_key,
            "method": "PCA",
            "pca_components": maximum_components,
            "neighbors": None,
            "min_dist": None,
            "seed": random_seed,
            **pca_metrics,
        }
    )

    tsne_key = _embedding_key("tsne", maximum_components, random_seed)
    perplexity = min(analysis.tsne_perplexity, (len(dataset.dosages) - 1) / 3)
    embeddings[tsne_key] = TSNE(
        n_components=2,
        perplexity=perplexity,
        init="pca",
        learning_rate="auto",
        max_iter=1_000,
        random_state=random_seed,
    ).fit_transform(pca_scores)
    tsne_metrics = _embedding_metrics(
        pca_scores,
        embeddings[tsne_key],
        dataset,
        analysis,
        random_seed,
    )
    records.append(
        {
            "key": tsne_key,
            "method": "t-SNE",
            "pca_components": maximum_components,
            "neighbors": None,
            "min_dist": None,
            "seed": random_seed,
            **tsne_metrics,
        }
    )

    for components, neighbors, min_dist, seed in itertools.product(
        analysis.pca_components,
        analysis.umap_neighbors,
        analysis.umap_min_dist,
        analysis.seeds,
    ):
        effective_components = min(components, maximum_components)
        reference = pca_scores[:, :effective_components]
        key = _embedding_key("umap", effective_components, seed, neighbors, min_dist)
        embeddings[key] = UMAP(
            n_components=2,
            n_neighbors=min(neighbors, len(reference) - 1),
            min_dist=min_dist,
            metric="euclidean",
            init="spectral",
            random_state=seed,
            n_jobs=1,
        ).fit_transform(reference)
        metrics = _embedding_metrics(reference, embeddings[key], dataset, analysis, seed)
        records.append(
            {
                "key": key,
                "method": "UMAP",
                "pca_components": effective_components,
                "neighbors": neighbors,
                "min_dist": min_dist,
                "seed": seed,
                **metrics,
            }
        )

    metrics_frame = pd.DataFrame.from_records(records)
    metrics_frame["seed_disparity"] = np.nan
    umap_rows = metrics_frame[metrics_frame["method"] == "UMAP"]
    for _, group in umap_rows.groupby(["pca_components", "neighbors", "min_dist"]):
        keys = group.sort_values("seed")["key"].tolist()
        if len(keys) < 2:
            disparity = 0.0
        else:
            reference = embeddings[keys[0]]
            disparities = [procrustes(reference, embeddings[key])[2] for key in keys[1:]]
            disparity = float(np.mean(disparities))
        metrics_frame.loc[metrics_frame["key"].isin(keys), "seed_disparity"] = disparity

    metadata = metrics_frame[
        ["key", "method", "pca_components", "neighbors", "min_dist", "seed"]
    ].copy()
    return StudyResults(
        standardized_genotypes=standardized,
        pca_scores=pca_scores,
        explained_variance_ratio=pca.explained_variance_ratio_,
        embeddings=embeddings,
        embedding_metadata=metadata,
        metrics=metrics_frame,
    )
