from __future__ import annotations

import hashlib
import json
from pathlib import Path

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator


class DataConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    vcf_url: str
    panel_url: str
    vcf_sha256: str | None = None
    panel_sha256: str | None = None
    cache_dir: Path = Path("data/cache")
    max_variants: int = Field(ge=50)
    min_maf: float = Field(gt=0, lt=0.5)
    max_missing_rate: float = Field(ge=0, lt=1)
    min_spacing_bp: int = Field(ge=0)


class AnalysisConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    pca_components: list[int]
    umap_neighbors: list[int]
    umap_min_dist: list[float]
    seeds: list[int]
    metric_neighbors: int = Field(ge=2)
    distance_pairs: int = Field(ge=100)
    tsne_perplexity: float = Field(gt=1)

    @model_validator(mode="after")
    def validate_grid(self) -> AnalysisConfig:
        if not self.pca_components or not self.umap_neighbors or not self.umap_min_dist:
            raise ValueError("analysis parameter lists cannot be empty")
        if not self.seeds:
            raise ValueError("at least one random seed is required")
        if any(value <= 1 for value in self.umap_neighbors):
            raise ValueError("UMAP neighbor counts must be greater than one")
        if any(value < 0 for value in self.umap_min_dist):
            raise ValueError("UMAP min_dist values cannot be negative")
        return self


class StudyConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    study_name: str
    random_seed: int
    data: DataConfig
    analysis: AnalysisConfig
    output_dir: Path = Path("artifacts")

    def stable_hash(self, length: int = 12) -> str:
        payload = json.dumps(self.model_dump(mode="json"), sort_keys=True).encode()
        return hashlib.sha256(payload).hexdigest()[:length]

    def quick(self) -> StudyConfig:
        return self.model_copy(
            update={
                "data": self.data.model_copy(update={"max_variants": 400}),
                "analysis": self.analysis.model_copy(
                    update={
                        "pca_components": [10],
                        "umap_neighbors": [15, 50],
                        "umap_min_dist": [0.1],
                        "seeds": self.analysis.seeds[:2],
                        "distance_pairs": 2_000,
                    }
                ),
            }
        )


def load_config(path: Path) -> StudyConfig:
    with path.open(encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    return StudyConfig.model_validate(raw)
