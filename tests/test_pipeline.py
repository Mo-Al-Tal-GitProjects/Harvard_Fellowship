from pathlib import Path

from genumap.config import AnalysisConfig, DataConfig, StudyConfig
from genumap.pipeline import run_study


def test_synthetic_pipeline_writes_complete_artifacts(tmp_path: Path) -> None:
    config = StudyConfig(
        study_name="test",
        random_seed=11,
        data=DataConfig(
            vcf_url="https://example.test/test.vcf.gz",
            panel_url="https://example.test/panel.tsv",
            cache_dir=tmp_path / "cache",
            max_variants=80,
            min_maf=0.05,
            max_missing_rate=0.02,
            min_spacing_bp=10_000,
        ),
        analysis=AnalysisConfig(
            pca_components=[5],
            umap_neighbors=[10],
            umap_min_dist=[0.1],
            seeds=[11, 12],
            metric_neighbors=5,
            distance_pairs=200,
            tsne_perplexity=10,
        ),
        output_dir=tmp_path / "artifacts",
    )
    run_dir = run_study(config, synthetic=True)
    assert (run_dir / "manifest.json").exists()
    assert (run_dir / "metrics.csv").exists()
    assert (run_dir / "results.md").exists()
    assert (run_dir / "figures" / "method_comparison.png").exists()
    assert (run_dir / "figures" / "umap_sensitivity.png").exists()
