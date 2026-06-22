from pathlib import Path

from genumap.config import load_config


def test_study_config_is_valid() -> None:
    config = load_config(Path("configs/study.yaml"))
    assert config.study_name == "1000g-chr22-umap-sensitivity"
    assert len(config.analysis.pca_components) == 3
    assert config.stable_hash() == config.stable_hash()


def test_quick_config_reduces_grid() -> None:
    config = load_config(Path("configs/study.yaml"))
    quick = config.quick()
    assert quick.data.max_variants == 400
    assert quick.analysis.pca_components == [10]
    assert len(quick.analysis.seeds) == 2
