from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from genumap.analysis import run_analysis
from genumap.config import StudyConfig
from genumap.data import load_empirical_dataset, make_synthetic_dataset
from genumap.reporting import write_artifacts


def run_study(config: StudyConfig, *, synthetic: bool = False, quick: bool = False) -> Path:
    resolved = config.quick() if quick else config
    if synthetic:
        dataset = make_synthetic_dataset(resolved.random_seed, resolved.data.max_variants)
        hashes: dict[str, str] = {}
        mode = "synthetic"
    else:
        dataset, hashes = load_empirical_dataset(resolved.data, resolved.random_seed)
        mode = "empirical"

    results = run_analysis(dataset, resolved.analysis, resolved.random_seed)
    timestamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    run_dir = resolved.output_dir / f"{timestamp}-{mode}-{resolved.stable_hash()}"
    write_artifacts(run_dir, resolved, dataset, results, hashes, synthetic, quick)
    return run_dir
