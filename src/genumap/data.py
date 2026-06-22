from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import ssl
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import certifi
import numpy as np
import pandas as pd

from genumap.config import DataConfig


@dataclass(frozen=True)
class GenotypeDataset:
    dosages: np.ndarray
    sample_ids: np.ndarray
    populations: np.ndarray
    super_populations: np.ndarray
    positions: np.ndarray
    variant_ids: np.ndarray
    source: str

    @property
    def n_samples(self) -> int:
        return self.dosages.shape[0]

    @property
    def n_variants(self) -> int:
        return self.dosages.shape[1]


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, destination: Path, expected_sha256: str | None = None) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if expected_sha256 is None or sha256_file(destination) == expected_sha256:
            return destination
        destination.unlink()

    partial = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "GENUMAP/0.1"})
    context = ssl.create_default_context(cafile=certifi.where())
    with urllib.request.urlopen(request, context=context) as response, partial.open("wb") as output:
        shutil.copyfileobj(response, output)
    if expected_sha256 is not None and sha256_file(partial) != expected_sha256:
        partial.unlink(missing_ok=True)
        raise ValueError(f"checksum mismatch for {url}")
    partial.replace(destination)
    return destination


def _processed_cache_path(config: DataConfig) -> Path:
    fields = {
        "max_variants": config.max_variants,
        "min_maf": config.min_maf,
        "max_missing_rate": config.max_missing_rate,
        "min_spacing_bp": config.min_spacing_bp,
        "vcf_url": config.vcf_url,
    }
    key = hashlib.sha256(json.dumps(fields, sort_keys=True).encode()).hexdigest()[:12]
    return config.cache_dir / f"processed_chr22_{key}.npz"


def _parse_gt(value: str) -> float:
    if value in {".", "./.", ".|."}:
        return np.nan
    alleles = value.replace("|", "/").split("/")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        return np.nan
    return float(int(alleles[0]) + int(alleles[1]))


def parse_phase3_vcf(
    vcf_path: Path,
    panel_path: Path,
    config: DataConfig,
    seed: int,
) -> GenotypeDataset:
    panel = pd.read_csv(panel_path, sep="\t", usecols=["sample", "pop", "super_pop"])
    panel = panel.set_index("sample")
    rng = np.random.default_rng(seed)

    selected: list[tuple[int, str, np.ndarray]] = []
    eligible_count = 0
    last_position = -config.min_spacing_bp
    sample_ids: list[str] | None = None

    with gzip.open(vcf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_ids = line.rstrip().split("\t")[9:]
                missing_panel = sorted(set(sample_ids) - set(panel.index))
                if missing_panel:
                    raise ValueError(f"panel is missing {len(missing_panel)} VCF samples")
                continue
            if line.startswith("#"):
                continue
            if sample_ids is None:
                raise ValueError("VCF sample header was not found")

            fields = line.rstrip().split("\t")
            position = int(fields[1])
            ref, alt, filter_value = fields[3], fields[4], fields[6]
            if (
                len(ref) != 1
                or len(alt) != 1
                or "," in alt
                or filter_value not in {"PASS", "."}
                or position - last_position < config.min_spacing_bp
            ):
                continue

            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                continue
            gt_index = format_fields.index("GT")
            dosage = np.fromiter(
                (_parse_gt(sample.split(":")[gt_index]) for sample in fields[9:]),
                dtype=np.float32,
                count=len(sample_ids),
            )
            missing_rate = float(np.isnan(dosage).mean())
            allele_frequency = float(np.nanmean(dosage) / 2)
            maf = min(allele_frequency, 1 - allele_frequency)
            if missing_rate > config.max_missing_rate or maf < config.min_maf:
                continue

            last_position = position
            eligible_count += 1
            record = (position, fields[2], dosage)
            if len(selected) < config.max_variants:
                selected.append(record)
            else:
                replacement = int(rng.integers(0, eligible_count))
                if replacement < config.max_variants:
                    selected[replacement] = record

    if sample_ids is None or not selected:
        raise ValueError("no variants passed the configured filters")

    selected.sort(key=lambda item: item[0])
    dosages = np.stack([item[2] for item in selected], axis=1)
    metadata = panel.loc[sample_ids]
    return GenotypeDataset(
        dosages=dosages,
        sample_ids=np.asarray(sample_ids),
        populations=metadata["pop"].to_numpy(str),
        super_populations=metadata["super_pop"].to_numpy(str),
        positions=np.asarray([item[0] for item in selected], dtype=np.int64),
        variant_ids=np.asarray([item[1] for item in selected]),
        source="1000 Genomes Project Phase 3 chromosome 22",
    )


def _save_dataset(path: Path, dataset: GenotypeDataset) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        dosages=dataset.dosages,
        sample_ids=dataset.sample_ids,
        populations=dataset.populations,
        super_populations=dataset.super_populations,
        positions=dataset.positions,
        variant_ids=dataset.variant_ids,
        source=np.asarray(dataset.source),
    )


def _load_dataset(path: Path) -> GenotypeDataset:
    with np.load(path, allow_pickle=False) as values:
        return GenotypeDataset(
            dosages=values["dosages"],
            sample_ids=values["sample_ids"],
            populations=values["populations"],
            super_populations=values["super_populations"],
            positions=values["positions"],
            variant_ids=values["variant_ids"],
            source=str(values["source"]),
        )


def load_empirical_dataset(config: DataConfig, seed: int) -> tuple[GenotypeDataset, dict[str, str]]:
    vcf_name = config.vcf_url.rsplit("/", maxsplit=1)[-1]
    panel_name = config.panel_url.rsplit("/", maxsplit=1)[-1]
    vcf_path = download_file(config.vcf_url, config.cache_dir / vcf_name, config.vcf_sha256)
    panel_path = download_file(
        config.panel_url,
        config.cache_dir / panel_name,
        config.panel_sha256,
    )
    hashes = {"vcf_sha256": sha256_file(vcf_path), "panel_sha256": sha256_file(panel_path)}

    processed_path = _processed_cache_path(config)
    if processed_path.exists():
        return _load_dataset(processed_path), hashes
    dataset = parse_phase3_vcf(vcf_path, panel_path, config, seed)
    _save_dataset(processed_path, dataset)
    return dataset, hashes


def make_synthetic_dataset(seed: int, n_variants: int = 400) -> GenotypeDataset:
    rng = np.random.default_rng(seed)
    super_names = np.asarray(["AFR", "AMR", "EAS", "EUR", "SAS"])
    pop_names = np.asarray(
        [f"{super_name}{index}" for super_name in super_names for index in (1, 2)]
    )
    samples_per_population = 20
    populations = np.repeat(pop_names, samples_per_population)
    super_populations = np.repeat(np.repeat(super_names, 2), samples_per_population)

    base = rng.beta(0.8, 0.8, size=n_variants)
    super_shift = rng.normal(0, 0.12, size=(len(super_names), n_variants))
    pop_shift = rng.normal(0, 0.035, size=(len(pop_names), n_variants))
    dosage_rows = []
    for pop_index in range(len(pop_names)):
        super_index = pop_index // 2
        probability = np.clip(base + super_shift[super_index] + pop_shift[pop_index], 0.02, 0.98)
        dosage_rows.append(
            rng.binomial(2, probability, size=(samples_per_population, n_variants)).astype(
                np.float32
            )
        )
    dosages = np.vstack(dosage_rows)
    return GenotypeDataset(
        dosages=dosages,
        sample_ids=np.asarray([f"SIM{index:04d}" for index in range(len(populations))]),
        populations=populations,
        super_populations=super_populations,
        positions=np.arange(1, n_variants + 1, dtype=np.int64) * 10_000,
        variant_ids=np.asarray([f"sim{index}" for index in range(n_variants)]),
        source="synthetic structured genotype fixture",
    )
