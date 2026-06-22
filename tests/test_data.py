from __future__ import annotations

import gzip
from pathlib import Path

from genumap.config import DataConfig
from genumap.data import parse_phase3_vcf


def test_parse_tiny_vcf(tmp_path: Path) -> None:
    panel = tmp_path / "panel.tsv"
    panel.write_text(
        "sample\tpop\tsuper_pop\tgender\n"
        "S1\tP1\tG1\tfemale\n"
        "S2\tP1\tG1\tmale\n"
        "S3\tP2\tG2\tfemale\n"
        "S4\tP2\tG2\tmale\n",
        encoding="utf-8",
    )
    vcf = tmp_path / "tiny.vcf.gz"
    with gzip.open(vcf, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
        )
        handle.write("22\t10000\trs1\tA\tG\t100\tPASS\t.\tGT\t0|0\t0|1\t1|1\t0|1\n")
        handle.write("22\t20000\trs2\tC\tT\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\n")

    config = DataConfig(
        vcf_url="https://example.test/tiny.vcf.gz",
        panel_url="https://example.test/panel.tsv",
        cache_dir=tmp_path,
        max_variants=50,
        min_maf=0.05,
        max_missing_rate=0.1,
        min_spacing_bp=1,
    )
    dataset = parse_phase3_vcf(vcf, panel, config, seed=7)
    assert dataset.dosages.shape == (4, 2)
    assert dataset.populations.tolist() == ["P1", "P1", "P2", "P2"]
