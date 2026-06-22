# Data handling

GENUMAP uses public-reference human genotype data, which still warrants careful
handling.

## Rules

- Download source data only from the configured International Genome Sample
  Resource URLs.
- Do not commit VCFs, sample-level genotype matrices, processed caches, or
  individual sample identifiers.
- Keep `data/cache/` and `artifacts/` ignored by Git.
- Record SHA-256 hashes in each empirical run manifest.
- Publish aggregate metrics and reviewed figures only.
- Retain the 1000 Genomes citation and data-reuse statements in any public result.

The pipeline writes population counts but omits sample IDs from result tables.
Embeddings remain local artifacts because each row represents an individual
reference sample even when the source identifiers are omitted.

All of Us participant-level data are not used by GENUMAP and must not be added to
this repository.
