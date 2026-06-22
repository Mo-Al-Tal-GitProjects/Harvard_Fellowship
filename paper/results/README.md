# Reviewed empirical result snapshot

This directory contains aggregate outputs from the complete prespecified run:

- 2,504 1000 Genomes Phase 3 samples;
- 2,500 filtered and deterministically sampled chromosome 22 SNPs;
- 108 UMAP fits across three PCA input dimensions, four neighborhood sizes,
  three `min_dist` values, and three random seeds; and
- PCA and t-SNE reference views.

`manifest.json` records the configuration hash, input SHA-256 fingerprints,
software versions, platform, and source commit. `metrics.csv` contains one row per
embedding and no sample-level data. The corresponding reviewed figures are in
`paper/figures/`.

The first full run was generated from the modernization working tree before its
first GENUMAP commit, so the manifest correctly identifies the historical base
commit and marks the worktree as dirty. The reviewed code and configuration in
this repository are the source snapshot for the result.

The original ignored artifact directory was
`artifacts/20260622T025810Z-empirical-4dba1289bbad/`.
