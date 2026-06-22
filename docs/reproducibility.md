# Reproducibility

## Exact environment

```bash
uv sync --frozen
uv run genumap check --config configs/study.yaml
uv run genumap run --config configs/study.yaml
```

`uv.lock` fixes the dependency graph. Each run stores the resolved configuration,
configuration hash, Git commit, input hashes, package versions, platform, metrics,
embedding manifest, and generated figures under a timestamped artifact directory.

## Quick verification

```bash
uv run genumap run --config configs/study.yaml --synthetic --quick
```

This exercises the complete workflow with deterministic simulated genotypes. It
does not validate empirical conclusions.

## Expected resources

- Initial download: approximately 206 MB.
- Full run: 2,504 samples, up to 2,500 variants, and 108 UMAP fits.
- The first UMAP invocation may spend extra time compiling Numba kernels.

Artifact directories are intentionally ignored by Git. A reviewed release should
publish selected aggregate metrics and figures with their manifest, not cached
genotypes or participant-level embeddings.
