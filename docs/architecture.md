# Architecture

GENUMAP is an experiment repository, not a general-purpose software platform.

```text
config YAML
    -> download/cache one fixed public dataset
    -> validate and select variants
    -> standardize dosages and compute PCA
    -> run the prespecified PCA/UMAP/t-SNE experiment
    -> calculate metrics
    -> write one versioned artifact directory
```

The package has four substantive modules:

- `data.py`: download, validation, VCF parsing, and synthetic test data;
- `analysis.py`: preprocessing, embeddings, stability, and fidelity metrics;
- `reporting.py`: figures, result summary, and provenance manifest; and
- `pipeline.py`: the single end-to-end workflow.

The CLI intentionally exposes only `check` and `run`. Data preparation, figure
generation, and reporting are implementation stages of the study and are not
presented as reusable user-facing tools.
