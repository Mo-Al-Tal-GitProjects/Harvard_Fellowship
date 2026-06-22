# Study design

## Research question

How stable and faithful are two-dimensional UMAP representations of common human
genetic variation under reasonable changes to input dimension, neighborhood
size, minimum distance, and random seed?

## Data

The empirical analysis uses the public 1000 Genomes Project Phase 3 chromosome
22 call set based on the 2013-05-02 sequence freeze. The release contains 2,504
individuals from 26 sampling populations and five super-populations. GENUMAP
downloads the call set and panel directly from the International Genome Sample
Resource; upstream data are never committed to this repository.

Chromosome 22 is a deliberately bounded first study. It makes the experiment
practical on a laptop while retaining real autosomal population variation. It
does not support claims about whole-genome effect sizes.

## Variant selection and preprocessing

The prespecified pipeline:

1. retains `PASS` biallelic SNPs;
2. requires minor allele frequency at least 0.05;
3. requires genotype missingness no greater than 0.02;
4. considers variants at least 10 kb apart;
5. uses deterministic reservoir sampling to retain at most 2,500 eligible SNPs
   across the chromosome;
6. mean-imputes the remaining missing dosages;
7. centers each SNP at `2p` and scales it by `sqrt(2p(1-p))`; and
8. computes randomized PCA with a fixed seed.

Physical spacing is a transparent approximation to linkage-disequilibrium
pruning, not a substitute for a population-aware LD analysis. This limitation
must accompany reported results.

## Prespecified experiment

UMAP is applied to 10, 25, and 50 principal components. The grid crosses:

- `n_neighbors`: 5, 15, 50, 100;
- `min_dist`: 0.0, 0.1, 0.5; and
- seeds: 2024, 2025, 2026.

This produces 108 UMAP fits. PCA and t-SNE are included as reference
visualizations, not as competitors in a universal benchmark.

## Metrics

- Trustworthiness of the two-dimensional representation.
- Mean overlap between input-space and embedding-space k-nearest neighbors.
- Spearman correlation of sampled pairwise distances.
- Population and super-population silhouette scores, treated descriptively.
- Jaccard overlap of within-population outliers between input and display space.
- Procrustes disparity across random seeds for a fixed parameter setting.

No single metric selects a “correct” embedding. Separation scores do not validate
the labels as natural or discrete biological categories.

## Interpretation boundary

1000 Genomes population labels describe sampling locations and recruitment
histories; super-populations are broad analytical groupings. Neither is a proxy
for race, and neither should be used to infer genetic determinism. GENUMAP studies
the behavior of a visualization algorithm, not the legitimacy of social
categories.
