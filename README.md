# udon

**U**pstream **D**istance in the **O**mnigenic **N**etwork scoring — an R package for
ranking proteins by their predicted influence on a trait of interest.

Given a list of "core" genes presumably associated with a trait (e.g. GWAS hits, rare-variant
associations), `udon` scores all ~2,900 proteins in a precomputed
gene-regulation network by how strongly their activity propagates to those
genes. The result is a ranked table that prioritises proteins upstream of
known trait-relevant biology.

## Installation

```r
# install.packages("remotes")
remotes::install_github("YOUR_ORG/udon")
```

The package ships precomputed network matrices (~62 MB); no external data
download is required.

## Quick start

### Without core genes' effect sizes to the trait (UDON_approx)

Provide gene names only — each is assigned equal weight:

```r
library(udon)

results <- udon_score(c("PCSK9", "LDLR", "ANGPTL3", "APOB", "SLC4A1", "ASGR1"))
head(results, 10)
#>        gene     score abs_score rank is_core
#> 1  SERPINA1 -3.235713  3.235713    1        FALSE
#> 2   ANGPTL3  2.608621  2.608621    2         TRUE
#> 3     ASGR1  2.141887  2.141887    3         TRUE
#> 4    SLC4A1  2.001943  2.001943    4         TRUE
#> 5     PCSK9  1.999791  1.999791    5         TRUE
#> 6      LDLR  1.999661  1.999661    6         TRUE
#> 7     APOA1 -1.959263  1.959263    7        FALSE
#> 8      APOB  1.007791  1.007791    8         TRUE
#> 9      LCAT -0.681997  0.681997    9        FALSE
#> 10    PCSK7 -0.503774  0.503774   10        FALSE
```

### With effect sizes (UDON)

Provide per-gene effect sizes (e.g. rare-variant LoF betas) via `lof_betas`.
Genes in `candidates` that are missing from `lof_betas` fall back to unit
weight automatically, so a partial beta vector is fine:

```r
results <- udon_score(
  core_genes = c("PCSK9", "LDLR", "ANGPTL3", "APOB", "SLC4A1", "ASGR1"),
  lof_betas  = c(PCSK9   = -0.985,
                 LDLR    =  0.284,
                 ANGPTL3 = -0.314,
                 APOB    = -1.765,
                 SLC4A1  = -1.026,
                 ASGR1   = -0.278)
)
head(results, 10)
#>        gene     score abs_score rank is_core
#> 1    SLC4A1  2.052516  2.052516    1         TRUE
#> 2     PCSK9  1.970492  1.970492    2         TRUE
#> 3      APOB  1.776362  1.776362    3         TRUE
#> 4  SERPINA1 -0.659056  0.659056    4        FALSE
#> 5   ANGPTL3  0.625289  0.625289    5         TRUE
#> 6      LDLR -0.569209  0.569209    6         TRUE
#> 7     ASGR1  0.534990  0.534990    7         TRUE
#> ...
```

## Output columns

| Column | Description |
|---|---|
| `gene` | Gene symbol |
| `score` | Signed UDON score from the causal network. Positive = activity propagates toward the candidate set; negative = anti-correlated propagation. |
| `abs_score` | `|score|` from the causal network — primary ranking key |
| `rank` | Rank 1 = strongest predicted influence. Primary: `abs_score` descending; tie-breaker: dense network score descending (for genes with causal score = 0). |
| `is_core` | `TRUE` if the gene is in the input list |

## Other functions

```r
# List all gene symbols present in the network (~2,900 proteins)
network_genes()
```

## Method

UDON scores gene *i* as:

$$\text{UDON}(i) = \sum_{c \in \text{candidates}} T[i, c] \cdot w(c)$$

where $T = M(I - M)^{-1}$ is the network propagation matrix (with diagonal
set to 1 so direct effects always count), and the weight is:

$$w(c) = \begin{cases} -\beta_\text{LoF}[c] & \text{if a LoF beta is provided} \\ 1 & \text{otherwise (UDON\_approx)} \end{cases}$$

The score $T[i,c]$ is computed from two gene-regulation networks:

- **Causal network** — sparse MVMR + rare-variant LoF edges (2,644 genes).
  This is the primary network; it determines the score and rank for any
  gene reachable from the candidates.
- **Dense network** — genome-wide raw protein-QTL betas (2,896 genes).
  Used only as a tie-breaker: genes that receive zero score from the causal
  network (not connected to the candidates) are ordered by their dense-network
  score.

Both networks are spectrally damped (largest eigenvalue capped at 0.95)
before inversion to ensure stable propagation.

## Notes

- Input genes not found in either network are skipped with a warning. You can
  check coverage with `intersect(your_genes, network_genes())`.
- The sign of `score` reflects the net direction of propagation. When
  candidate genes have mixed effect directions (e.g. some LoF betas positive,
  some negative), ranking by `abs_score` is appropriate.
- The network matrices are loaded once per R session on the first `udon_score()`
  call (~5–10 seconds). Subsequent calls with different gene lists are fast.

## Citation

If you use this package, please cite:

> [Paper title, authors, journal, year — add when published]

The precomputed network matrices are derived from UK Biobank protein and
genetic data. Details of the network construction are described in the paper.
