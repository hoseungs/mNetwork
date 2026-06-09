# mNetwork

`mNetwork` provides R functions for **conditional, nonparametric,
inference-based** network construction, designed in particular for microbiome
co-occurrence networks.

This repository provides the R functions for the paper:

> **Constructing microbiome co-occurrence networks with confidence: A
> conditional, nonparametric, inference-based approach.**

> Hoseung Song, Yunhua Xiang et al.

The method is built on the **Scaled Expected Conditional Covariance (SEcov)**:

- **Assumption-free.** SEcov makes no parametric assumption about the data
  distribution or the form of the dependence between taxa. It captures
  non-linear conditional associations and reduces to partial correlation when
  the data are Gaussian.
- **A global, robust measure.** Unlike partial or conditional correlation,
  SEcov is a population-level quantity that does not vary with OTU abundance,
  yielding a more robust network. It does not require estimating the precision
  matrix, so it scales to high-dimensional data.
- **Formal inference.** The SEcov estimator is asymptotically normal, so the
  method produces a p-value for every taxon pair. Edges are selected by
  statistical significance rather than an arbitrary correlation cutoff or
  regularization parameter, and the false discovery rate can be controlled.
- **Off-the-shelf machine learning.** The required conditional means are
  estimated with random forests (via the [`ranger`](https://cran.r-project.org/package=ranger)
  package), but any flexible regression method can be substituted.

---

## Files in this repository

| File | Description |
|------|-------------|
| `mNetwork.R`        | The main function `mnetwork()`. |
| `make_ADJ.R`        | The helper function `make_adj()`, used by `mnetwork()`. |
| `Sample_Data.RData` | Example clr-transformed OTU table (102 samples × 22 OTUs). |

---

## Getting started

These are the core R functions. Download the files, load them with `source()`,
and make sure the `ranger` package is installed.

```r
# install ranger once (required: mnetwork uses random forests)
# install.packages("ranger")
library(ranger)

# load the functions (both files are needed)
source("mNetwork.R")
source("make_ADJ.R")
```

---

## Main function: `mnetwork()`

```r
mnetwork(dat, d, e)
```

| Argument | Description |
|----------|-------------|
| `dat` | An `n × d` data set used for network construction (e.g. an OTU table with `n` samples and `d` taxa). For microbiome counts we recommend normalizing first, e.g. with a centered log-ratio (clr) transformation — see below. |
| `d`   | The number of variables (columns) in `dat`. |
| `e`   | The number of edges to keep in the constructed network. We recommend `e = d`. |

It returns a list with three elements:

| Element | Description |
|---------|-------------|
| `p_mat`        | A `d × d` matrix of p-values for each taxon pair, obtained from the asymptotic distribution of SEcov. |
| `adj`          | The estimated adjacency matrix. Exactly `e` edges (the `e` most significant pairs) are connected. |
| `partial_corr` | The estimated SEcov partial-correlation matrix (edge weights; positive = positive association, negative = negative association). |

---

## Quick example

Using the included `Sample_Data` (a clr-transformed OTU table, 102 samples × 22 OTUs):

```r
library(ranger)
source("mNetwork.R")
source("make_ADJ.R")

load("Sample_Data.RData")
OTU <- Sample_Data           # 102 samples x 22 OTUs (already clr-transformed)

d <- ncol(OTU)               # number of OTUs (= 22)
e <- d                       # number of edges to keep (recommended: e = d)

set.seed(1)                  # random forests are stochastic — set a seed for reproducibility
fit <- mnetwork(OTU, d, e)

fit$p_mat[1:5, 1:5]          # SEcov p-values
fit$adj[1:5, 1:5]            # adjacency matrix (e edges connected)
fit$partial_corr[1:5, 1:5]   # SEcov partial-correlation matrix (edge weights)
```

---

## Using your own data

`mnetwork()` expects a normalized, sample-by-taxon matrix. For raw microbiome
counts, transform the data first. The paper uses a centered log-ratio (clr)
transformation with a pseudocount to handle zeros:

```r
# raw_counts: an n (samples) x d (taxa) matrix or data frame of integer counts

clr <- function(counts, pseudocount = 1) {
  x <- log(counts + pseudocount)
  sweep(x, 1, rowMeans(x), "-")   # subtract each sample's log-geometric-mean
}

OTU <- clr(raw_counts)
d   <- ncol(OTU)

set.seed(1)
fit <- mnetwork(OTU, d, e = d)
```

> **Note on run time.** `mnetwork()` fits two random forests for each of the
> `d(d − 1)/2` taxon pairs, so the computation grows quickly with the number of
> taxa. For large `d`, expect longer run times.

---

## Visualizing the network (optional)

The adjacency and partial-correlation matrices plug directly into
[`igraph`](https://cran.r-project.org/package=igraph). Coloring edges by the
sign of the conditional association reproduces the green (positive) / red
(negative) convention used in the paper.

```r
# install.packages("igraph")
library(igraph)

g <- graph_from_adjacency_matrix(fit$adj, mode = "undirected", diag = FALSE)

el <- as_edgelist(g, names = FALSE)
E(g)$color <- ifelse(fit$partial_corr[el] > 0, "forestgreen", "tomato")

plot(g, vertex.size = 6, vertex.label = NA, edge.width = 2)
```

---

## Helper function: `make_adj()`

`mnetwork()` calls `make_adj()` internally, but you can also use it on its own
to turn any p-value matrix into an adjacency matrix with a chosen number of
edges (the pairs with the smallest p-values).

```r
make_adj(pv_mat, idx)
```

| Argument | Description |
|----------|-------------|
| `pv_mat` | A square matrix whose entries are p-values. |
| `idx`    | The number of edges to connect. |

It returns a list with `est_adj` (the adjacency matrix) and `cut` (the p-value
cutoff that yields `idx` edges).

```r
res <- make_adj(pv_mat = fit$p_mat, idx = 30)  # keep the 30 most significant edges
res$est_adj
res$cut
```

---

## References

- Yunhua Xiang and Noah Simon (2020). *A flexible framework for non-parametric
  graphical modeling that accommodates machine learning.* International
  Conference on Machine Learning (ICML), PMLR, 10442–10451.
- Marvin N. Wright and Andreas Ziegler (2017). *ranger: A fast implementation of
  random forests for high-dimensional data in C++ and R.* Journal of Statistical
  Software, 77(1), 1–17.

---

## Maintainer

Hoseung Song (`hoseung@kaist.ac.kr`)
