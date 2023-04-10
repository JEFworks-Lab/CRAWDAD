# CRAWDAD

<!-- badges: start -->

<!-- badges: end -->

`CRAWDAD`: Cell Relationship Analysis Workflow Done Across Distances

<img src=""/>

The overall approach is now published in [link]

## Overview

<img src=""/>

## Example

```{r}

data(sim)

## visualize
plt <- crawdad::vizAllClusters(cells = sim,
                               coms = sim$type,
                               title = "sim",
                               axisAdj = 1, s = 6, a = 0.5) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 1))
plt

## convert to SP
cells <- crawdad::toSP(pos = sim[,c("x", "y")],
                        celltypes = sim$type)
cells

```

```{r}

oldw <- getOption("warn")
options(warn = -1)

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                          resolutions = c(100, 200, 500, 800, 1000),
                          perms = 1,
                          ncores = ncores,
                          seed = 1,
                          verbose = TRUE)

options(warn = oldw)

```

```{r}

oldw <- getOption("warn")
options(warn = -1)

## find trends, passing background as parameter
results <- crawdad::findTrends(cells,
                        dist = 100,
                        shuffle.list = shuffle.list,
                        ncores = ncores,
                        verbose = TRUE)

options(warn = oldw)

```

```{r}

plt <- crawdad::vizTrends(dat = results) +
  ggplot2::scale_x_log10()
plt 

plt <- crawdad::vizTrends.heatmap(dat = results, annotation = TRUE,
                                  z_limit = 20, # Z score limits (default +/- 20). Winsorize values to these limits
                                  sig.thresh = 1.96, # Z score significance threshold (default: 1.96). Non-significant values < 1.96 and > -1.96 colored white
                                  ncols = 4 # specify number of columns in facet wrap (default: 4)
                                  )
plt 

```

<img src=""/>

More details can be found in the tutorials.

## Tutorials
- [Getting started with `CRAWDAD`](tutorial.md)
- [Simulated data analysis](1_simulations.md)
- [Slide-seq analysis](2_slideseq.md)
- [seqFISH analysis](2_seqfish.md)
- [Spleen analysis](3_spleen.md)

## Installation

To install `CRAWDAD`, we recommend using `remotes`:

``` r
require(remotes)
remotes::install_github('JEFworks-Lab/CRAWDAD')
```

## Contributing

We welcome any bug reports, enhancement requests, general questions, and other contributions. To submit a bug report or enhancement request, please use the `CRAWDAD` GitHub issues tracker. For more substantial contributions, please fork this repo, push your changes to your fork, and submit a pull request with a good commit message.
