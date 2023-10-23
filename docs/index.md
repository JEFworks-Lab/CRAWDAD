# CRAWDAD

<!-- badges: start -->

<!-- badges: end -->

`CRAWDAD` enables the characterization of cell-type relationships across different scales.

<p align="center">
  <img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/img/CRAWDAD_logo.png?raw=true" height="200"/>
</p>

## Overview

`CRAWDAD` is a statistical framework that uses labeled spatial omics data to identify the colocalization or separation of cell types at different scales. CRAWDAD identifies regions where multiple cells colocalize, the scale of such colocalization, and also subsets the cell types based on their proximity to others. Therefore, CRAWDAD is a powerful tool for tissue characterization and comparison.

<img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/img/CRAWDAD_workflow.png?raw=true"/>

## Installation

To install `CRAWDAD`, we recommend using `remotes`:

``` r
require(remotes)
remotes::install_github('JEFworks-Lab/CRAWDAD')
```

## Example

``` r
library(crawdad)
library(tidyverse)
## load the spleen data of the pkhl sample 
data('pkhl')
## convert dataframe to spatial points (SP)
cells <- crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
## define the scales to analyze the data
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
## shuffle cells to create null background
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 3,
                                            ncores = 7,
                                            seed = 1,
                                            verbose = TRUE)
## calculate the zscore for the cell-type pairs at different scales
results <- crawdad::findTrends(cells,
                               dist = 100,
                               shuffle.list = shuffle.list,
                               ncores = 7,
                               verbose = TRUE,
                               returnMeans = FALSE)
dat <- crawdad::meltResultsList(results, withPerms = TRUE)
## calculate the zscore for the multiple-test correction
ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)
## summary visualization
vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig) +
  theme(axis.text.x = element_text(angle = 35, h = 0))
```

<img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/img/coloc.png?raw=true" height="600"/>

``` r
## visualize trend for one cell-type pair
dat %>% 
  filter(reference == 'Podoplanin') %>% 
  filter(neighbor == 'CD4 Memory T cells') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
```

<img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/img/trend.png?raw=true" height="350"/>

More details can be found in the tutorials.

## Tutorials
- [`CRAWDAD` applied to simulated data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/1_simulations.md)
- [`CRAWDAD` applied to a mouse embryo seqfish data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/2_seqfish.md)
- [`CRAWDAD` applied to a mouse cerebellum slideseq data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/2_slideseq.md)
- [`CRAWDAD` applied to a human spleen codex data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/3_spleen.md)

## Installation

To install `CRAWDAD`, we recommend using `remotes`:

``` r
require(remotes)
remotes::install_github('JEFworks-Lab/CRAWDAD')
```

## Contributing

We welcome any bug reports, enhancement requests, general questions, and other contributions. To submit a bug report or enhancement request, please use the `CRAWDAD` GitHub issues tracker. For more substantial contributions, please fork this repo, push your changes to your fork, and submit a pull request with a good commit message.
