# Cell-type Relationship Analysis Workflow Done Across Distances

`CRAWDAD` enables the characterization of spatial cell-type relationships across different length scales.

## Overview

`CRAWDAD` is a statistical framework that uses cell-type labeled spatial omics data to identify the colocalization or separation of cell types at different length scales. CRAWDAD identifies the spatial relationship of cell types in the tissue and the length scale in which they reach significance. It identifies groups of cells with similar spatial relationships patterns. Also, CRAWDAD subsets the cell types based on their proximity to others. Lastly, CRAWDAD compares multiple tissue samples using their cell-type spatial relationship patterns. Therefore, CRAWDAD is a powerful tool for tissue characterization and comparison.

\
<img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/img/overview_github.png?raw=true"/>

\
More information can be found in our publication: [**Characterizing cell-type spatial relationships across length scales in spatially resolved omics data.** Rafael dos Santos Peixoto, Brendan F. Miller, Maigan A. Brusko, Gohta Aihara, Lyla Atta, Manjari Anant, Mark A. Atkinson, Todd M. Brusko, Clive H. Wasserfall, Jean Fan.](https://doi.org/10.1038/s41467-024-55700-1)

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
cells <- crawdad::toSF(pos = pkhl[,c("x", "y")], cellTypes = pkhl$celltypes)
## define the scales to analyze the data
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
## shuffle cells to create null background
shuffle_list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 3,
                                            ncores = 7,
                                            seed = 1,
                                            verbose = TRUE)
## calculate the zscore for the cell-type pairs at different scales
results <- crawdad::findTrends(cells,
                               neighDist = 50,
                               shuffleList = shuffle_list,
                               ncores = 7,
                               verbose = TRUE,
                               returnMeans = FALSE)
dat <- crawdad::meltResultsList(results, withPerms = TRUE)
## calculate the zscore for the multiple-test correction
zsig <- correctZBonferroni(dat)
## summary visualization
vizColocDotplot(dat, zSigThresh = zsig, zScoreLimit = 2*zsig, 
                dotSizes = c(3,15)) +
  theme(axis.text.x = element_text(angle = 35, h = 0))
```

<img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/img/coloc.png?raw=true" height="510"/>

``` r
## visualize trend for one cell-type pair
dat %>% 
  mutate(id = 'pkhl') %>% 
  filter(reference == 'Podoplanin') %>% 
  filter(neighbor == 'CD4 Memory T cells') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig)
```

<img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/img/trend.png?raw=true" height="325"/>

More details can be found in the tutorials.

## Tutorials
- [`CRAWDAD` applied to simulated data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/1_simulations.md)
- [`CRAWDAD` applied to a mouse embryo seqfish data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/2_seqfish.md)
- [`CRAWDAD` applied to a mouse cerebellum slideseq data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/2_slideseq.md)
- [`CRAWDAD` applied to a human spleen codex data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/3_spleen.md)

## Instalation

CRAWDAD was tested on R version 4.3.0 (2023-04-21) -- "Already Tomorrow".
Installation time is 1 minute in a Mac Pro Computer.
