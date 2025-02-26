## Cell-type Relationship Analysis Workflow Done Across Distances

`CRAWDAD` enables the characterization of spatial cell-type relationships across different length scales.

## Overview

`CRAWDAD` is a statistical framework that uses cell-type labeled spatial omics data to identify the colocalization or separation of cell types at different length scales. CRAWDAD identifies the spatial relationship of cell types in the tissue and the length scale in which they reach significance. It identifies groups of cells with similar spatial relationships patterns. Also, CRAWDAD subsets the cell types based on their proximity to others. Lastly, CRAWDAD compares multiple tissue samples using their cell-type spatial relationship patterns. Therefore, CRAWDAD is a powerful tool for tissue characterization and comparison.

<img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/img/overview_github.png?raw=true"/>

## Installation

To install `CRAWDAD`, we recommend using `remotes`:

``` r
require(remotes)
remotes::install_github('JEFworks-Lab/CRAWDAD')
```

## Visualization of cells in spleen pkhl and ksfb samples

``` r

library(crawdad)
library(tidyverse)
library(gridExtra)

## Load the spleen data of the pkhl sample 
data('pkhl')
## Load the spleen data of the ksfb sample 
data('ksfb')

all_celltypes <- unique(c(pkhl$celltype, ksfb$celltype))
colors <- rainbow(length(all_celltypes))
names(colors) <- all_celltypes

## Plot 1: pkhl dataset
p1 <- pkhl %>% 
  ggplot(aes(x,y,color = celltype)) +
  geom_point(size = .01) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(size = 1))) + 
  labs(title = 'pkhl') + 
  theme_minimal()

## Plot 2: ksfb dataset
p2 <- ksfb %>% 
  ggplot(aes(x,y,color = celltype)) +
  geom_point(size = .01) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(size = 1))) + 
  labs(title = 'ksfb') + 
  theme_minimal()

options(repr.plot.width = 20, repr.plot.height = 9)
grid.arrange(p1, p2, nrow = 1)
```
<img src="https://github.com/rafaeldossantospeixoto/sdk_analysis/blob/main/spleen/visualization.png?raw=true" height="510"/>

## Determine neighboring distance

For choosing neighboring distance for CRAWDAD analysis, please refer to the [CRAWDAD documentation](https://github.com/JEFworks-Lab/CRAWDAD/blob/main/docs/choosing_neighDist.md).

## Individual sample analysis

```r
## pkhl sample

## Convert dataframe to spatial points (SP)
cells1 <- crawdad::toSF(pos = pkhl[,c("x", "y")], cellTypes = pkhl$celltypes)
## Define the scales to analyze the data
scales1 <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
## Shuffle cells to create null background
shuffle_list1 <- crawdad:::makeShuffledCells(cells1,
                                             scales = scales1,
                                             perms = 3,
                                             ncores = 11,
                                             seed = 1,
                                             verbose = TRUE)
## Calculate the z-score for the cell-type pairs at different scales
results1 <- crawdad::findTrends(cells1,
                                neighDist = 50,
                                shuffleList = shuffle_list1,
                                ncores = 7,
                                verbose = TRUE,
                                returnMeans = FALSE)
dat_pkhl <- crawdad::meltResultsList(results1, withPerms = TRUE)

## Calculate the z-score for the multiple-test correction
zsig1 <- correctZBonferroni(dat_pkhl)
## Summary visualization
vizColocDotplot(dat_pkhl, zSigThresh = zsig1, zScoreLimit = 2*zsig1, 
                dotSizes = c(3,15)) +
  theme(axis.text.x = element_text(angle = 35, h = 0))

## ksfb sample

## Convert dataframe to spatial points (SP)
cells2 <- crawdad::toSF(pos = ksfb[,c("x", "y")], cellTypes = ksfb$celltypes)
## Define the scales to analyze the data
scales2 <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
## Shuffle cells to create null background
shuffle_list2 <- crawdad:::makeShuffledCells(cells2,
                                             scales = scales2,
                                             perms = 3,
                                             ncores = 7,
                                             seed = 1,
                                             verbose = TRUE)
## Calculate the z-score for the cell-type pairs at different scales
results2 <- crawdad::findTrends(cells2,
                                neighDist = 50,
                                shuffleList = shuffle_list2,
                                ncores = 7,
                                verbose = TRUE,
                                returnMeans = FALSE)
dat_ksfb <- crawdad::meltResultsList(results2, withPerms = TRUE)

## Calculate the z-score for the multiple-test correction
zsig2 <- correctZBonferroni(dat_ksfb)
## Summary visualization

vizColocDotplot(dat_ksfb, zSigThresh = zsig2, zScoreLimit = 2*zsig2, 
                dotSizes = c(3,15)) +
  theme(axis.text.x = element_text(angle = 35, h = 0))
```

<img src="https://github.com/rafaeldossantospeixoto/sdk_analysis/blob/main/spleen/CRAWDAD.png?raw=true" height="510"/>

## Multiple sample analysis
```r
## Add id to datasets
dat_pkhl$id <- 'pkhl'
dat_ksfb$id <- 'ksfb'
## Calculate auc
auc_samples <- calculateAUC(list(dat_pkhl, dat_ksfb))

options(repr.plot.width = 10, repr.plot.height = 9)
vizVarianceSamples(auc_samples)
```
<img src="https://github.com/rafaeldossantospeixoto/sdk_analysis/blob/main/spleen/variance.png?raw=true" height="510"/>

## Visualize trends for one cell-type pair
```r
reference_var <- 'Ki67 proliferating'
neighbor_var <- 'Fol B cells'

d <- dplyr::bind_rows(list(dat_ksfb, dat_pkhl))
d %>%
  filter(reference == reference_var) %>%
  filter(neighbor == neighbor_var) %>%
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig1)

```
<img src="https://github.com/rafaeldossantospeixoto/sdk_analysis/blob/main/spleen/Ki67 proliferating-Fol B cells_relationships_plot.png?raw=true" height="510"/>
  
```r

reference_var <- 'B cells, red pulp'
neighbor_var <- 'Fol B cells'

d <- dplyr::bind_rows(list(dat_ksfb, dat_pkhl))
d %>%
  filter(reference == reference_var) %>%
  filter(neighbor == neighbor_var) %>%
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig1)

```
<img src="https://github.com/rafaeldossantospeixoto/sdk_analysis/blob/main/spleen/B cells, red pulp-Fol B cells_relationships_plot.png?raw=true" height="510"/>
  
## Visualize cell-type relationships
```r
plot_celltypes <- function(celltypes) {
  options(repr.plot.width = 18, repr.plot.height = 9)
  
  ## Combine unique cell types from both datasets
  all_celltypes <- unique(c(pkhl$celltypes, ksfb$celltypes))
  
  ## Define consistent colors for all datasets
  possible_colors <- c('darkblue', 'red', 'cyan', 'yellow', 'maroon')
  colors <- setNames(
    rep('#F0F0F0', length(all_celltypes)),
    all_celltypes
  )
  for (i in seq_along(celltypes)) {
    if (i <= length(possible_colors)) {
      colors[celltypes[i]] <- possible_colors[i]
    }
  }
  
  ## Plot for pkhl
  p1 <- pkhl %>%
    mutate(celltype_highlight = celltypes %in% celltypes) %>%
    arrange(celltype_highlight, celltypes) %>%
    ggplot() +
    geom_point(aes(x, y, color = celltypes), size = .01) +
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size = 1))) + 
    labs(title = 'pkhl') + 
    theme_minimal()
  
  ## Plot for ksfb
  p2 <- ksfb %>%
    mutate(celltype_highlight = celltypes %in% celltypes) %>%
    arrange(celltype_highlight, celltypes) %>%
    ggplot() +
    geom_point(aes(x, y, color = celltypes), size = .01) +
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size = 1))) + 
    labs(title = 'ksfb') + 
    theme_minimal()
  
  ## Arrange plots
  grid.arrange(p1, p2, nrow = 1)
}

## Example usage
reference_cell <- 'Ki67 proliferating'
neighbor_cell <- 'Fol B cells'

plot_celltypes(c(reference_cell, neighbor_cell))
```
<img src="https://github.com/rafaeldossantospeixoto/sdk_analysis/blob/main/spleen/celltype_relationship_visualizations/Ki67 proliferating_Fol B cells_plot.png?raw=true" height="510"/>


```r
reference_cell <- 'B cells, red pulp'
neighbor_cell <- 'Fol B cells'

plot_celltypes(c(reference_cell, neighbor_cell))
```
<img src="https://github.com/rafaeldossantospeixoto/sdk_analysis/blob/main/spleen/celltype_relationship_visualizations/B cells, red pulp_Fol B cells_plot.png?raw=true" height="510"/>
  

More details can be found in the tutorials.

## Tutorials
- [`CRAWDAD` applied to simulated data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/1_simulations.md)
- [`CRAWDAD` applied to a mouse embryo seqfish data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/2_seqfish.md)
- [`CRAWDAD` applied to a mouse cerebellum slideseq data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/2_slideseq.md)
- [`CRAWDAD` applied to a human spleen codex data](https://github.com/JEFworks/CRAWDAD/blob/main/docs/3_spleen.md)

## Instalation

CRAWDAD was tested on R version 4.3.0 (2023-04-21) -- "Already Tomorrow".
Installation time is 1 minute in a Mac Pro Computer.
