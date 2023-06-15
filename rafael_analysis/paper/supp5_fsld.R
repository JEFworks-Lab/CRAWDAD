
# Load package and data ---------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 1

fsld <- read.csv2(file = paste0(here::here(), "/data/spleen/FSLD.meta.csv.gz"), row.names = 1)
fsld <- fsld[,c("x", "y", "celltypes_folBcombined")]
## make sure the coordinates are numeric
fsld <- fsld %>%
  dplyr::mutate_at(vars(x, y), as.numeric)

colnames(fsld) <- c('x', 'y', 'celltypes')

sort(table(fsld$celltypes))
set.seed(8) 
cols <- sample(rainbow(length(unique(fsld$celltypes))))
names(cols) <- unique(fsld$celltypes)
cols['indistinct'] <- 'grey'

table(fsld$celltypes)
vi <- fsld$celltypes %in% c('Podoplanin', 'CD4 Memory T cells', 'Fol B cells')
ggplot(fsld, aes(x=x, y=y, col=celltypes)) + geom_point(size=0.2, alpha=0.5) + theme_void() +
  scale_color_manual(values=cols[c('Podoplanin', 'CD4 Memory T cells', 'Fol B cells')],
                     na.value = 'lightgrey') + theme(legend.position = 'none')



# CRAWDAD -----------------------------------------------------------------

## read in Brendan's old results
# results <- readRDS("/Users/jeanfan/Downloads/fsld.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds")

## convert to SP
cells <- crawdad::toSP(pos = fsld[,c("x", "y")],
                       celltypes = fsld$celltypes)

scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 3000, 6000)
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)

## find trends, passing background as parameter
results <- crawdad::findTrends(cells,
                               dist = 100,
                               shuffle.list = shuffle.list,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)

save(results, shuffle.list, file="rafael_analysis/paper/supp5_fsld.RData")
load("rafael_analysis/paper/supp5_fsld.RData")

dat <- crawdad::meltResultsList(results, withPerms = TRUE)
head(dat)

vizAllClusters(cells, fsld$celltypes)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')




# Subsetting --------------------------------------------------------------

binomMat <- crawdad::binomialTestMatrix(cells,
                                        neigh.dist = 100,
                                        ncores = ncores,
                                        verbose = TRUE)

head(binomMat)

saveRDS(binomMat, file="rafael_analysis/paper/supp5_fsld_binomMat.RDS")
binomMat <- readRDS("rafael_analysis/paper/supp5_fsld_binomMat.RDS")

subset.list <- crawdad::selectSubsets(binomMat,
                                      cells$celltypes,
                                      sub.type = "near",
                                      sub.thresh = 0.05,
                                      ncores = ncores,
                                      verbose = TRUE)

saveRDS(subset.list, file="rafael_analysis/paper/supp5_fsld_subset.RDS")
subset.list <- readRDS("rafael_analysis/paper/supp5_fsld_subset.RDS")

results.subsets <- crawdad::findTrends(cells,
                                       dist = 100,
                                       shuffle.list = shuffle.list,
                                       subset.list = subset.list,
                                       ncores = ncores,
                                       verbose = TRUE,
                                       returnMeans = FALSE)
## 8.0865 hours to run
results.subsets
saveRDS(results.subsets, file="rafael_analysis/paper/supp5_fsld_results.subsets.RDS")
results.subsets <- readRDS("rafael_analysis/paper/supp5_fsld_results.subsets.RDS")




# Viz subsetting ----------------------------------------------------------

## cts of interest
plt <- crawdad::vizAllClusters(cells = fsld,
                               coms = as.factor(fsld$celltypes),
                               ofInterest = c("CD4 Memory T cells", 
                                              "Fol B cells", "Podoplanin"),
                               nacol = crawdad:::transparentCol(color = "white", 
                                                                percent = 100),
                               axisAdj = 1, s = 3, a = 0.5) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), 
                                                 ncol = 1))
plt

## CD$ Near Fol B
## make a temporary annotation factor with annotations for cells that are of a specific subsets
annots_temp <- crawdad::selectLabels(df = fsld,
                                     com = fsld$celltypes,
                                     subset_list = subset.list,
                                     cellIDs = c("Fol B cells", "Podoplanin"),
                                     subsetIDs = c("CD4 Memory T cells_near_Fol B cells"))

## visualize the subset only
plt <- crawdad::vizAllClusters(cells = fsld,
                               coms = annots_temp,
                               ofInterest = NULL,
                               title = "fsld",
                               nacol = crawdad:::transparentCol(color = "white", percent = 100),
                               axisAdj = 1, s = 3, a = 0.5) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 1))
plt



# Figure 3f ---------------------------------------------------------------

## CD4+ memory near follicle B cells colocalize with podoplanin trend
dats <- crawdad::meltResultsList(results.subsets, withPerms = TRUE)

d1 <- dats[grepl(pattern = "CD4 Memory T cells_near_Fol B cells", 
                 dats$reference) & dats$neighbor %in% c("Fol B cells", "Podoplanin"),]
plt <- vizTrends(dat = d1, facet = FALSE, id = "neighbor", title = "CD4 Memory T cells near Fol B cells") +
  ggplot2::scale_x_log10()
# ggplot2::theme(legend.position="none")
plt

dat_filter <- d1 %>% 
  filter(neighbor == 'Podoplanin')
p <- vizTrends(dat_filter, lines = T, withPerms = T)
p
pdf('rafael_analysis/paper/supp5_fsld_follB_trend.pdf')
p 
dev.off()


# Figure 3g ---------------------------------------------------------------

## CD4+ memory near red bulb colocalize with podoplanin trend
d2 <- dats[grepl(pattern = "CD4 Memory T cells_near_B cells, red pulp", 
                 dats$reference) & dats$neighbor %in% c("Fol B cells", "Podoplanin"),]
plt <- vizTrends(dat = d2, facet = FALSE, id = "neighbor", title = "CD4 Memory T cells near Neutrophils/Monocytes") +
  ggplot2::scale_x_log10()
# ggplot2::theme(legend.position="none")
plt

dat_filter <- d2 %>% 
  filter(neighbor == 'Podoplanin')
p <- vizTrends(dat_filter, lines = T, withPerms = T)
p
pdf('rafael_analysis/paper/supp5_fsld_redpulb_trend.pdf')
p 
dev.off()



