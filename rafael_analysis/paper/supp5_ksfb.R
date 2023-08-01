
# Load package and data ---------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 7

ksfb <- read.csv2(file = paste0(here::here(), "/data/spleen/KSFB.meta.csv.gz"), row.names = 1)
ksfb <- ksfb[,c("x", "y", "celltypes_folBcombined")]
## make sure the coordinates are numeric
ksfb <- ksfb %>%
  dplyr::mutate_at(vars(x, y), as.numeric)

colnames(ksfb) <- c('x', 'y', 'celltypes')

sort(table(ksfb$celltypes))
set.seed(8) 
cols <- sample(rainbow(length(unique(ksfb$celltypes))))
names(cols) <- unique(ksfb$celltypes)
cols['indistinct'] <- 'grey'

table(ksfb$celltypes)
vi <- ksfb$celltypes %in% c('Podoplanin', 'CD4 Memory T cells', 'Fol B cells')
ggplot(ksfb, aes(x=x, y=y, col=celltypes)) + geom_point(size=0.2, alpha=0.5) + theme_void() +
  scale_color_manual(values=cols[c('Podoplanin', 'CD4 Memory T cells', 'Fol B cells')],
                     na.value = 'lightgrey') + theme(legend.position = 'none')



# CRAWDAD -----------------------------------------------------------------

## read in Brendan's old results
# results <- readRDS("/Users/jeanfan/Downloads/ksfb.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds")

## convert to SP
cells <- crawdad::toSP(pos = ksfb[,c("x", "y")],
                       celltypes = ksfb$celltypes)

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

save(results, shuffle.list, file="rafael_analysis/paper/supp5_ksfb.RData")
load("rafael_analysis/paper/supp5_ksfb.RData")

dat <- crawdad::meltResultsList(results, withPerms = TRUE)
head(dat)



## Multiple-test correction
ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)



vizAllClusters(cells, ksfb$celltypes)

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

saveRDS(binomMat, file="rafael_analysis/paper/supp5_ksfb_binomMat.RDS")
binomMat <- readRDS("rafael_analysis/paper/supp5_ksfb_binomMat.RDS")

subset.list <- crawdad::selectSubsets(binomMat,
                                      cells$celltypes,
                                      sub.type = "near",
                                      sub.thresh = 0.05,
                                      ncores = ncores,
                                      verbose = TRUE)

saveRDS(subset.list, file="rafael_analysis/paper/supp5_ksfb_subset.RDS")
subset.list <- readRDS("rafael_analysis/paper/supp5_ksfb_subset.RDS")

results.subsets <- crawdad::findTrends(cells,
                                       dist = 100,
                                       shuffle.list = shuffle.list,
                                       subset.list = subset.list,
                                       ncores = ncores,
                                       verbose = TRUE,
                                       returnMeans = FALSE)
## 8.0865 hours to run
results.subsets
saveRDS(results.subsets, file="rafael_analysis/paper/supp5_ksfb_results.subsets.RDS")
results.subsets <- readRDS("rafael_analysis/paper/supp5_ksfb_results.subsets.RDS")



## subsets
dats <- crawdad::meltResultsList(results.subsets, withPerms = TRUE)

## Multiple-test correction
ntestss <- length(unique(dats$reference)) * length(unique(dats$neighbor))
psigs <- 0.05/ntestss
zsigs <- round(qnorm(psigs/2, lower.tail = F), 2)




# Viz subsetting ----------------------------------------------------------

## cts of interest
plt <- crawdad::vizAllClusters(cells = ksfb,
                               coms = as.factor(ksfb$celltypes),
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
annots_temp <- crawdad::selectLabels(df = ksfb,
                                     com = ksfb$celltypes,
                                     subset_list = subset.list,
                                     cellIDs = c("Fol B cells", "Podoplanin"),
                                     subsetIDs = c("CD4 Memory T cells_near_Fol B cells"))

## visualize the subset only
plt <- crawdad::vizAllClusters(cells = ksfb,
                               coms = annots_temp,
                               ofInterest = NULL,
                               title = "ksfb",
                               nacol = crawdad:::transparentCol(color = "white", percent = 100),
                               axisAdj = 1, s = 3, a = 0.5) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 1))
plt



# Figure 3f ---------------------------------------------------------------

## CD4+ memory near follicle B cells colocalize with podoplanin trend
d1 <- dats[grepl(pattern = "CD4 Memory T cells_near_Fol B cells", 
                 dats$reference) & dats$neighbor %in% c("Fol B cells", "Podoplanin"),]
plt <- vizTrends(dat = d1, facet = FALSE, id = "neighbor", title = "CD4 Memory T cells near Fol B cells") +
  ggplot2::scale_x_log10()
# ggplot2::theme(legend.position="none")
plt

dat_filter <- d1 %>% 
  filter(neighbor == 'Podoplanin')
p <- vizTrends(dat_filter, lines = T, withPerms = T, sig.thresh = zsigs)
p
pdf('rafael_analysis/paper/supp5_ksfb_follB_trend.pdf')
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
p <- vizTrends(dat_filter, lines = T, withPerms = T, sig.thresh = zsigs)
p
pdf('rafael_analysis/paper/supp5_ksfb_redpulb_trend.pdf')
p 
dev.off()



# Figure 3d ---------------------------------------------------------------

p <- vizColocDotplot(dat_filtered, reorder = TRUE, zsig.thresh = zsig, zscore.limit = zsig*2) + 
  scale_x_discrete(position = 'bottom') +
  scale_y_discrete(position = 'right') +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 90, h = 1),
        legend.box = 'vertical')
p
pdf('rafael_analysis/paper/supp5c_ksfb.pdf', height = 7.7, width = 6.05)
p 
dev.off()

## reorder according to pkhl
ct_ordered <- readRDS("rafael_analysis/paper/ct_ordered_spleen.RDS")
p <- vizColocDotplot(dat_filtered, reorder = TRUE, zsig.thresh = zsig, zscore.limit = zsig*2) + 
  scale_x_discrete(position = 'bottom', limits=ct_ordered) +
  scale_y_discrete(position = 'right', limits=ct_ordered) +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 90, h = 1),
        legend.box = 'vertical')
p
pdf('rafael_analysis/paper/supp5c_ksfb_ordered_pkhl.pdf', height = 7.7, width = 6.05)
p 
dev.off()



# Figure 3e ---------------------------------------------------------------

results.subsets <- readRDS("rafael_analysis/paper/supp5_ksfb_results.subsets.RDS")
subset.list <- readRDS("rafael_analysis/paper/supp5_ksfb_subset.RDS")
subcells <- subset.list
ssample <- ksfb

## visualize names
names(subcells)
names(subcells)[grepl('CD4 Memory T cells_near', names(subcells))]
names(subcells)[grepl('Podoplanin_near_', names(subcells))]

## select subsets
c_cd4_folb <- subcells[["CD4 Memory T cells_near_Fol B cells"]]
c_cd4_bcell <- subcells[["CD4 Memory T cells_near_B cells, red pulp"]]
c_podoplanin <- which(ssample$celltypes == "Podoplanin")

## rewrite types
ssample[c_cd4_folb, ]$celltypes <- 'CD4 Memory T cells_near_Fol B cells'
ssample[c_cd4_bcell, ]$celltypes <- 'CD4 Memory T cells_near_B cells, red pulp'

## define cts of interest and colors
cts_interest <- ssample$celltypes %in% c('Podoplanin', 
                                         'CD4 Memory T cells_near_Fol B cells', 
                                         'CD4 Memory T cells_near_B cells, red pulp')
cols_ssample <- c("#FFFF00", "#FF0080", "#0000FF")
names(cols_ssample) <-  c('Podoplanin', 
                          'CD4 Memory T cells_near_Fol B cells', 
                          'CD4 Memory T cells_near_B cells, red pulp')

df_cts <- ssample %>% 
  filter(celltypes %in% names(cols_ssample))

df_bkg <-  ssample %>% 
  filter(!(celltypes %in% names(cols_ssample)))

p <- ggplot() + 
  geom_point(data = df_bkg, aes(x=x, y=y), 
             color = 'lightgrey', size=1.5, alpha=0.5) + 
  
  geom_point(data = df_cts[df_cts$celltypes == 'Podoplanin', ], 
             aes(x=x, y=y), color = '#FFFF00', size=1.5, alpha=0.5) + 
  
  geom_point(data = df_cts[df_cts$celltypes == 'CD4 Memory T cells_near_Fol B cells', ], 
             aes(x=x, y=y), color = '#FF0080', size=1.5, alpha=0.5) + 
  
  geom_point(data = df_cts[df_cts$celltypes == 'CD4 Memory T cells_near_B cells, red pulp', ], 
             aes(x=x, y=y), color = '#0000FF', size=1.5, alpha=0.5) + 
  theme_void() + 
  theme(legend.position="none")
p
png('rafael_analysis/paper/supp5_ksfb_spat_cts.png', height = 1024, width = 1024)
p 
dev.off()





# Replication across samples ----------------------------------------------

## Podoplanin near follicle B cells (white pupl) colocalize with CD4 trend
d1 <- dats[grepl(pattern = "Podoplanin_near_Fol B cells", 
                 dats$reference) & dats$neighbor %in% c("CD4 Memory T cells"),]
plt <- vizTrends(dat = d1, facet = FALSE, id = "neighbor", title = "Podoplanin cells near Fol B cells") +
  ggplot2::scale_x_log10()
# ggplot2::theme(legend.position="none")
plt

dat_filter <- d1 %>% 
  filter(neighbor == 'CD4 Memory T cells')
p <- vizTrends(dat_filter, lines = T, withPerms = T, sig.thresh = zsigs)
p
pdf('rafael_analysis/paper/replication/ksfb_podofolB_CD4.pdf')
p 
dev.off()



## Podoplanin near follicle B cells (red pupl) colocalize with CD4 trend
d1 <- dats[grepl(pattern = "Podoplanin_near_B cells, red pulp", 
                 dats$reference) & dats$neighbor %in% c("CD4 Memory T cells"),]
plt <- vizTrends(dat = d1, facet = FALSE, id = "neighbor", title = "Podoplanin cells near B cells, red pulp") +
  ggplot2::scale_x_log10()
# ggplot2::theme(legend.position="none")
plt

dat_filter <- d1 %>% 
  filter(neighbor == 'CD4 Memory T cells')
p <- vizTrends(dat_filter, lines = T, withPerms = T, sig.thresh = zsigs)
p
pdf('rafael_analysis/paper/replication/ksfb_podoredpulp_CD4.pdf')
p 
dev.off()



# Subsets Podoplanin ------------------------------------------------------

results.subsets <- readRDS("rafael_analysis/paper/supp5_ksfb_results.subsets.RDS")
subset.list <- readRDS("rafael_analysis/paper/supp5_ksfb_subset.RDS")
subcells <- subset.list
ssample <- ksfb

## visualize names
names(subcells)[grepl('Podoplanin_near_', names(subcells))]

## select subsets
c_podo_folb <- subcells[["Podoplanin_near_Fol B cells"]]
c_podo <- which(ssample$celltypes == "Podoplanin")
c_podo_notfolb <- c_podo[!c_podo %in% c_podo_folb]
c_cd4 <- which(ssample$celltypes == "CD4 Memory T cells")

## rewrite types
ssample[c_podo_folb, ]$celltypes <- 'Podoplanin near Fol B cells'
ssample[c_podo_notfolb, ]$celltypes <- 'Podoplanin not near Fol B cells'

## define cts of interest and colors
cts_interest <- ssample$celltypes %in% c('CD4 Memory T cells', 
                                         'Podoplanin near Fol B cells', 
                                         'Podoplanin not near Fol B cells')

cols_ssample <- c("#FFFF00", "#FF0080", "#0000FF")
names(cols_ssample) <-  c('CD4 Memory T cells', 
                          'Podoplanin near Fol B cells', 
                          'Podoplanin not near Fol B cells')

df_cts <- ssample %>% 
  filter(celltypes %in% names(cols_ssample))

df_bkg <-  ssample %>% 
  filter(!(celltypes %in% names(cols_ssample)))

p <- ggplot() + 
  geom_point(data = df_bkg, aes(x=x, y=y), 
             color = 'lightgrey', size=1.5, alpha=0.5) + 
  
  geom_point(data = df_cts[df_cts$celltypes == 'CD4 Memory T cells', ], 
             aes(x=x, y=y), color = '#FF0080', size=1.5, alpha=0.5) + 
  
  geom_point(data = df_cts[df_cts$celltypes == 'Podoplanin near Fol B cells', ], 
             aes(x=x, y=y), color = '#FFFF00', size=1.5, alpha=0.5) + 
  
  geom_point(data = df_cts[df_cts$celltypes == 'Podoplanin not near Fol B cells', ], 
             aes(x=x, y=y), color = '#0000FF', size=1.5, alpha=0.5) +
  theme_void() + 
  theme(legend.position="none")
p
png('rafael_analysis/paper/subsets/podo_ksfb_spat_cts.png', height = 1024, width = 1024)
p 
dev.off()






# Subsetting 500 --------------------------------------------------------------

cells <- crawdad::toSP(pos = ksfb[,c("x", "y")],
                       celltypes = ksfb$celltypes)
load("rafael_analysis/paper/supp5_ksfb.RData")


binomMat <- crawdad::binomialTestMatrix(cells,
                                        neigh.dist = 500,
                                        ncores = ncores,
                                        verbose = TRUE)

head(binomMat)

saveRDS(binomMat, file="rafael_analysis/paper/supp5_ksfb_binomMat500.RDS")
binomMat <- readRDS("rafael_analysis/paper/supp5_ksfb_binomMat500.RDS")

subset.list <- crawdad::selectSubsets(binomMat,
                                      cells$celltypes,
                                      sub.type = "near",
                                      sub.thresh = 0.05,
                                      ncores = ncores,
                                      verbose = TRUE)

saveRDS(subset.list, file="rafael_analysis/paper/supp5_ksfb_subset500.RDS")
subset.list <- readRDS("rafael_analysis/paper/supp5_ksfb_subset500.RDS")

results.subsets <- crawdad::findTrends(cells,
                                       dist = 100,
                                       shuffle.list = shuffle.list,
                                       subset.list = subset.list,
                                       ncores = ncores,
                                       verbose = TRUE,
                                       returnMeans = FALSE)
## 8.0865 hours to run
results.subsets
saveRDS(results.subsets, file="rafael_analysis/paper/supp5_ksfb_results500.subsets.RDS")
results.subsets <- readRDS("rafael_analysis/paper/supp5_ksfb_results500.subsets.RDS")



## subsets
dats <- crawdad::meltResultsList(results.subsets, withPerms = TRUE)

## Multiple-test correction
ntestss <- length(unique(dats$reference)) * length(unique(dats$neighbor))
psigs <- 0.05/ntestss
zsigs <- round(qnorm(psigs/2, lower.tail = F), 2)




# Spat cts 500 ------------------------------------------------------

subset.list <- readRDS("rafael_analysis/paper/supp5_ksfb_subset500.RDS")
subcells <- subset.list
ssample <- ksfb

## visualize names
names(subcells)[grepl('Podoplanin_near_', names(subcells))]

## select subsets
c_podo_folb <- subcells[["Podoplanin_near_Fol B cells"]]
c_podo <- which(ssample$celltypes == "Podoplanin")
c_podo_notfolb <- c_podo[!c_podo %in% c_podo_folb]
c_cd4 <- which(ssample$celltypes == "CD4 Memory T cells")

## rewrite types
ssample[c_podo_folb, ]$celltypes <- 'Podoplanin near Fol B cells'
ssample[c_podo_notfolb, ]$celltypes <- 'Podoplanin not near Fol B cells'

## define cts of interest and colors
cts_interest <- ssample$celltypes %in% c('CD4 Memory T cells', 
                                         'Podoplanin near Fol B cells', 
                                         'Podoplanin not near Fol B cells')

cols_ssample <- c("#FFFF00", "#FF0080", "#0000FF")
names(cols_ssample) <-  c('CD4 Memory T cells', 
                          'Podoplanin near Fol B cells', 
                          'Podoplanin not near Fol B cells')

df_cts <- ssample %>% 
  filter(celltypes %in% names(cols_ssample))

df_bkg <-  ssample %>% 
  filter(!(celltypes %in% names(cols_ssample)))

p <- ggplot() + 
  geom_point(data = df_bkg, aes(x=x, y=y), 
             color = 'lightgrey', size=1.5, alpha=0.5) + 
  
  geom_point(data = df_cts[df_cts$celltypes == 'CD4 Memory T cells', ], 
             aes(x=x, y=y), color = '#FF0080', size=1.5, alpha=0.5) + 
  
  geom_point(data = df_cts[df_cts$celltypes == '
                           ', ], 
             aes(x=x, y=y), color = '#FFFF00', size=1.5, alpha=0.5) + 
  
  geom_point(data = df_cts[df_cts$celltypes == 'Podoplanin not near Fol B cells', ], 
             aes(x=x, y=y), color = '#0000FF', size=1.5, alpha=0.5) +
  theme_void() + 
  theme(legend.position="none")
p
png('rafael_analysis/paper/subsets500/podo_ksfb_spat_cts_500.png', height = 1024, width = 1024)
p 
dev.off()



# Trends podo 500 ----------------------------------------------

load("rafael_analysis/paper/supp5_ksfb.RData")
cells <- crawdad::toSP(pos = ksfb[,c("x", "y")],
                       celltypes = ksfb$celltypes)
subset.list <- readRDS("rafael_analysis/paper/supp5_ksfb_subset500.RDS")

## change cell types
subset_names <- names(subset.list)
subset_names[grepl("Podoplanin_near_Fol B cells", subset_names)] <- "temp"
subset_names[grepl("^Podoplanin_near_", subset_names)] <- "Podoplanin_not_near_Fol B cells"
subset_names[grepl("temp", subset_names)] <- "Podoplanin_near_Fol B cells"
names(subset.list) <- subset_names


results.subsets <- crawdad::findTrends(cells,
                                       dist = 100,
                                       shuffle.list = shuffle.list,
                                       subset.list = subset.list,
                                       ncores = ncores,
                                       verbose = TRUE,
                                       returnMeans = FALSE)
## 8.0865 hours to run
saveRDS(results.subsets, file="rafael_analysis/paper/supp5_ksfb_podo_results500.subsets.RDS")
results.subsets <- readRDS("rafael_analysis/paper/supp5_ksfb_podo_results500.subsets.RDS")



## subsets
dats <- crawdad::meltResultsList(results.subsets, withPerms = TRUE)

## Multiple-test correction
ntestss <- length(unique(dats$reference)) * length(unique(dats$neighbor))
psigs <- 0.05/ntestss
zsigs <- round(qnorm(psigs/2, lower.tail = F), 2)



## Podoplanin near follicle B cells
d1 <- dats[grepl(pattern = "Podoplanin_near_Fol B cells", 
                 dats$reference) & dats$neighbor %in% c("CD4 Memory T cells"),]
plt <- vizTrends(dat = d1, facet = FALSE, id = "neighbor", 
                 title = "Podoplanin cells near Fol B cells") +
  ggplot2::scale_x_log10()
plt
dat_filter <- d1 %>% 
  filter(neighbor == 'CD4 Memory T cells')
p <- vizTrends(dat_filter, lines = T, withPerms = T, sig.thresh = zsigs)
p
pdf('rafael_analysis/paper/subsets500/ksfb_podo_near_folB_CD4.pdf')
p 
dev.off()



## Podoplanin not near near follicle B cells
d1 <- dats[grepl(pattern = "Podoplanin_not_near_Fol B cells", 
                 dats$reference) & dats$neighbor %in% c("CD4 Memory T cells"),]
plt <- vizTrends(dat = d1, facet = FALSE, id = "neighbor", 
                 title = "Podoplanin cells not near Fol B cells") +
  ggplot2::scale_x_log10()
plt
dat_filter <- d1 %>% 
  filter(neighbor == 'CD4 Memory T cells')
p <- vizTrends(dat_filter, lines = T, withPerms = T, sig.thresh = zsigs)
p
pdf('rafael_analysis/paper/subsets500/ksfb_podo_not_near_folB_CD4.pdf')
p 
dev.off()
