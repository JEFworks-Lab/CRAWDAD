
# Load package and data ---------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 28

data("pkhl")
names(pkhl)

library(ggplot2)
sort(table(pkhl$celltypes))
set.seed(8) 
cols <- sample(rainbow(lengsth(unique(pkhl$celltypes))))
names(cols) <- unique(pkhl$celltypes)
cols['indistinct'] <- 'grey'
ggplot(pkhl, aes(x=x, y=y, col=celltypes)) + geom_point(size=0.2, alpha=0.5) + theme_void() +
  scale_color_manual(values=cols)

ggplot(pkhl, aes(x=x, y=y, col=celltypes)) + geom_point(size=0.2, alpha=0.5) + theme_void() +
  scale_color_manual(values=cols) + theme(legend.position = 'none')

table(pkhl$celltypes)
vi <- pkhl$celltypes %in% c('Podoplanin', 'CD4 Memory T cells', 'Fol B cells')
ggplot(pkhl, aes(x=x, y=y, col=celltypes)) + geom_point(size=0.2, alpha=0.5) + theme_void() +
  scale_color_manual(values=cols[c('Podoplanin', 'CD4 Memory T cells', 'Fol B cells')],
                     na.value = 'lightgrey') + theme(legend.position = 'none')



# CRAWDAD -----------------------------------------------------------------

## read in Brendan's old results
# results <- readRDS("/Users/jeanfan/Downloads/pkhl.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds")

## convert to SP
cells <- crawdad::toSP(pos = pkhl[,c("x", "y")],
                       celltypes = pkhl$celltypes)

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

save(results, shuffle.list, file="rafael_analysis/paper/fig3.RData")
load("rafael_analysis/paper/fig3.RData")

dat <- crawdad::meltResultsList(results, withPerms = TRUE)
head(dat)

vizAllClusters(cells, pkhl$celltypes)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

# Figure 1d ---------------------------------------------------------------

# ct_order <- rev(c('Podoplanin', 'CD4 Memory T cells', 'Fol B cells', 
#                   'Macrophages', 'CD8 Memory T cells', 'Ki67 proliferating', 
#                   'Myeloid cells', 'B cells, red pulp', 'Blood endothelial',
#                   'Sinusoidal cells', 'Neutrophils/Monocytes'))
# p <- vizColocDotplot(dat) +
#   scale_x_discrete(limits = ct_order) +
#   scale_y_discrete(limits = ct_order, position = 'right') +
#   theme(legend.position='bottom',
#         axis.text.x = element_text(angle = 90, h = 1),
#         legend.box = 'vertical')
p <- vizColocDotplot(dat, s) +
  scale_y_discrete(position = 'right') +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 90, h = 1),
        legend.box = 'vertical')
p
spdf('rafael_analysis/paper/3d.pdf', height = 7, width = 5.5)
p 
dev.off()




# Subsetting --------------------------------------------------------------

binomMat <- crawdad::binomialTestMatrix(cells,
                                        neigh.dist = 100,
                                        ncores = ncores,
                                        verbose = TRUE)

head(binomMat)

saveRDS(binomMat, file="rafael_analysis/paper/fig3_binomMat.RDS")
binomMat <- readRDS("rafael_analysis/paper/fig3_binomMat.RDS")

subset.list <- crawdad::selectSubsets(binomMat,
                                      cells$celltypes,
                                      sub.type = "near",
                                      sub.thresh = 0.05,
                                      ncores = ncores,
                                      verbose = TRUE)

saveRDS(subset.list, file="rafael_analysis/paper/fig3_subset.RDS")
subset.list <- readRDS("rafael_analysis/paper/fig3_subset.RDS")

results.subsets <- crawdad::findTrends(cells,
                                       dist = 100,
                                       shuffle.list = shuffle.list,
                                       subset.list = subset.list,
                                       ncores = ncores,
                                       verbose = TRUE,
                                       returnMeans = FALSE)
## 8.0865 hours to run
results.subsets
saveRDS(results.subsets, file="rafael_analysis/paper/fig3_results.subsets.RDS")
results.subsets <- readRDS("rafael_analysis/paper/fig3_results.subsets.RDS")




# Viz subsetting ----------------------------------------------------------

## cts of interest
plt <- crawdad::vizAllClusters(cells = pkhl,
                               coms = as.factor(pkhl$celltypes),
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
annots_temp <- crawdad::selectLabels(df = pkhl,
                                     com = pkhl$celltypes,
                                     subset_list = subset.list,
                                     cellIDs = c("Fol B cells", "Podoplanin"),
                                     subsetIDs = c("CD4 Memory T cells_near_Fol B cells"))

## visualize the subset only
plt <- crawdad::vizAllClusters(cells = pkhl,
                               coms = annots_temp,
                               ofInterest = NULL,
                               title = "PKHL",
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
pdf('rafael_analysis/paper/3f_trend.pdf')
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
pdf('rafael_analysis/paper/3g_trend.pdf')
p 
dev.off()



# Other -------------------------------------------------------------------

## order
mat <- reshape2::dcast(data = dat, 
                       formula = reference~neighbor,
                       fun.aggregate = sum,
                       value.var = "Z")
rownames(mat) <- mat[,1]
mat <- mat[,-1]
mat <- as.matrix(mat)
quantile(mat)
head(mat)

vi <- which(rownames(mat) == 'indistinct' )
mat <- mat[-vi, -vi]

#mat[is.na(mat)] <- 0
hc <- hclust(dist(t(mat)), method='ward.D')
heatmap(mat, scale='none', 
        Rowv = as.dendrogram(hc),
        Colv = as.dendrogram(hc))

## rafael's code for visualizing
library(tidyverse)
sig_dat <- dat %>%
  filter(abs(Z) >= 1.96) %>% 
  group_by(neighbor, reference) %>% 
  filter(resolution == min(resolution, na.rm = TRUE))
head(sig_dat)

sig_dat$Z[sig_dat$Z > 5] <- 5
sig_dat$Z[sig_dat$Z < -5] <- -5

sig_dat %>% 
  ggplot() +
  geom_point(aes(x=reference, y=neighbor, 
                 color=Z, size=rank(1/resolution))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_gradient2(
    low = 'blue',
    mid = 'lightgrey',
    high = 'red',
    na.value = "#eeeeee"
  ) + 
  scale_size_continuous(range = c(1, 5)) + 
  scale_x_discrete(position = "top", limits=rev(hc$labels[hc$order]))  + 
  scale_y_discrete(limits=rev(hc$labels[hc$order]))  + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) 

vizTrends(dat[dat$neighbor == 'CD4 Memory T cells',], points=FALSE, id='reference', colors=cols) + 
  ggplot2::facet_grid() 

vizTrends(dat[dat$reference == 'CD4 Memory T cells',], points=FALSE, id='neighbor', colors=cols) + 
  ggplot2::facet_grid() 

vizTrends(dat[dat$neighbor == 'CD4 Memory T cells' & dat$reference %in% c('Fol B cells',  'Podoplanin'),], points=FALSE, id='neighbor') + 
  ggplot2::facet_grid(reference~neighbor) 

vizTrends(dat[dat$reference == 'CD4 Memory T cells' & dat$neighbor %in% c('Fol B cells',  'Podoplanin'),], points=FALSE, id='neighbor') + 
  ggplot2::facet_grid(neighbor~reference) 

### triplet results
subset <- readRDS("/Users/jeanfan/Downloads/pkhl.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds")

sub <- subset[['100']][names(subset[['100']])[grepl('CD4 Memory T cells_near', names(subset[['100']]))]]
dat <- crawdad::meltResultsList(sub)
head(dat)

mat <- reshape2::dcast(data = dat,formula = reference~neighbor,fun.aggregate = sum,value.var = "Z")
rownames(mat) <- mat[,1]
mat <- mat[,-1]
mat <- as.matrix(mat)
quantile(mat)
head(mat)

dat <- dat[dat$neighbor %in% c('Fol B cells', 'CD4 Memory T cells','Podoplanin'),]

## rafael's code for visualizing
library(tidyverse)
sig_dat <- dat %>%
  filter(abs(Z) >= 1.96) %>% 
  group_by(neighbor, reference) %>% 
  filter(resolution == min(resolution, na.rm = TRUE))
head(sig_dat)

sig_dat$Z[sig_dat$Z > 5] <- 5
sig_dat$Z[sig_dat$Z < -5] <- -5

sig_dat %>% 
  ggplot() +
  geom_point(aes(x=reference, y=neighbor, 
                 color=Z, size=rank(1/resolution))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_gradient2(
    low = 'blue',
    mid = 'lightgrey',
    high = 'red',
    na.value = "#eeeeee"
  ) + 
  scale_size_continuous(range = c(1, 5)) + 
  #scale_x_discrete(position = "top", limits=colnames(mat)[hc$order])  + 
  #scale_y_discrete(limits=hc$labels[hc$order])  + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) 


##### visualize subsets
subcells <- readRDS('~/Downloads/pkhl.subsets.near.subdist100.rds')
names(subcells)
names(subcells)[grepl('CD4 Memory T cells_near', names(subcells))]

a <- subcells[["CD4 Memory T cells_near_Fol B cells, ctr"]]
b <- subcells[["CD4 Memory T cells_near_Fol B cells, out"]]
c <- subcells[["CD4 Memory T cells_near_B cells, red pulp"]]
#other <- setdiff(rownames(pkhl)[pkhl$celltypes == 'CD4 Memory T cells'], c(a,b,c))

foo <- pkhl
foo[a,]$celltypes <- 'CD4 Memory T cells_near_Fol B cells'
foo[b,]$celltypes <- 'CD4 Memory T cells_near_Fol B cells'
foo[c,]$celltypes <- "CD4 Memory T cells_near_B cells, red pulp"
#foo[other,]$celltypes <- 'CD4 Memory T cells_notnear_Fol B cells'
vi <- foo$celltypes %in% c('Podoplanin', 'CD4 Memory T cells_near_Fol B cells', 'CD4 Memory T cells_near_B cells, red pulp')
colsfoo <- c("#FFFF00", "#FF0080", "#0000FF")
names(colsfoo) <-  c('Podoplanin', 'CD4 Memory T cells_near_Fol B cells', "CD4 Memory T cells_near_B cells, red pulp")
ggplot(foo, aes(x=x, y=y, col=celltypes)) + geom_point(size=0.2, alpha=0.5) + theme_void() + 
  scale_color_manual(values=colsfoo,
                     na.value = 'lightgrey') 
## barplot
freq <- sapply(subcells[names(subcells)[grepl('CD4 Memory T cells_near', names(subcells))]], length)
freq['CD4 Memory T cells_near_Fol B cells'] <- freq['CD4 Memory T cells_near_Fol B cells, ctr'] + freq['CD4 Memory T cells_near_Fol B cells, out']
freq <- freq[!(names(freq) %in% c('CD4 Memory T cells_near_Fol B cells, ctr', 'CD4 Memory T cells_near_Fol B cells, out'))]
freq <- freq[sort(names(freq))]
barplot(freq, las=2)

## avoid mutual exclusivity issue
foo <- pkhl
a <- subcells[["CD4 Memory T cells_near_Fol B cells, ctr"]]
b <- subcells[["CD4 Memory T cells_near_Fol B cells, out"]]
foo[a,]$celltypes <- 'CD4 Memory T cells_near_Fol B cells'
foo[b,]$celltypes <- 'CD4 Memory T cells_near_Fol B cells'
table(foo$celltypes)
vi <- foo$celltypes %in% c('Podoplanin', 'CD4 Memory T cells_near_Fol B cells', 'Fol B cells')
ggplot(foo, aes(x=x, y=y, col=celltypes)) + geom_point(size=0.2, alpha=0.5) + theme_void() +
  scale_color_manual(values=cols[c('Podoplanin', 'CD4 Memory T cells', 'Fol B cells')],
                     na.value = 'lightgrey') + theme(legend.position = 'none')



