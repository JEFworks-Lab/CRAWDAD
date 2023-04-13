## test reading single_cells.csv
## Hickey_Colon_CODEX

data <- read.csv('data/single_cells.csv.gz', row.names=1)
dim(data)
head(data)

## note there are multiple samples
table(data$donor, data$region_num)
table(data$array, data$region_num)
## multiple samples from same patient

vi <- data$array == 'B010A' & data$region_num == 'reg001'

library(ggplot2)
ggplot(data[vi,], aes(x=x, y=y, col=cell_type)) + geom_point(size=0.01)
ggplot(data[vi,], aes(x=x, y=y, col=Neighborhood)) + geom_point(size=0.01)
ggplot(data[vi,], aes(x=x, y=y, col=Community)) + geom_point(size=0.01)
ggplot(data[vi,], aes(x=x, y=y, col=Tissue.Segment)) + geom_point(size=0.01)

## try out crawdad?
## take smaller region
data.sub <- data[vi,]
dim(data.sub)
ggplot(data.sub, aes(x=x, y=y, col=cell_type)) + geom_point(size=0.01) + theme_minimal()
#vii <- data.sub$y > 5000 & data.sub$x > 5000
#data.sub <- data.sub[vii,]
#dim(data.sub)
#ggplot(data.sub, aes(x=x, y=y, col=cell_type)) + geom_point(size=0.01)

## try crawdad
# install.packages('/home/rafael/Desktop/jefworks/projects/multiscale_celltype_colocalization_analysis/crawdad_0.1.0.tar.gz')
library(crawdad)
cells <- crawdad:::toSP(pos = data.sub[,c("x", "y")],
                        celltypes = data.sub[,'cell_type'])
plot(cells, pch=".")

## generate background
library(parallel)
# shuffle.list <- crawdad:::makeShuffledCells(cells,
#                                             resolutions = c(100, 200, 300, 400, 
#                                                             500, 600, 700, 800, 
#                                                             900, 1000, 1250, 
#                                                             1500, 1750, 2000,
#                                                             2250, 2500),
#                                             perms = 1,
#                                             ncores = 16,
#                                             seed = 1,
#                                             verbose = TRUE)
# saveRDS(shuffle.list, 'data/shuffle.list.RDS')
shuffle.list <- readRDS('data/shuffle.list.RDS')

## find trends, passing background as parameter
# results <- crawdad::findTrends(cells,
#                                dist = 100,
#                                shuffle.list = shuffle.list,
#                                ncores = 16,
#                                verbose = TRUE)
# saveRDS(results, 'data/trends.RDS')
results <- readRDS('data/trends.RDS')

## plot
dat <- crawdad::meltResultsList(results, id = "dist_100")
saveRDS(results, 'data/dat.RDS')
grDevices::pdf(file = "results/Hickey_Colon_CODEX_pairwise_trends.pdf", width = 48, height = 48)
crawdad::plotTrends(results = dat, idcol = "id")
dev.off()

sort(table(data.sub$cell_type))

## double check results
ct1 <- 'Enterocyte' 
ct2 <- 'Lymphatic'
foo <- data.sub$cell_type
foo[!(foo %in% c(ct1, ct2))] <- NA
ggplot(data.sub, aes(x=x, y=y, col=foo)) + 
  geom_point(size=0.01) + 
  scale_colour_manual(
    values = c('blue', 'red'),
    na.value = "#eeeeee"
  ) + theme_minimal()

