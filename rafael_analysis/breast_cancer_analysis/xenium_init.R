## Code to help you get started with the Xenium data
## Oct 27, 2022
## @Jean

library(tidyverse)

########## read in the data
## read in cell metainfo
cellinfo <- read.csv('data/Xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.csv.gz')
head(cellinfo)
pos <- cbind(cellinfo$x_centroid, cellinfo$y_centroid)
rownames(pos) <- cellinfo$cell_id


cellinfo %>% ggplot() +
    geom_point(aes(x=x_centroid, y=y_centroid),
               size=.1)

## read in cell gene expression counts
library(rhdf5)
file <- 'data/Xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5'
rhdf5::h5ls(file)
data <- rhdf5::h5read(file, 'matrix')
h5closeAll()
names(data)

counts <- data$data
indices <- data$indices
indptr <- data$indptr
shp <- data$shape
features <- data$features
barcodes <- data$barcodes
library(Matrix)
cd <- sparseMatrix(i = indices[] + 1, p = indptr[],
                   x = as.numeric(x = counts[]), dims = shp[], repr = "T")
rownames(cd) <- features$name
colnames(cd) <- barcodes
class(cd)
cd <- Matrix(cd, sparse=TRUE)

libsize <- colSums(cd)
col <- colorRampPalette(c('white', 'red'))(100)[round(libsize/max(libsize)*100)]
plot(pos, col=col, pch=".")

####### read in previous analysis results from 10X
## umap
emb <- read.csv('data/Xenium/projection.csv.gz', row.names = 1)
head(emb)
plot(emb, col=col, pch=".", main='UMAP')

## clustering
com <- read.csv('data/Xenium/clusters.csv.gz', row.names=1)
head(com)
clusters <- com$Cluster
names(clusters) <- rownames(com)
clusters <- factor(clusters)
head(clusters)

## map from factors to color hues
groupcol <- rainbow(length(levels(clusters)))[clusters]
names(groupcol) <- names(clusters)
par(mfrow=c(1,2))
plot(pos, col=groupcol[rownames(pos)], pch=".")
plot(emb, col=groupcol[rownames(emb)], pch=".", main='UMAP')


# Saving pos and clusters -------------------------------------------------

df_pos_type.csv <- cbind(pos, clusters[rownames(pos)], groupcol[rownames(pos)])
colnames(df_pos_type.csv) <- c("x", "y", "type", "color")
head(df_pos_type.csv)
write.csv(df_pos_type.csv, "data/df_pos_type.csv")

####### exploring the data
## log10(counts per million normalization + 1)
mat <- log10(Matrix::t(Matrix::t(cd)/Matrix::colSums(cd)) * 1e6 + 1)

g <- 'MS4A1' ## CD20 a B cell marker
gexp <- mat[g,]
hist(gexp)

## helper function to map from numeric to color shades
map2col <- function (x, pal = colorRampPalette(c("lightgrey", "red"))(100),
          na.col = "lightgrey", limits = NULL) {
  original <- x
  x <- na.omit(x)
  if (is.null(limits))
    limits = range(x)
  y <- pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) +
                                 1), all.inside = TRUE)]
  names(y) <- names(x)
  colors <- rep(na.col, length(original))
  names(colors) <- names(original)
  colors[names(y)] <- y
  return(colors)
}

gexpcol <- map2col(gexp)
plot(pos, col=gexpcol[rownames(pos)], pch=".", main=g)
plot(emb, col=gexpcol[rownames(emb)], pch=".", main='UMAP')

####### can also use packages
library(MERINGUE)
par(mfrow=c(1,2))
MERINGUE::plotEmbedding(pos, groups=clusters, cex=0.1)
MERINGUE::plotEmbedding(emb, groups=clusters, cex=0.1, mark.clusters = TRUE, mark.cluster.cex = 1)

pdf("clusters.pdf")
dev.off()

####### questions to work through
## 1. Which pairs of cell clusters are spatially co-localized?
## 2. For pairs of cell clusters that are spatially co-localized, are there specific subsets
## based on co-localization to a cells of a third cell cluster that is distinct?
## 3. Based on our knowledge of marker genes or other data-driven approaches, annotate these cell clusters
## 4. Based on these cell cluster annotations, is there any primary literature
## supporting the biological relevance of your identified spatial co-localization patterns?
