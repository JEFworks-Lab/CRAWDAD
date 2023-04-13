# Load Xenium data --------------------------------------------------------

corner <- function(x, n=5){
    x[1:n, 1:n]
}

## Code obtained from the xenium_init file.

library(tidyverse)
library(gridExtra)
set.seed(42)

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

####### read in previous analysis results from 10X
## umap
emb <- read.csv('data/Xenium/projection.csv.gz', row.names = 1)
head(emb)

## clustering
com <- read.csv('data/Xenium/clusters.csv.gz', row.names=1)
head(com)
clusters <- com$Cluster
names(clusters) <- rownames(com)
clusters <- factor(clusters)
head(clusters)

## select only clustered cells
cd <- cd[, names(clusters)]

## filter genes
all_genes <- rownames(cd)
good_genes <- all_genes[str_detect(all_genes, 'BLANK_*', negate = TRUE)]
good_genes <- good_genes[str_detect(good_genes, 'NegControlCodeword_*', negate = TRUE)]
good_genes <- good_genes[str_detect(good_genes, 'NegControlProbe_*', negate = TRUE)]
good_genes <- good_genes[str_detect(good_genes, 'antisense_*', negate = TRUE)]
length(good_genes)

## select only good genes
dim(cd)
cd <- cd[which(rownames(cd) %in% good_genes), ]
dim(cd)

## colors
df <- read.csv("data/df_pos_type.csv", row.names = 1)
df$type <- factor(df$type)
df$color <- rainbow(length(levels(df$type)))[df$type]
df$x <- -df$x # for the figures to look like the ones in the paper
df$cid <- rownames(df)
head(df)

library(MERINGUE)

## apply normalization
mtx <- normalizeCounts(counts = cd,
                       log=FALSE,
                       verbose=TRUE)
colnames(mtx) <- colnames(cd)
rownames(mtx) <- rownames(cd)
dim(mtx)
