## Code to compare clusters 14 and 6
## Feb, 21, 2023
## @Rafael


# Load Xenium data --------------------------------------------------------

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

## select only clustered cells
cd <- cd[, names(clusters)]

## colors
df <- read.csv("data/df_pos_type.csv", row.names = 1)
df$type <- factor(df$type)
df$color <- rainbow(length(levels(df$type)))[df$type]
df$x <- -df$x # for the figures to look like the ones in the paper
head(df)



# Cluster analysis --------------------------------------------------------

## all clusters
p_emb_clus <- emb %>%
    filter(rownames(emb) %in% names(clusters)) %>%
    mutate(cluster = clusters) %>%
    ggplot() +
    geom_point(aes(UMAP.1, UMAP.2, color = cluster),
               size = .1) +
    scale_colour_manual(values=setNames(df$color, df$type)) +
    guides(colour = guide_legend(override.aes = list(size=5)))

p_pos_clus <- cellinfo %>%
    filter(cell_id %in% names(clusters)) %>%
    mutate(cluster = clusters) %>%
    ggplot() +
    geom_point(aes(-x_centroid, y_centroid, color = cluster),
               size = .1) +
    scale_colour_manual(values=setNames(df$color, df$type)) +
    guides(colour = guide_legend(override.aes = list(size=5)))

grid.arrange(p_emb_clus, p_pos_clus, nrow = 1)


# Differential genes ------------------------------------------------------

library(MERINGUE)

## apply normalization
mtx <- normalizeCounts(counts = cd,
                       log=FALSE,
                       verbose=TRUE)
colnames(mtx) <- colnames(cd)
rownames(mtx) <- rownames(cd)

## identify significantly differentially upregulated genes
dg <- getDifferentialGenes(mtx, clusters)
dg.sig <- lapply(dg, function(x) {
    x <- x[x$p.adj < 0.05,]
    x <- na.omit(x)
    x <- x[x$Z > 0, ]
    # x <- x[x$highest,]
    rownames(x)
})
print(lapply(dg.sig, length))

## highest gene exp across populations
dgh.sig <- lapply(dg, function(x) {
    x <- x[x$p.adj < 0.05,]
    x <- na.omit(x)
    x <- x[x$Z > 0, ]
    x <- x[x$highest,]
    rownames(x)
})
print(lapply(dgh.sig, length))

## heatmap
# dg.genes <- unlist(dg.sig)
# ggroup <- unlist(lapply(1:length(dg.sig), function(i) {
#     rep(names(dg.sig)[i], length(dg.sig[[i]]))
# }))
# names(ggroup) <- dg.genes
# ggroup <- factor(ggroup)
#
# ccol <- rainbow(length(levels(clusters)))[clusters]
# names(ccol) <- names(clusters) # column colors
# gcol <- rainbow(length(levels(ggroup)), v=0.5)[ggroup]
# names(gcol) <- names(ggroup) # row colors
#
# m <- as.matrix(mtx[dg.genes, names(sort(clusters))])
# m <- winsorize(t(scale(t(m))))
# heatmap(m, scale="none",
#         Colv=NA, Rowv=NA, labRow=NA, labCol=NA,
#         ColSideColors=ccol[colnames(m)],
#         RowSideColors=gcol[rownames(m)],
#         col=colorRampPalette(c('blue', 'white', 'red'))(100)
# )



# Plotting functions ------------------------------------------------------

plot_pos_cluster <- function(selected_cluster){
    cellinfo %>%
        filter(cell_id %in% names(clusters)) %>%
        mutate(cluster = clusters) %>%
        ggplot() +
        geom_point(aes(-x_centroid, y_centroid,
                       color = cluster == selected_cluster),
                   size = .1) +
        scale_color_manual(values = c("lightgray", "blue")) +
        labs(color = selected_cluster)
}

plot_emb_cluster <- function(selected_cluster){
    emb %>%
        filter(rownames(emb) %in% names(clusters)) %>%
        mutate(cluster = clusters) %>%
        ggplot() +
        geom_point(aes(UMAP.1, UMAP.2, color = cluster == selected_cluster),
                   size = .1) +
        scale_color_manual(values = c("lightgray", "blue")) +
        labs(color = selected_cluster)
}

plot_pos_gene <- function(selected_gene, log10_norm){
    cellinfo %>%
        filter(cell_id %in% names(clusters)) %>%
        mutate(cluster = clusters) %>%
        mutate(gene = if (log10_norm) log10(cd[selected_gene, ] + 1) else cd[selected_gene, ]) %>%
        ggplot() +
        geom_point(aes(-x_centroid, y_centroid, color = gene),
                   size = .1) +
        scale_color_continuous(low = "lightgray", high = "black") +
        labs(color = selected_gene)
}

plot_emb_gene <- function(selected_gene, log10_norm){
    emb %>%
        filter(rownames(emb) %in% names(clusters)) %>%
        mutate(cluster = clusters) %>%
        mutate(gene = if (log10_norm) log10(cd[selected_gene, ] + 1) else cd[selected_gene, ]) %>%
        ggplot() +
        geom_point(aes(UMAP.1, UMAP.2, color = gene),
                   size = .1) +
        scale_color_continuous(low = "lightgray", high = "black") +
        labs(color = selected_gene)
}

plot_gene_cluster <- function(gene, marked_cluster, log10_norm = TRUE){

    ppg <- plot_pos_gene(gene, log10_norm)
    peg <- plot_emb_gene(gene, log10_norm)
    ppc <- plot_pos_cluster(marked_cluster)
    pec <- plot_emb_cluster(marked_cluster)

    grid.arrange(ppc, ppg, pec, peg)
}

## plot co-localized clusters in pos and emb
plot_pos_clusters <- function(selected_clusters, colors){
    cellinfo %>%
        filter(cell_id %in% names(clusters)) %>%
        mutate(cluster = sapply(clusters, function(c){if (c %in% selected_clusters) c else "other"}) ) %>%
        ggplot() +
        geom_point(aes(-x_centroid, y_centroid,
                       color = cluster),
                   size = .1) +
        scale_color_manual(values = colors) +
        guides(color = guide_legend(override.aes = list(size = 10)))
}

plot_emb_clusters <- function(selected_clusters, colors){
    emb %>%
        filter(rownames(emb) %in% names(clusters)) %>%
        mutate(cluster = sapply(clusters, function(c){if (c %in% selected_clusters) c else "other"}) ) %>%
        ggplot() +
        geom_point(aes(UMAP.1, UMAP.2, color = cluster),
                   size = .1) +
        scale_color_manual(values = colors) +
        guides(color = guide_legend(override.aes = list(size = 10)))
}



# Comparing clusters 14 and 6 ---------------------------------------------

colors <- c("14" = "blue", "6" = "red", "other" = "lightgray")
grid.arrange(plot_pos_clusters(names(colors), colors),
             plot_emb_clusters(names(colors), colors),
             nrow = 1)

## dg with highest expression in cluster
dgh.sig[[6]] ## empty
dgh.sig[[14]]
# [1] "ANKRD30A" "C6orf132" "CD9"      "CLDN4"    "DAPK3"    "DNTTIP1"  "GATA3"    "HOOK2"    "HPX"
# [10] "KLF5"     "KLRF1"    "KRT8"     "LYPD3"    "MLPH"     "SCD"      "SDC4"     "SEC24A"   "SMS"
# [19] "TFAP2A"   "TRAPPC3"  "VOPP1"

## get mean expression across different clusters
## the others column represents the all clusters expect the selected_clusters
mean_exp <- function(selected_clusters, selected_genes){

    other_clusters <- unique(clusters)[!(unique(clusters) %in% selected_clusters)]
    list_clusters <- append(as.list(selected_clusters), list(other_clusters))
    names(list_clusters) <- append(selected_clusters, "others")

    list_exps <- list()
    for (name_cg in names(list_clusters)) {
        cells_cluster <- names(clusters[clusters %in% list_clusters[[name_cg]]])
        mgexp <- sapply(selected_genes, function(g){
            mean(mtx[g, cells_cluster])
        })
        list_exps[[name_cg]] <- mgexp
    }

    return(do.call(cbind, list_exps))
}

## get mean count across different clusters
## the others column represents the all clusters expect the selected_clusters
mean_count <- function(selected_clusters, selected_genes){

    other_clusters <- unique(clusters)[!(unique(clusters) %in% selected_clusters)]
    list_clusters <- append(as.list(selected_clusters), list(other_clusters))
    names(list_clusters) <- append(selected_clusters, "others")

    list_exps <- list()
    for (name_cg in names(list_clusters)) {
        cells_cluster <- names(clusters[clusters %in% list_clusters[[name_cg]]])
        mgexp <- sapply(selected_genes, function(g){
            mean(mtx[g, cells_cluster] > 0)
        })
        list_exps[[name_cg]] <- mgexp
    }

    return(do.call(cbind, list_exps))
}

## did not work
mean_exp2 <- function(selected_clusters, selected_genes){

    other_clusters <- unique(clusters)[!(unique(clusters) %in% selected_clusters)]
    list_clusters <- append(as.list(selected_clusters), list(other_clusters))
    names(list_clusters) <- append(selected_clusters, "others")

    comb <- expand.grid(list_clusters, selected_genes)
    purrr::walk2(comb$Var1, comb$Var2, function(g, c){
        mean(mtx[g, c])
    })
}

## heatmap of gene expression
mean_exp(c(14, 6), dgh.sig[[14]]) %>% heatmap()
## ? aren't all clusters supposed to be highest in cluster 14

## fancy heatmap
mean_exp(c(14, 6), dgh.sig[[14]]) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene") %>%
    pivot_longer(-gene, names_to = "cluster", values_to = "mean") %>%
    ggplot(aes(x = cluster, y = gene, fill = mean)) +
    geom_tile() +
    scale_fill_gradient(low="lightgray", high="red")

## heatmap for 14
mean_exp(c(14, 6), dg.sig[[14]]) %>% heatmap()

## heatmap for 6
mean_exp(c(14, 6), dg.sig[[6]]) %>% heatmap()

## Compare 6 to 14
length(dg.sig[[6]])
length(dg.sig[[14]])
length(intersect(dg.sig[[6]], dg.sig[[14]]))
dg.sig[[6]][!(dg.sig[[6]] %in% dg.sig[[14]])]
## [1] "DMKN"    "ENAH"    "NOSTRIN" "TENT5C"

## DMKN
plot_gene_cluster("DMKN", 6) ## from gene panel: breast glandular cells
plot_gene_cluster("ENAH", 6) ## from gene panel: breast glandular cells
plot_gene_cluster("NOSTRIN", 6) ## from gene panel: adipocytes
plot_gene_cluster("TENT5C", 6) ## from gene panel: T cells

## Comparing 14 and 6
cells <- which(clusters %in% c(14, 6))
dg2 <- getDifferentialGenes(mtx[, cells], clusters[cells])
dg2.sig <- lapply(dg2, function(x) {
    x <- x[x$p.adj < 0.05,]
    x <- na.omit(x)
    x <- x[x$Z > 0, ]
    # x <- x[x$highest,]
    rownames(x)
})
print(lapply(dg2.sig, length))

## highest gene exp across populations
dgh2.sig <- lapply(dg2, function(x) {
    x <- x[x$p.adj < 0.05,]
    x <- na.omit(x)
    x <- x[x$Z > 0, ]
    x <- x[x$highest,]
    rownames(x)
})
print(lapply(dgh2.sig, length))

## Compare 6 to 14
length(dgh2.sig[[6]])
length(dgh2.sig[[14]])
length(intersect(dgh2.sig[[6]], dgh2.sig[[14]]))
# [1] 0

mean_exp(c(14, 6), dgh2.sig[[14]]) %>% heatmap()
mean_exp(c(14, 6), dgh2.sig[[6]]) %>% heatmap()

plot_gene_cluster("GATA3", 14) ## from gene panel:

plot_gene_cluster("AGR3", 6) ## from gene panel: endoplasmic reticulum ?
plot_gene_cluster("BACE2", 6) ## from gene panel:
plot_gene_cluster("DMKN", 6) ## from gene panel: breast glandular cells

setdiff(unique(clusters), c(10)) %>% mean_exp(dg.sig[[6]]) %>% heatmap()

## mean exp
c(6, 14) %>% mean_exp(c('GATA3', 'CEACAM6', 'TACSTD2')) %>% heatmap(main="mean exp")
## mean count
## this code does not work because I am checking the percentage of cells
## expressing the gene at all cluster, not in the region, as mentioned in the paper
c(6, 14) %>% mean_count(c('MZB1', 'AGR3', 'MKI67')) %>% heatmap(main="cell fraction")

## comparing this last visualization with the mean count
## 14 is DCIS #1 and 6 is DCIS #2
## although the computation is not regorous, the results agree with the paper
## https://www.biorxiv.org/content/biorxiv/early/2022/10/07/2022.10.06.510405/F6.large.jpg
colors <- c("14" = "blue", "6" = "red", "other" = "lightgray")
grid.arrange(plot_pos_clusters(names(colors), colors),
             plot_emb_clusters(names(colors), colors),
             nrow = 1)



# DCIS #2 -----------------------------------------------------------------

grid.arrange(plot_pos_cluster(6),
             plot_pos_cluster(11),
             plot_pos_cluster(12),
             nrow = 2)

## significant genes
dgh.sig[[11]]
dgh.sig[[12]]
mean_exp(c(6, 11, 12), c(dgh.sig[[11]], dgh.sig[[12]])) %>% heatmap()
mean_exp(c(6, 11, 12), dg.sig[[11]]) %>% heatmap()

dg.sig[[11]]
dg.sig[[12]]

plot_gene_cluster("ERBB2", 6, log10_norm = FALSE)
plot_gene_cluster("GATA3", 6, log10_norm = FALSE)
plot_gene_cluster("TACSTD2", 6, log10_norm = FALSE)

colors <- c("6" = "blue", "11" = "red", "12" = "yellow", "other" = "lightgray")
grid.arrange(plot_pos_clusters(names(colors), colors),
             plot_emb_clusters(names(colors), colors),
             nrow = 2)

## reduced myoepithelial markers KRT15, KRT23, and ALDH1A3 which could potentially be associated with increased invasiveness suggested by the higher expression of invasive markers found in DCIS #2
setdiff(unique(clusters), c(10)) %>% mean_exp(c('KRT15', 'KRT23', 'ALDH1A3')) %>% heatmap()
## these genes are highly expressed in cluster 18, which is myoepithelial ?


# DCIS #1 -----------------------------------------------------------------

colors <- c("14" = "blue", "17" = "red", "18" = "yellow", "22" = "purple", "other" = "lightgray")
grid.arrange(plot_pos_clusters(names(colors), colors),
             plot_emb_clusters(names(colors), colors),
             nrow = 1)

## significant genes
dgh.sig[[14]]
dgh.sig[[17]]
dgh.sig[[18]]
dgh.sig[[22]]
mean_exp(c(14, 17, 18, 22),
         c(dgh.sig[[14]], dgh.sig[[17]], dgh.sig[[18]], dgh.sig[[22]])) %>%
    heatmap()

## Comparing 14 and 17
cells <- which(clusters %in% c(14, 17))
dg_14_17 <- getDifferentialGenes(mtx[, cells], clusters[cells])
dg_14_17.sig <- lapply(dg_14_17, function(x) {
    x <- x[x$p.adj < 0.05,]
    x <- na.omit(x)
    x <- x[x$Z > 0, ]
    # x <- x[x$highest,]
    rownames(x)
})
print(lapply(dg_14_17.sig, length))

## highest gene exp across populations
dgh_14_17.sig <- lapply(dg_14_17, function(x) {
    x <- x[x$p.adj < 0.05,]
    x <- na.omit(x)
    x <- x[x$Z > 0, ]
    x <- x[x$highest,]
    rownames(x)
})
print(lapply(dgh_14_17.sig, length))

## Compare 17 to 14
length(dgh_14_17.sig[[17]])
length(dgh_14_17.sig[[14]])
length(intersect(dgh_14_17.sig[[17]], dgh_14_17.sig[[14]]))
# [1] 0

mean_exp(c(14, 17), dgh_14_17.sig[[14]]) %>% heatmap()
mean_exp(c(14, 17), dgh_14_17.sig[[17]]) %>% heatmap()



# Comparison 11 and 18 ----------------------------------------------------

colors <- c("11" = "blue", "18" = "red", "other" = "lightgray")
grid.arrange(plot_pos_clusters(names(colors), colors),
             plot_emb_clusters(names(colors), colors),
             nrow = 1)

## significant genes
dgh.sig[[11]]
dgh.sig[[18]]
mean_exp(c(11, 18),
         c(dgh.sig[[11]], dgh.sig[[18]])) %>%
    heatmap()

## Comparing 14 and 17
cells <- which(clusters %in% c(14, 17))
dg_14_17 <- getDifferentialGenes(mtx[, cells], clusters[cells])
dg_14_17.sig <- lapply(dg_14_17, function(x) {
    x <- x[x$p.adj < 0.05,]
    x <- na.omit(x)
    x <- x[x$Z > 0, ]
    # x <- x[x$highest,]
    rownames(x)
})
print(lapply(dg_14_17.sig, length))

## highest gene exp across populations
dgh_14_17.sig <- lapply(dg_14_17, function(x) {
    x <- x[x$p.adj < 0.05,]
    x <- na.omit(x)
    x <- x[x$Z > 0, ]
    x <- x[x$highest,]
    rownames(x)
})
print(lapply(dgh_14_17.sig, length))

## Compare 17 to 14
length(dgh_14_17.sig[[17]])
length(dgh_14_17.sig[[14]])
length(intersect(dgh_14_17.sig[[17]], dgh_14_17.sig[[14]]))
# [1] 0

mean_exp(c(14, 17), dgh_14_17.sig[[14]]) %>% heatmap()
mean_exp(c(14, 17), dgh_14_17.sig[[17]]) %>% heatmap()
