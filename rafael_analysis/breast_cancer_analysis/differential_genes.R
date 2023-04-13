## Code to identify de genes
## Feb, 09, 2023
## @Rafael


# Load Xenium data --------------------------------------------------------

## Code obtained from the xenium_init file.

library(tidyverse)
library(gridExtra)

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
<- <- sparseMatrix(i = indices[] + 1, p = indptr[],
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

## cluster 10 seems to correspond with cells with low gene count
## should it be removed?
p_low_count <- cellinfo %>%
    ggplot() +
    geom_point(aes(-x_centroid, y_centroid,
                   color = transcript_counts < 25),
               size = .1) +
    scale_color_manual(values = c("lightgray", "red")) +
    labs(color = '< 25 counts')

p_c10 <- cellinfo %>%
    filter(cell_id %in% names(clusters)) %>%
    mutate(cluster = clusters) %>%
    ggplot() +
    geom_point(aes(-x_centroid, y_centroid, color = cluster == 10),
               size = .1) +
    scale_color_manual(values = c("lightgray", "red")) +
    labs(color = 'cluster 10')

grid.arrange(p_low_count, p_c10, nrow = 1)

## most of the cells in cluster 10 have now gene count, but not all
cellinfo %>%
    filter(cell_id %in% names(clusters)) %>%
    mutate(cluster = clusters) %>%
    filter(cluster == 10) %>%
    select(transcript_counts) %>%
    summary()

## filter cells with clusters
length(clusters)

dim(cellinfo)
cellinfo <- cellinfo %>%
    filter(cell_id %in% names(clusters))
dim(cellinfo)

dim(cd)
cd <- cd[, names(clusters)]
dim(cd)



# Differential genes ------------------------------------------------------

## Using Meringue
## https://jef.works/MERINGUE/mOB_analysis

library(MERINGUE)

## should I apply clean counts?
# counts <- cleanCounts(counts = cd,
#                       min.reads = 100,
#                       min.lib.size = 100,
#                       plot=TRUE,
#                       verbose=TRUE)

## apply normalization?
mat <- normalizeCounts(counts = cd,
                       log=FALSE,
                       verbose=TRUE)

# Identify significantly differentially upregulated genes
# in each identified cluster by Wilcox test
dg <- getDifferentialGenes(mat, clusters)
dg.sig <- lapply(dg, function(x) {
    x <- x[x$p.adj < 0.05,]
    x <- na.omit(x)
    x <- x[x$Z > 0, ]
    # x <- x[x$highest,]
    rownames(x)
})
print(lapply(dg.sig, length))

## they got the gene information from https://www.science.org/doi/10.1126/sciadv.abh2169
dg.sig[2]

all_de_genes <- flatten(dg.sig)

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
        mutate(gene = cd[gene, ]) %>%
        ggplot() +
        geom_point(aes(UMAP.1, UMAP.2, color = cluster == selected_cluster),
                   size = .1) +
        scale_color_manual(values = c("lightgray", "blue")) +
        labs(color = selected_cluster)
}

plot_pos_gene <- function(selected_gene){
    cellinfo %>%
        filter(cell_id %in% names(clusters)) %>%
        mutate(cluster = clusters) %>%
        mutate(gene = cd[selected_gene, ]) %>%
        ggplot() +
        geom_point(aes(-x_centroid, y_centroid, color = gene),
                   size = .1) +
        scale_color_continuous(low = "lightgray", high = "black") +
        labs(color = selected_gene)
}

plot_emb_gene <- function(selected_gene){
    emb %>%
        filter(rownames(emb) %in% names(clusters)) %>%
        mutate(cluster = clusters) %>%
        mutate(gene = cd[selected_gene, ]) %>%
        ggplot() +
        geom_point(aes(UMAP.1, UMAP.2, color = gene),
                   size = .1) +
        scale_color_continuous(low = "lightgray", high = "black") +
        labs(color = selected_gene)
}


plot_gene_cluster <- function(gene, marked_cluster){

    ppg <- plot_pos_gene(gene)
    peg <- plot_emb_gene(gene)
    ppc <- plot_pos_cluster(marked_cluster)
    pec <- plot_emb_cluster(marked_cluster)

    grid.arrange(ppc, ppg, pec, peg)
}


# Cluster 14 --------------------------------------------------------------

dg.sig[14]
dg[[14]][dg.sig[[14]], ]

grid.arrange(plot_pos_cluster(14), plot_emb_cluster(14), nrow = 1)

## ANKRD30A: enriched in breast cancer
## https://www.proteinatlas.org/ENSG00000148513-ANKRD30A
## This gene encodes a DNA-binding transcription factor that is uniquely expressed in mammary epithelium and the testis. Altered expression levels have been associated with breast cancer progression.
## https://www.ncbi.nlm.nih.gov/gene/91074
gene <- "ANKRD30A"
plot_gene_cluster(gene, 14)

## Cell type enriched (Adipose visceral - Mesothelial cells, Liver - Cholangiocyte, Skin - Keratinocyte (other), Thyroid - Thyroid glandular cells)
## Low cancer specificity
## Cancer enhanced (Bile duct cancer)
## https://www.proteinatlas.org/ENSG00000188112-C6orf132
gene <- "C6orf132"
plot_gene_cluster(gene, 14)

## CD9
## Cancer-related genes
## ? https://www.proteinatlas.org/ENSG00000010278-CD9 the sources do not seem to agree
gene <- "CD9"
plot_gene_cluster(gene, 14)

## CLDN4: Diseases associated with CLDN4 include Ovarian Cystadenoma and Breast Carcinoma In Situ.
## https://www.genecards.org/cgi-bin/carddisp.pl?gene=CLDN4
gene <- "CLDN4"
plot_gene_cluster(gene, 14)

## Death-associated protein kinase 3 (DAPK3) induces morphological changes in apoptosis when overexpressed in mammalian cells. These results suggest that DAPK3 may play a role in the induction of apoptosis.
## https://www.ncbi.nlm.nih.gov/gene/1613
gene <- "DAPK3"
plot_gene_cluster(gene, 14)

## DNTTIP1: Increases DNTT terminal deoxynucleotidyltransferase activity (in vitro) 1. Also acts as a transcriptional regulator, binding to the consensus sequence 5'-GNTGCATG-3' following an AT-tract. Associates with RAB20 promoter and positively regulates its transcription. Binds DNA and nucleosomes; may recruit HDAC1 complexes to nucleosomes or naked DNA.
## https://www.proteinatlas.org/ENSG00000101457-DNTTIP1
## DNTTIP1 promotes nasopharyngeal carcinoma metastasis via recruiting HDAC1 to DUSP2 promoter and activating ERK signaling pathway
## https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(22)00281-X/fulltext
gene <- "DNTTIP1"
plot_gene_cluster(gene, 14)


## ERBB2: Breast cancer - Unknown function (mainly)
## https://www.proteinatlas.org/ENSG00000141736-ERBB2
## Amplification and/or overexpression of this gene has been reported in numerous cancers, including breast and ovarian tumors.
## https://www.ncbi.nlm.nih.gov/gene/2064
gene <- "ERBB2"
plot_gene_cluster(gene, 14)

## It is also proposed to be a clinically important marker for various types of cancer, particularly those of the breast. However, the role, if any, of GATA3 in the development of these cancers is under study and remains unclear.
## https://journals.lww.com/anatomicpathology/Fulltext/2013/09000/Value_of_GATA3_Immunostaining_in_Tumor_Diagnosis_.6.aspx
gene <- "GATA3"
plot_gene_cluster(gene, 14)

## The FHF complex may function to promote vesicle trafficking and/or fusion via the homotypic vesicular protein sorting complex (the HOPS complex). Contributes to the establishment and maintenance of centrosome function. May function in the positioning or formation of aggresomes, which are pericentriolar accumulations of misfolded proteins, proteasomes and chaperones. FHF complex promotes the distribution of AP-4 complex to the perinuclear area of the cell.
## https://www.proteinatlas.org/ENSG00000095066-HOOK2
gene <- "HOOK2"
plot_gene_cluster(gene, 14)

## Cancer enhanced (Breast cancer, Liver cancer)
## https://www.proteinatlas.org/ENSG00000110169-HPX
gene <- "HPX"
plot_gene_cluster(gene, 14)

## Expression of this gene may be changed in a variety of different cancers and in cardiovascular disease.
## https://www.proteinatlas.org/ENSG00000102554-KLF5
gene <- "KLF5"
plot_gene_cluster(gene, 14)


## ERBB2, GATA3, HPX
fc_clust <- 14
cells_clust <- rownames(cellinfo) %in% names(clusters)[clusters == fc_clust]
fc_genes <- c('ERBB2', 'GATA3', 'HPX')
fc_results <- sapply(fc_genes, function(g){
    a <- median(mat[g, cells_clust])
    b <- median(mat[g, !cells_clust])
    (a/b)
})




# Other marker genes ------------------------------------------------------

## FOXA1: Breast Glandular Cell Marker
gene <- "FOXA1"
which(sapply(dg.sig, function(gl){gene %in% gl}))
marked_cluster <- 19
dg.sig[marked_cluster]
dg[[marked_cluster]][dg.sig[[marked_cluster]], ]
plot_gene_cluster(gene, marked_cluster)


gene <- "CSN2"
which(sapply(dg.sig, function(gl){gene %in% gl}))


gene <- "ESR1"
which(sapply(dg.sig, function(gl){gene %in% gl}))
marked_cluster <- 17
plot_gene_cluster(gene, marked_cluster)



# Co-localization analysis ------------------------------------------------

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

selected_clusters <- c(14, 11, 6)
colors <- c("14" = "blue", "11" = "red", "6" = "yellow", "other" = "lightgray")
grid.arrange(plot_pos_clusters(selected_clusters, colors), plot_emb_clusters(selected_clusters, colors), nrow = 1)

selected_clusters <- c(14, 11, 18)
colors <- c("14" = "blue", "11" = "red", "18" = "yellow", "other" = "lightgray")
grid.arrange(plot_pos_clusters(selected_clusters, colors), plot_emb_clusters(selected_clusters, colors), nrow = 1)

selected_clusters <- c(14, 11, 18)
colors <- c("14" = "blue", "11" = "red", "18" = "yellow", "other" = "lightgray")
grid.arrange(plot_pos_clusters(selected_clusters, colors), plot_emb_clusters(selected_clusters, colors), nrow = 1)

selected_clusters <- c(14, 18)
colors <- c("14" = "blue", "18" = "red", "other" = "lightgray")
grid.arrange(plot_pos_clusters(selected_clusters, colors), plot_emb_clusters(selected_clusters, colors), nrow = 1)

## cluster 11
dg.sig[11]
dg[[11]][dg.sig[[11]], ]

grid.arrange(plot_pos_cluster(11), plot_emb_cluster(11), nrow = 1)

## ACTG2: Group enriched (Breast myoepithelial cells, Smooth muscle cells)
## https://www.proteinatlas.org/ENSG00000163017-ACTG2
## Breast myoepithelial cells
## Their table of annotation, but it does not agree with the explorer
plot_gene_cluster("ACTG2", 11)
plot_gene_cluster("ACTA2", 11)
plot_gene_cluster("KRT15", 11)
plot_gene_cluster("KRT15", 18)

## Squamous epithelial cells - Keratinization (mainly) ## not much info
## https://www.proteinatlas.org/ENSG00000096696-DSP
plot_gene_cluster("DSP", 11)

## cluster 6
grid.arrange(plot_pos_cluster(6), plot_emb_cluster(6), nrow = 1)

dg.sig[6]
head((dg[[6]]))
dg[[6]][dg.sig[[6]], ] %>%
    arrange(desc(Z)) %>%
    head()

## Compare 6 to 14
length(dg.sig[[6]])
length(dg.sig[[14]])
length(intersect(dg.sig[[6]], dg.sig[[14]]))
dg.sig[[6]][!(dg.sig[[6]] %in% dg.sig[[14]])]
## [1] "DMKN"    "ENAH"    "NOSTRIN" "TENT5C"

## DMKN
plot_gene_cluster("DMKN", 6) ## breast glandular cells
plot_gene_cluster("ENAH", 6) ## breast glandular cells
plot_gene_cluster("NOSTRIN", 6) ## adipocytes
plot_gene_cluster("TENT5C", 6) ## T cells

## TACSTD2: This intronless gene encodes a carcinoma-associated antigen. This antigen is a cell surface receptor that transduces calcium signals.
## Cell type enhanced (Basal respiratory cells, Squamous epithelial cells, Suprabasal keratinocytes, Ionocytes, Club cells, Ductal cells, Basal squamous epithelial cells, Basal keratinocytes)

## KRT8

## CEACAM6

## S100A14

## CDH1

## SCD

## see what clusters are colocalized with 6
selected_clusters <- c(5, 2)
colors <- c("5" = "red", "2" = "blue", "other" = "lightgray")
grid.arrange(plot_pos_clusters(selected_clusters, colors), plot_emb_clusters(selected_clusters, colors), nrow = 1)

selected_clusters <- c(5, 2)
colors <- c("5" = "red", "2" = "blue", "other" = "lightgray")
grid.arrange(plot_pos_clusters(selected_clusters, colors), plot_emb_clusters(selected_clusters, colors), nrow = 1)

selected_clusters <- c(12, 11)
colors <- c("12" = "red", "11" = "blue", "other" = "lightgray")
grid.arrange(plot_pos_clusters(selected_clusters, colors), plot_emb_clusters(selected_clusters, colors), nrow = 1)
