install.packages("/home/rafael/Desktop/jefworks/projects/multiscale_celltype_colocalization_analysis/msctcla_0.1.0.tar.gz", repos = NULL, type="source", build_vignettes=FALSE)

library(msctcla)
library(tidyverse)

# Set up parameters -------------------------------------------------------

## full path to the subset dataset list rds
subsetDataPath <- "data/subsets.near.rds"

## test if a cell is "near" or "away" a given cell type
## by testing if it's neighbors are significantly enriched or depleted in a cell type via binomial test
subsetType <- "near"

## distance to define subsets; ie the distance to define neighbors for this test
subdist <- 100


## for "away", loosely test for enrichment, then take the cells that still fail this test
## so assume that they must be highly depleted
if(subsetType == "near"){
    subThresh <- 0.05
} else if(subsetType == "away"){
    subThresh <- 0.5
} else (
    stop("subsetType must be either 'near' or 'away'")
)

## number of cores
ncs <- 16

## the `pos` data.frame, but these lines load it back in properly
meta <- read.csv2(file = "data/df_pos_type.csv", row.names = 1, sep=',')

## filter for clusters
head(meta)
meta <- meta[which(meta$type %in% c(14, 18)), ]
head(meta)

meta <- meta[,c("x", "y", "type")]
## make sure the coordinates are numeric
meta <- meta %>%
    dplyr::mutate_at(vars(x, y), as.numeric)


## create the spatial dataframe:
pos <- meta[,c("x", "y")]
celltypes <- meta[,3]

if(length(levels(celltypes)) == 0){
    message("Warning: `celltypes` does not have levels. Creating levels from values")
    celltypes <- factor(celltypes)
    names(celltypes) <- rownames(pos)
}

cells <- sp::SpatialPointsDataFrame(
    coords = as.data.frame(pos),
    data = data.frame(
        celltypes = celltypes
        #name=rownames(pos)
    ))
cells <- sf::st_as_sf(cells)

## Change rowname assignments of cells to integers.
## Solution to keep rows in same order later on when
## randomly shuffling cell labels
rownames(cells) <- as.character(1:dim(cells)[1])
## changed to the actual names
rownames(cells) <- rownames(meta)

# make assumption that cell type attribute is constant throughout the geometries of each cell
## it removed the warning that keep popping up, which says this assumption is made anyways
sf::st_agr(cells) <- "constant"




# get the subsets ---------------------------------------------------------

## if subset data exists, set up parameters to load it in, otherwise set them up to save it once made
if(file.exists(subsetDataPath)){
    loadSubset <- subsetDataPath
    saveSubset <- NA
} else {
    loadSubset <- NA
    saveSubset <- subsetDataPath
}

if(assertthat::is.string(loadSubset)){
    subset.list <- readRDS(file = loadSubset)
} else {

    ## parallel shuffling of the grids in each resolution
    subset.list <- getSubsets(cells = cells,
                              sub.dist = subdist,
                              sub.type = subsetType,
                              sub.thresh = subThresh,
                              ncores = ncs,
                              verbose = TRUE)
    if(assertthat::is.string(saveSubset)){
        saveRDS(object = subset.list, file = saveSubset)
    }
}



# finding trends parameters -----------------------------------------------

## full path to the shuffled data rds
## if it exists, it will be loaded into the function, if not, will be created and saved at this location
shuffledDataPath <- paste0("data/subsets.shuffled.rds")
outfile_subsets_triplet <- paste0("data/subsets.triplet.near.binom.rds")


## full path to the subset dataset list rds
subsetDataPath <- paste0("data/subsets.near.rds")

## test if a cell is "near" or "away" a given cell type
## by testing if it's neighbors are significantly enriched or depleted in a cell type via binomial test
subsetType <- "near"

## distance to define subsets; ie the distance to define neighbors for this test
subdist <- 100


## for "away", loosely test for enrichment, then take the cells that still fail this test
## so assume that they must be highly depleted
if(subsetType == "near"){
    subThresh <- 0.05
} else if(subsetType == "away"){
    subThresh <- 0.5
} else (
    stop("subsetType must be either 'near' or 'away'")
)


## number of cores
ncs <- 16

## number of permutations
prms <- 1

## seed
sd <- 1

## don't count neighbor cells more than once
removedups <- TRUE

## for subset, use one neighbor distance
distances <- c(100)
## shuffling resolutions
resolutions <- c(50, 250, 500, 1000, 2500, 5000)



# find trends -------------------------------------------------------------

## if subset data exists, set up parameters to load it in, otherwise set them up to save it once made
if(file.exists(subsetDataPath)){
    loadSubset <- subsetDataPath
    saveSubset <- NA
} else {
    loadSubset <- NA
    saveSubset <- subsetDataPath
}

triplet_results <- lapply(distances, function(d){

    if(file.exists(shuffledDataPath)){

        results <- findTrends(pos = meta[,c("x", "y")],
                                celltypes = meta[,3],
                                resolutions,
                                dist = d,
                                sub.dist = subdist,
                                sub.type = subsetType,
                                sub.thresh = subThresh,
                                perms = prms,
                                ncores = ncs,
                                verbose = TRUE,
                                loadShuffleFile = shuffledDataPath,
                                loadSubsetFile = loadSubset,
                                saveSubsetFile = saveSubset,
                                seed = sd)
    } else {
        results <- findTrends(pos = meta[,c("x", "y")],
                                celltypes = meta[,3],
                                resolutions,
                                dist = d,
                                sub.dist = subdist,
                                sub.type = subsetType,
                                sub.thresh = subThresh,
                                perms = prms,
                                ncores = ncs,
                                verbose = TRUE,
                                saveShuffleFilePath = shuffledDataPath,
                                loadSubsetFile = loadSubset,
                                saveSubsetFile = saveSubset,
                                seed = sd)
    }
    return(results)

})

names(triplet_results) <- distances

saveRDS(object = triplet_results, file = outfile_subsets_triplet)

results <- triplet_results[["100"]]
df_results <- meltResultsList(resultsList = results)
head(df_results)
write.csv(df_results, "results/df_results_100.csv")

# Visualization -----------------------------------------------------------

library(gridExtra)
plot_pos_clusters <- function(selected_clusters, colors){
    cellinfo %>%
        mutate(cluster = sapply(type, function(c){if (c %in% selected_clusters) c else "other"}) ) %>%
        ggplot() +
        geom_point(aes(-x, y,
                       color = cluster),
                   size = .1) +
        scale_color_manual(values = colors) +
        guides(color = guide_legend(override.aes = list(size = 10)))
}

plot_emb_clusters <- function(selected_clusters, colors){
    emb %>%
        filter(rownames(emb) %in% rownames(df_ac)) %>%
        mutate(cluster = sapply(df_ac[which(rownames(df_ac) %in% rownames(emb)), ]$type,
                                function(c){if (c %in% selected_clusters) c else "other"})) %>%
        ggplot() +
        geom_point(aes(UMAP.1, UMAP.2, color = cluster),
                   size = .1) +
        scale_color_manual(values = colors) +
        guides(color = guide_legend(override.aes = list(size = 10)))
}

pdf(file = "results/tripet/grid.pdf", width = 12, height = 12)
df_results %>% ggplot() +
    geom_line(aes(x=resolution, y=Z)) +
    facet_grid(reference ~ neighbor,
               labeller = labeller(.rows = label_both, .cols = label_both),
               scales="free")
dev.off()

## Visualize Z scores
df_results$Z %>%
    hist()

## Convert to factor to plot the data
df_results$neighbor <- as.factor(df_results$neighbor)
df_results$reference <- as.factor(df_results$reference)

## individual references
library(ggrepel)
ref_plots <- list()
for (ref in unique(df_results$reference)){
    ref_plots[[ref]] <-
        df_results %>%
        filter(reference == ref) %>%
        ggplot(aes(x=resolution, y=Z, color=neighbor)) +
        geom_line() +
        # scale_colour_manual(values=setNames(df$color, df$type)) +
        labs(title = ref) +
        theme(legend.position = "none") +
        geom_text_repel(aes(label = neighbor),
                        nudge_x = 1,
                        na.rm = TRUE,
                        force = 1,
                        box.padding = 1,
                        segment.alpha = .5,
                        data = df_results %>%
                            filter(reference == ref) %>%
                            group_by(neighbor) %>%
                            filter(resolution == max(resolution))) +
        geom_text_repel(aes(label = neighbor),
                        nudge_x = 1,
                        na.rm = TRUE,
                        force = 1,
                        box.padding = 1,
                        segment.alpha = .5,
                        data = df_results %>%
                            filter(reference == ref) %>%
                            group_by(neighbor) %>%
                            filter(resolution == min(resolution)))

    # pdf(file = paste0("results/tripet/types/",ref,".pdf"), width = 7, height = 7)
    # print(ref_plots[[ref]])
    # dev.off()
}

grid.arrange(grobs = ref_plots)

## investigate subsets
names(subset.list)
subset.list[['14_near_14']][1:15]

## check if the indexes match for the first rows
df_ac <- read.csv2(file = "data/df_pos_type.csv", row.names = 1, sep=',')
head(meta)
df_ac[which(df_ac$type == 14),] %>% head()

## the cells in each subset are unique to that subset
length(unlist(subset.list))
sum(sapply(subset.list, length))
## there are less cells in the subset than in the whole data
dim(meta)

## the cells in the subsets are not in the orginal data
all(unlist(subset.list) %in% rownames(meta))
unlist(subset.list)[1] %in% rownames(meta)
max(as.numeric(unlist(subset.list)))

## the cells are different
hist(as.numeric(unlist(subset.list)))
hist(as.numeric(rownames(meta)))
hist(as.numeric(rownames(df_ac)))



## plot all cells
library(tidyverse)
df_ac <- read.csv2(file = "data/df_pos_type.csv", row.names = 1, sep=',')
df_ac$subset <- "other"
head(df_ac)

subset.list <- readRDS("data/subsets.near.rds")
names(subset.list)

for (ss in names(subset.list)) {
    df_ac[subset.list[[ss]], "subset"] <- ss
}
df_ac$subset <- as.factor(df_ac$subset)
head(df_ac)

## make sure the coordinates are numeric
df_ac <- df_ac %>%
    dplyr::mutate_at(vars(x, y), as.numeric)
summary(df_ac)

plot_pos_tripet <- df_ac %>%
    ggplot() +
    geom_point(aes(-x, y, color = subset),
               size = .1) +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    scale_colour_manual(values=c("red", "purple", "darkblue", "yellow", "lightgray"))
dev.off()

cellinfo <- df_ac
colors <- c("14" = "red", "18" = "yellow", "other" = "lightgray")


pdf("results/tripet/subset_space.pdf", width = 14, height = 6)
grid.arrange(plot_pos_tripet,
             plot_pos_clusters(names(colors), colors),
             nrow = 1)
dev.off()

## umap
emb <- read.csv('data/Xenium/projection.csv.gz', row.names = 1)
head(emb)
plot(emb, pch=".", main='UMAP')

plot_emb_tripet <- emb %>%
    filter(rownames(emb) %in% rownames(df_ac)) %>%
    mutate(subset = df_ac[which(rownames(df_ac) %in% rownames(emb)), ]$subset) %>%
    ggplot() +
    geom_point(aes(UMAP.1, UMAP.2, color = subset),
               size = .1) +
    scale_color_manual(values = c("red", "purple", "darkblue", "yellow", "lightgray")) +
    guides(color = guide_legend(override.aes = list(size = 10)))

plot_emb_clusters(names(colors), colors)

pdf("results/tripet/subset_emb.pdf", width = 14, height = 6)
grid.arrange(plot_emb_tripet,
             plot_emb_clusters(names(colors), colors),
             nrow = 1)
dev.off()



# Adjust subset -----------------------------------------------------------

## full path to the subset dataset list rds
subsetDataPath <- "data/subsets.near.1000.rds"
subsetType <- "near"
subdist <- 1000
subThresh <- 0.05
ncs <- 16
meta <- read.csv2(file = "data/df_pos_type.csv", row.names = 1, sep=',')

## filter for clusters
head(meta)
meta <- meta[which(meta$type %in% c(14, 18)), ]
head(meta)

meta <- meta[,c("x", "y", "type")]
meta <- meta %>%
    dplyr::mutate_at(vars(x, y), as.numeric)
pos <- meta[,c("x", "y")]
celltypes <- meta[,3]
if(length(levels(celltypes)) == 0){
    message("Warning: `celltypes` does not have levels. Creating levels from values")
    celltypes <- factor(celltypes)
    names(celltypes) <- rownames(pos)
}
cells <- sp::SpatialPointsDataFrame(
    coords = as.data.frame(pos),
    data = data.frame(
        celltypes = celltypes
    ))
cells <- sf::st_as_sf(cells)
rownames(cells) <- rownames(meta)
sf::st_agr(cells) <- "constant"

if(file.exists(subsetDataPath)){
    loadSubset <- subsetDataPath
    saveSubset <- NA
} else {
    loadSubset <- NA
    saveSubset <- subsetDataPath
}

if(assertthat::is.string(loadSubset)){
    subset.list <- readRDS(file = loadSubset)
} else {
    ## parallel shuffling of the grids in each resolution
    subset.list <- getSubsets(cells = cells,
                              sub.dist = subdist,
                              sub.type = subsetType,
                              sub.thresh = subThresh,
                              ncores = ncs,
                              verbose = TRUE)
    if(assertthat::is.string(saveSubset)){
        saveRDS(object = subset.list, file = saveSubset)
    }
}

## Visualize
## get pos data
library(tidyverse)
df_ac <- read.csv2(file = "data/df_pos_type.csv", row.names = 1, sep=',')
df_ac$subset <- "other"
head(df_ac)

subset.list <- readRDS("data/subsets.near.1000.rds")
names(subset.list)

for (ss in names(subset.list)) {
    df_ac[subset.list[[ss]], "subset"] <- ss
}
df_ac$subset <- as.factor(df_ac$subset)
head(df_ac)

df_ac <- df_ac %>%
    dplyr::mutate_at(vars(x, y), as.numeric)
summary(df_ac)

## plot pos
plot_pos_tripet <- df_ac %>%
    ggplot() +
    geom_point(aes(-x, y, color = subset),
               size = .1) +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    scale_colour_manual(values=c("red", "purple", "darkblue", "yellow", "lightgray"))
dev.off()

cellinfo <- df_ac
colors <- c("14" = "red", "18" = "yellow", "other" = "lightgray")


pdf("results/tripet/subset_space_1000.pdf", width = 14, height = 6)
grid.arrange(plot_pos_tripet,
             plot_pos_clusters(names(colors), colors),
             nrow = 1)
dev.off()




# DE Analysis -------------------------------------------------------------

library(tidyverse)
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

library(MERINGUE)

## apply normalization
mtx <- normalizeCounts(counts = cd,
                       log=FALSE,
                       verbose=TRUE)
colnames(mtx) <- colnames(cd)
rownames(mtx) <- rownames(cd)

## filter cells from clusters
dim(mtx)
subset.list <- readRDS("data/subsets.near.rds")
length(unlist(subset.list))
mtx <- mtx[, unlist(subset.list)]
dim(mtx)
clusters_subset <- sapply(names(subset.list),
                          function(ns){
                              vs <- rep(ns, length(subset.list[[ns]]))
                              names(vs) <- subset.list[[ns]]
                              return(vs)
                              })

clusters_subset <- purrr::flatten(clusters_subset)
all(names(clusters_subset) == colnames(mtx))

## create factor compatible with function
tcs <- as.character(clusters_subset)
names(tcs) <- names(clusters_subset)
tcs <- factor(tcs)
clusters_subset <- tcs
head(clusters_subset)

## identify significantly differentially upregulated genes
dg <- getDifferentialGenes(mtx, clusters_subset)
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

## visualize hdg
mean_exp <- function(selected_clusters, selected_genes){

    list_exps <- list()
    for (name_cg in selected_clusters) {
        cells_cluster <- names(clusters_subset[clusters_subset %in% name_cg])
        mgexp <- sapply(selected_genes, function(g){
            mean(mtx[g, cells_cluster] > 0)
        })
        list_exps[[name_cg]] <- mgexp
    }

    return(do.call(cbind, list_exps))
}

## fix figures
pdf('results/tripet/s100_dgh14n18.pdf', height = 9, width = 7)
mean_exp(names(subset.list), dgh.sig[['14_near_18']]) %>% heatmap()
dev.off()

## get spatial and emb exp
pdf('results/tripet/s100_dgh14n14.pdf', height = 18, width = 14)
mean_exp(names(subset.list), dgh.sig[['14_near_14']]) %>% heatmap()
dev.off()

## get spatial and emb exp
pdf('results/tripet/s100_dgh18n14.pdf', height = 18, width = 14)
mean_exp(names(subset.list), dgh.sig[['18_near_14']]) %>% heatmap()
dev.off()

## get spatial and emb exp
pdf('results/tripet/s100_dgh18n18.pdf', height = 18, width = 14)
mean_exp(names(subset.list), dgh.sig[['18_near_18']]) %>% heatmap()
dev.off()
