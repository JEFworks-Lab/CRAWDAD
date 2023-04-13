
# Load data ---------------------------------------------------------------

## load all data
source('load_data.R')

## filter for clusters 14 and 18
df <- df[which(df$type %in% c(14, 18)), ]
head(df)

df %>% ggplot(aes(x, y, color = type)) +
    geom_point(size = .5)


# DBSCAN ------------------------------------------------------------------

library(dbscan)
set.seed(42)

pos <- df[, c('x', 'y')]
head(pos)

## manually adjust eps and minPts
dbscan_res <- dbscan(pos, eps = 75, minPts = 25)

df %>% ggplot(aes(x, y, color = as.factor(dbscan_res$cluster))) +
    geom_point(size = .1) +
    guides(color = guide_legend(override.aes = list(size = 7)))

df %>%
    mutate(dcluster = as.factor(dbscan_res$cluster)) %>%
    ggplot(aes(x, y)) +
    geom_point(size = .5) +
    facet_wrap(dcluster ~ .)

df %>% ggplot(aes(x, y, color = as.factor(dbscan_res$cluster) == 6)) +
    geom_point(size = .1) +
    guides(color = guide_legend(override.aes = list(size = 10)))



# Dclusters cell types ----------------------------------------------------

df <- df %>%
    mutate(dcluster = as.factor(dbscan_res$cluster)) %>%
    filter(dcluster != 0)
head(df)

dclusters <- unique(df$dcluster)

pct_18 <- sapply(dclusters, function(dc){
    ct <- df[which(df$dcluster == dc), 'type']
    pct <- sum(ct == 18)/length(ct)
    return(pct)
})
names(pct_18) <- dclusters
pct_18

pure_18 <- names(which(pct_18 > .9))
pure_14 <- names(which(pct_18 < .1))
mixed_types <- names(which((pct_18 >= .1) & (pct_18 <= .9)))
mixed_types

head(df)
## map the cell types based on purity and also compare to the original type
## to see if it is a pure cell type in its pure cluster.
map_ct <- as.character(sapply(rownames(df), function(c){
    ct <- df[c, 'type']
    dc <- df[c, 'dcluster']
    if (dc %in% pure_18) {
        if (ct == 18) {
            return('pure18')
        } else {
            return('other')
        }
    } else if (dc %in% pure_14) {
        if (ct == 14) {
            return('pure14')
        } else {
            return('other')
        }
    } else if (dc %in% mixed_types) {
        if (ct == 18) {
            return('mixed18')
        } else if (ct == 14) {
            return('mixed14')
        } else {
            return('other')
        }
    }
}))
map_ct[1:5]
unique(map_ct)

df <- df %>%
    mutate(pure_ct = map_ct)
df$pure_ct <- as.factor(df$pure_ct)
head(df)

df %>%
    ggplot() +
    geom_point(aes(x, y, color = pure_ct), size = .1) +
    guides(color = guide_legend(override.aes = list(size = 10)))



# DEA ---------------------------------------------------------------------

## filter for the data of interest
perform_dge <- function(selected_clusters, cluster_column){
    ## select cells and filter matrix
    cells <- df[which(df[, cluster_column] %in% selected_clusters), 'cid']
    mtx_dg <- mtx[, cells]

    ## filter metadata and get clusters
    df_dg <- df[which(df$cid %in% cells), ]
    clusters_dg <- df_dg[, cluster_column]
    names(clusters_dg) <- df_dg$cid

    print("Unique clusters")
    print(unique(df_dg$pure_ct))
    print(unique(clusters_dg))

    ## perform dea
    dg <- getDifferentialGenes(mtx_dg, clusters_dg)
    dg.sig <- lapply(dg, function(x) {
        x <- x[x$p.adj < 0.05,]
        x <- na.omit(x)
        x <- x[x$Z > 0, ]
        rownames(x)
    })
    print("Genes per cluster")
    print(lapply(dg.sig, length))
    ## highest gene exp across populations
    dgh.sig <- lapply(dg, function(x) {
        x <- x[x$p.adj < 0.05,]
        x <- na.omit(x)
        x <- x[x$Z > 0, ]
        x <- x[x$highest,]
        rownames(x)
    })
    print("Genes with highest expression per cluster")
    print(lapply(dgh.sig, length))

    return(list('dg' = dg, 'dg.sig' = dg.sig, 'dgh.sig' = dgh.sig))
}

dgs <- perform_dge(c('pure18', 'mixed18'), 'pure_ct')
dg.sig <- dgs[['dg.sig']]
dgh.sig <- dgs[['dgh.sig']]
dg <- dgs[['dg']]

# Plotting ----------------------------------------------------------------

## calculate mean expression of selected genes in selected cluster
mean_exp <- function(selected_clusters, cluster_column, selected_genes) {
    list_exps <- list()
    for (name_cg in selected_clusters) {
        cells_cluster <- df[which(df[, cluster_column] == name_cg), 'cid']
        mgexp <- sapply(selected_genes, function(g){
            mean(mtx[g, cells_cluster])
        })
        list_exps[[name_cg]] <- mgexp
        }
    return(do.call(cbind, list_exps))
}

dgh.sig[['pure18']]
pdf('results/dbscan/dgh_pure18.pdf', height = 18, width = 14)
mean_exp(c('pure18', 'mixed18'), 'pure_ct', dgh.sig[['pure18']]) %>%
    heatmap()
dev.off()

dgh.sig[['mixed18']]
pdf('results/dbscan/dgh_mixed18.pdf', height = 18, width = 14)
mean_exp(c('pure18', 'mixed18'), 'pure_ct', dgh.sig[['mixed18']]) %>%
    heatmap()
dev.off()

## check previous code for mean > 0

source('plotting.R')

selected_clusters <- c('pure18', 'mixed18')
selected_cells <- df[which(df[, 'pure_ct'] %in% selected_clusters), 'cid']

plot_pos_cluster('pure_ct', selected_cells)
plot_pos_gene('CD9', selected_cells)

dg[['pure18']] %>%
    filter(rownames(dg[['pure18']]) %in% dgh.sig[['pure18']]) %>%
    arrange(desc(M)) %>%
    head(15)
top4 <- dg[['pure18']] %>%
    filter(rownames(dg[['pure18']]) %in% dgh.sig[['pure18']]) %>%
    arrange(desc(M)) %>%
    head(4) %>%
    rownames()

lp <- lapply(top4, FUN = function(g){
    plot_pos_gene(g, selected_cells)
})
grid.arrange(grobs = lp)

dg[['mixed18']] %>%
    filter(rownames(dg[['mixed18']]) %in% dgh.sig[['mixed18']]) %>%
    arrange(desc(M)) %>%
    head(15)
top4 <- dg[['mixed18']] %>%
    filter(rownames(dg[['mixed18']]) %in% dgh.sig[['mixed18']]) %>%
    arrange(desc(M)) %>%
    head(4) %>%
    rownames()

lp <- lapply(top4, FUN = function(g){
    plot_pos_gene(g, selected_cells)
})
grid.arrange(grobs = lp)



# For cluster 14 ----------------------------------------------------------

dgs <- perform_dge(c('pure14', 'mixed14'), 'pure_ct')
dg.sig <- dgs[['dg.sig']]
dgh.sig <- dgs[['dgh.sig']]
dg <- dgs[['dg']]

dgh.sig[['pure14']]
pdf('results/dbscan/dgh_pure14.pdf', height = 18, width = 14)
mean_exp(c('pure14', 'mixed14'), 'pure_ct', dgh.sig[['pure14']]) %>%
    heatmap()
dev.off()

dgh.sig[['mixed14']]
pdf('results/dbscan/dgh_mixed14.pdf', height = 18, width = 14)
mean_exp(c('pure14', 'mixed14'), 'pure_ct', dgh.sig[['mixed14']]) %>%
    heatmap()
dev.off()

selected_clusters <- c('pure14', 'mixed14')
selected_cells <- df[which(df[, 'pure_ct'] %in% selected_clusters), 'cid']

plot_pos_cluster('pure_ct', selected_cells)

dg[['pure14']] %>%
    filter(rownames(dg[['pure14']]) %in% dgh.sig[['pure14']]) %>%
    arrange(desc(M))

lp <- lapply(dgh.sig[['pure14']], FUN = function(g){
    plot_pos_gene(g, selected_cells)
})
grid.arrange(grobs = lp)

dg[['mixed14']] %>%
    filter(rownames(dg[['mixed14']]) %in% dgh.sig[['mixed14']]) %>%
    arrange(desc(M))
top4 <- dg[['mixed14']] %>%
    filter(rownames(dg[['mixed14']]) %in% dgh.sig[['mixed14']]) %>%
    arrange(desc(M)) %>%
    head(4) %>%
    rownames()

lp <- lapply(top4, FUN = function(g){
    plot_pos_gene(g, selected_cells)
})
grid.arrange(grobs = lp)

