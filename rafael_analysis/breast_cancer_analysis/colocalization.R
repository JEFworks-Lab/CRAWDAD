
library(msctcla)
library(tidyverse)

# Load data ---------------------------------------------------------------

df <- read.csv("df.csv", row.names = 1)
head(df)
dim(df)

df <- df[, c('center_x', 'center_y', 'celltype')]
colnames(df) <- c('x', 'y', 'type')
head(df)

# Cell-type co-localization analysis --------------------------------------

## parameters fir computing trends
neighbor_distance <- 100 # each distance to define neighbors
resolutions <- c(100, 250, 500, 1000) # different shuffling resolutions
permutations <- 1 # number of permutations for shuffling;
seed <- 1
n_cores <- 16
positions <- df
Sys.time()
results <- findTrends(pos = positions[,c("x", "y")],
                      celltypes = positions$type,
                      resolutions,
                      dist = neighbor_distance,
                      perms = permutations,
                      ncores = n_cores,
                      verbose = TRUE,
                      saveShuffleFile = 'test.rds',
                      seed = seed)
Sys.time()
## process data
df_results <- meltResultsList(resultsList = results)

head(df_results)
write.csv(df_results, "crawdad_results.csv")

## Convert to factor to plot the data
df_results$neighbor <- as.factor(df_results$neighbor)
df_results$reference <- as.factor(df_results$reference)

### visualize
## colorblind-friendly palette: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                        "#D55E00", "#CC79A7", "#000000", "#F0E442")

# plotTrends(results = df_results, figPath = 'results/trends.pdf',
#            width = 24, height = 24)

## The Z scale of the figure hides the trends because there is a Z-score of >1000
## On the simulated data, they ranged from -15 to +15
pdf(file = "grid.pdf", width = 24, height = 24)
df_results %>% ggplot() +
    geom_line(aes(x=resolution, y=Z, color=neighbor), linewidth=1.3) +
    scale_color_manual(values = cbPalette) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    facet_grid(reference ~ .,
               labeller = labeller(.rows = label_both, .cols = label_both),
               scales="free") +
    geom_hline(yintercept=c(2,-2), linetype="dashed",
               color = "blue", size=.3)

dev.off()

## Clusters
df %>%
    ggplot() +
    geom_point(aes(x, y, color = as.factor(type)),
               size = .5) +
    labs(x='x', y='y',
         title='Clusters in space') +
    facet_wrap(type ~ .) +
    scale_color_manual(values = cbPalette) +
    guides(color = guide_legend(override.aes = list(size=5)))

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
        scale_colour_manual(values=setNames(df$color, df$type)) +
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

    pdf(file = paste0("results/types/",ref,".pdf"), width = 7, height = 7)
    print(ref_plots[[ref]])
    dev.off()
}

## all references
## grid arrange
pdf(file = "results/references_grid.pdf", width = 7*5, height = 7*5)
grid.arrange(grobs = ref_plots)
dev.off()

## facet grid
pdf(file = "results/references_list.pdf", width = 7, height = 7*length(unique(df_results$reference)))
df_results %>%
    filter(reference == reference) %>%
    ggplot(aes(x=resolution, y=Z, color=neighbor)) +
    geom_line() +
    scale_colour_manual(values=setNames(df$color, df$type)) +
    theme(legend.position = "none") +
    geom_text_repel(aes(label = neighbor),
                    nudge_x = 1,
                    na.rm = TRUE,
                    force = 1,
                    box.padding = 1,
                    segment.alpha = .5,
                    data = df_results %>%
                        filter(reference == reference) %>%
                        group_by(neighbor) %>%
                        filter(resolution == max(resolution))) +
    geom_text_repel(aes(label = neighbor),
                    nudge_x = 1,
                    na.rm = TRUE,
                    force = 1,
                    box.padding = 1,
                    segment.alpha = .5,
                    data = df_results %>%
                        filter(reference == reference) %>%
                        group_by(neighbor) %>%
                        filter(resolution == min(resolution))) +
    facet_grid(reference ~ .,
               scales="free")
dev.off()



# Visualize cell types ----------------------------------------------------

df %>%
    filter(type %in% c(12, 11, 14, 6)) %>%
    ggplot() +
    geom_point(aes(x=x, y=y, color=type), size=.01) +
    scale_colour_manual(values=setNames(df$color, df$type)) +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    facet_grid()
