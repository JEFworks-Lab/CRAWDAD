## visualize hdg
plot_pos_gene <- function(selected_gene, selected_cells){
    df %>%
        filter(cid %in% selected_cells) %>%
        mutate(gene = log1p(cd[selected_gene, selected_cells])) %>%
        ggplot() +
        geom_point(aes(x, y, color = gene),
                   size = .25) +
        scale_color_continuous(low = "lightgray", high = "black") +
        labs(color = selected_gene)
}

plot_pos_cluster <- function(cluster_column, selected_cells) {
    df %>%
        filter(cid %in% selected_cells) %>%
        ggplot(aes(x, y, color = .data[[cluster_column]])) +
        geom_point(size = .1) +
        guides(color = guide_legend(override.aes = list(size = 10)))
}
