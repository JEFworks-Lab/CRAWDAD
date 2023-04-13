cellinfo %>%
    filter(cell_id %in% names(clusters)) %>%
    mutate(cluster = clusters) %>%
    ggplot() +
    geom_point(aes(-x_centroid, y_centroid,
                   color = (cluster %in% c(14))),
               size = .1) +
    scale_color_manual(values = c("#ffc6e8", "#ed008c")) +
    cowplot::theme_nothing() +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x = NULL, y = NULL)
