library(tidyverse)

source('init.R')

head(dat)
## select the smallest resolution in which the z score becomes significant 
## it only selects the positive Z scores and disregard the negative ones
## so, there are many results that do not show
sig_dat <- dat %>%
  filter(Z >= 2) %>% 
  group_by(neighbor, reference) %>% 
  filter(resolution == min(resolution, na.rm = TRUE))
head(sig_dat)

sig_dat %>% 
  ggplot() +
  geom_tile(aes(x=reference, y=neighbor, fill=resolution)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## because of the filter in the previous part, there are many NAs in the matrix
## the heatmap cannot be created
# sig_mat <- sig_dat %>% 
#   select(!c(Z, id)) %>% 
#   pivot_wider(names_from = neighbor, values_from = resolution) %>% 
#   as.data.frame()
# rownames(sig_mat) <- sig_mat$reference
# sig_mat <- sig_mat %>% 
#   subset(select = -reference) %>% 
#   as.matrix()
# sig_mat %>% heatmap()

## plotting Z and resolution
sig_dat %>% 
  ggplot() +
  geom_point(aes(x=reference, y=neighbor, 
                 color=resolution, size=Z)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## smooth
dat %>% 
  filter(reference == 'B',
         neighbor == 'Endothelial') %>% 
  ggplot(aes(x=resolution, y=Z)) + 
  geom_line() + 
  geom_point() +
  geom_smooth() +
  labs(title = 'Ref B vs Neigh Endothelial')

## smooth function
cts <- unique(dat$reference)
x_vals = seq(100, 2500, 1)
sig_res <- list()
for (ref in cts) {
  for (nei in cts) {
    lo_df <- dat %>% 
      filter(reference == ref,
             neighbor == nei)
    ## why are there nan values
    if (!is.nan(lo_df$Z[1])) {
      lo <- loess(lo_df$Z~lo_df$resolution)
      
      y_pred = predict(lo, x_vals)
      res_sig = which(y_pred >= 2)[1] + 100
      sig_res[[ref]][[nei]] = res_sig
    }
  }
}

## create df
sig_df <- as.data.frame(do.call(rbind, sig_res)) %>% 
  tibble::rownames_to_column(var = 'reference') %>% 
  pivot_longer(!reference, names_to = 'neighbor', values_to = 'resolution') %>% 
  mutate(resolution = as.integer(resolution)) %>% 
  drop_na(resolution) %>% 
  as.data.frame()
head(sig_df)

## plot tiles
sig_df %>% 
  ggplot() +
  geom_tile(aes(x=reference, y=neighbor, fill=resolution)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## check if not wrong: it is right!
# colnames(sig_df) <- c('neighbor', 'reference', 'resolution')
# sig_df %>% 
#   ggplot() +
#   geom_tile(aes(x=reference, y=neighbor, fill=resolution)) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## investigate some data
unique(dat$reference)
ct1 <- 'B' 
ct2 <- 'Endothelial'
dat %>% 
  filter(reference == ct1,
         neighbor == ct2) %>% 
  ggplot(aes(x=resolution, y=Z)) + 
  geom_line() + 
  geom_point() +
  geom_smooth() +
  labs(title = paste('Ref', ct1, 'vs Neigh', ct2))

foo <- data.sub$cell_type
foo[!(foo %in% c(ct1, ct2))] <- NA
ggplot(data.sub, aes(x=x, y=y, col=foo)) + 
  geom_point(size=0.1) + 
  scale_colour_manual(
    values = c('blue', 'red'),
    na.value = "#eeeeee"
  ) + theme_minimal()