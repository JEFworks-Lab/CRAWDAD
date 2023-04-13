## run intro.rdm first


# Subcluster --------------------------------------------------------------

library(tidyverse)
library(sp)
library(sf)

df <- do.call(rbind, st_geometry(pkhl)) %>%
  as_tibble() %>% setNames(c('x', 'y'))
df$type <- pkhl$celltypes
head(df)

df %>% ggplot(aes(x, y, color = type)) +
  geom_point(size = .1)

# DBSCAN ------------------------------------------------------------------

library(dbscan)
set.seed(42)

pos <- df[, c('x', 'y')]
head(pos)

## manually adjust eps and minPts
dbscan_res <- dbscan(pos, eps = 75, minPts = 25)

ct_dbscan_res <- lapply(unique(df$type), function(ct){
  df %>% 
    filter(type == ct) %>% 
    select(x, y) %>% 
    dbscan(eps = 75, minPts = 25)
})

df %>% ggplot(aes(x, y, color = as.factor(dbscan_res$cluster))) +
  geom_point(size = .1) +
  guides(color = guide_legend(override.aes = list(size = 7)))

df %>%
  mutate(dcluster = as.factor(dbscan_res$cluster)) %>%
  ggplot(aes(x, y)) +
  geom_point(size = .01) +
  facet_wrap(dcluster ~ .)

df %>% ggplot(aes(x, y, color = as.factor(dbscan_res$cluster) == 0)) +
  geom_point(size = .1) +
  guides(color = guide_legend(override.aes = list(size = 10)))

df %>% ggplot(aes(x, y, color = as.factor(dbscan_res$cluster) == 1)) +
  geom_point(size = .1) +
  guides(color = guide_legend(override.aes = list(size = 10)))

df %>% ggplot(aes(x, y, color = as.factor(dbscan_res$cluster) == 3)) +
  geom_point(size = .1) +
  guides(color = guide_legend(override.aes = list(size = 10)))
