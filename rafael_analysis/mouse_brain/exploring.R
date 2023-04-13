library(tidyverse)

# Load the data -----------------------------------------------------------

df <- read.csv('data/S2R1_with_structure_id_name_and_cell_type_clustering_11-18-2021.csv.gz')
dim(df)
df[1:5, 1:20]

meta <- df[, 1:17]
head(meta)
gexp <- df[, 18:ncol(df)]
gexp[1:5, 1:5]

meta %>% 
  ggplot() +
  geom_point(aes(x = center_x, y = center_y, color = leiden), 
             size = 0.1)

meta %>% 
  filter(str_detect(leiden, 'Astrocytes')) %>% 
  ggplot() +
  geom_point(aes(x = center_x, y = center_y, color = leiden == 'Astrocytes(1)'), 
             size = 0.1)

# Subcluster --------------------------------------------------------------

df <- meta[, c('center_x', 'center_y', 'leiden')] %>% 
  filter(leiden == 'Astrocytes(1)')
colnames(df) <- c('x', 'y', 'type')

df %>% ggplot(aes(x, y, color = type)) +
  geom_point(size = .1)

# DBSCAN ------------------------------------------------------------------

library(dbscan)
set.seed(42)

pos <- df[, c('x', 'y')]
head(pos)

## manually adjust eps and minPts
dbscan_res <- dbscan(pos, eps = 75, minPts = 5)

df %>% ggplot(aes(x, y, color = as.factor(dbscan_res$cluster))) +
  geom_point(size = .1) +
  guides(color = guide_legend(override.aes = list(size = 7)))

## testing optics
optics_res <- optics(pos, eps = 75, minPts = 25)
optics_res <- extractDBSCAN(optics_res, eps_cl = 0.4)

df %>% ggplot(aes(x, y, color = as.factor(optics_res$cluster))) +
  geom_point(size = .1) +
  guides(color = guide_legend(override.aes = list(size = 7)))
