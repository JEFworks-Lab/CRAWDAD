library(tidyverse)
library(crawdad)
library(parallel)

data <- read.csv('data/single_cells.csv.gz', row.names=1)
dim(data)
head(data)

## note there are multiple samples
table(data$donor, data$region_num)
## multiple samples from same patient
table(data$array, data$region_num)

## variables looping
sample_ids <- data %>% group_by(array, region_num) %>% summarise()
plot_list <- list()

## looping through samples
for (i in 1:nrow(sample_ids)) {
  
  sample_array = sample_ids[i, 1] %>% as.character()
  sample_region_num = sample_ids[i, 2] %>% as.character()
  
  vi <- data$array == sample_array & data$region_num == sample_region_num
  data.sub <- data[vi,]
  cells <- crawdad:::toSP(pos = data.sub[,c("x", "y")],
                          celltypes = data.sub[,'cell_type'])
  
  ## generate background
  shuffle.list <- crawdad:::makeShuffledCells(cells,
                                              resolutions = c(100, 200, 300, 400,
                                                              500, 600, 700, 800,
                                                              900, 1000, 1250,
                                                              1500, 1750, 2000,
                                                              2250, 2500),
                                              perms = 1,
                                              ncores = 16,
                                              seed = 1,
                                              verbose = TRUE)
  file_name <- paste0('results/data/shuffle.list_', 
                      sample_array, '_', sample_region_num, '.RDS')
  saveRDS(shuffle.list, file_name)
  shuffle.list <- readRDS(file_name)
  
  ## find trends, passing background as parameter
  results <- crawdad::findTrends(cells,
                                 dist = 100,
                                 shuffle.list = shuffle.list,
                                 ncores = 16,
                                 verbose = TRUE)
  file_name <- paste0('results/data/trends_', 
                      sample_array, '_', sample_region_num, '.RDS')
  saveRDS(results, file_name)
  results <- readRDS(file_name)
  dat <- crawdad::meltResultsList(results, id = "dist_100")
  
  ## plotting
  sig_dat <- dat %>%
    filter(Z >= 2) %>% 
    group_by(neighbor, reference) %>% 
    filter(resolution == min(resolution, na.rm = TRUE))
  plot_list[[sample_array]][[sample_region_num]] <- 
    sig_dat %>% 
    ggplot() +
    geom_tile(aes(x=reference, y=neighbor, fill=resolution)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

