library('crawdad')
ncores <- 7



# Generate data -----------------------------------------------------------

#### simulate
## cells
set.seed(1)
size <- 4000
x <- runif(size, min = 0, max = 1000)
y <- runif(size, min = 0, max = 1000)
p <- data.frame(x = x, y = y, type='D')
rownames(p) <- paste0('cell', 1:size)

## structures

## large A circles
as <- c(250, 250, 750, 750)*2
bs <- c(250, 750, 250, 750)*2
invisible(sapply(1:4, function(i) {
  a <- as[i]
  b <- bs[i]
  ro <- 150*2
  co <- 'A'
  po <- 1
  c1o <- rownames(p[((p$x-a)^2 + (p$y - b)^2 < ro^2),])
  p[c1o,]$type <<- sample(co, size = length(c1o), replace = TRUE, prob = po)
}))

## B blobs inside A blobs
invisible(sapply(1:4, function(i) {
  ro <- 50*2 # 80
  # co <- 'B'
  co <- c('B', 'C')
  # po <- 1
  po <- c(0.5, 0.5)
  
  ## inside structure
  a <- as[i]
  b <- bs[i]
  c1o <- rownames(p[((p$x-a)^2 + (p$y - b)^2 < ro^2),])
  p[c1o,]$type <<- sample(co, size = length(c1o), replace = TRUE, prob = po)
}))

## visualize
plt <- vizAllClusters(cells = p,
                      coms = p$type,
                      title = "sim",
                      axisAdj = 1, s = 6, a = 0.5) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 1))
plt

#### convert to SP
cells <- crawdad::toSP(pos = p[,c("x", "y")],
                       celltypes = p$type)
plot(cells)



# Analyze with CRAWDAD ----------------------------------------------------

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           #scales = 10^(seq(1.7, 3.7, by = 0.2)),
                                           scales = seq(100, 1000, by=50),
                                           perms = 10,
                                           ncores = ncores,
                                           seed = 1,
                                           verbose = TRUE)
## find trends, passing background as parameter
results <- crawdad::findTrends(cells,
                               dist = 100,
                               shuffle.list = shuffle.list,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)

dat <- crawdad::meltResultsList(results, withPerms = T)



# Visualize with CRAWDAD --------------------------------------------------

vizColocDotplot(dat)

library(tidyverse)
dat_filter <- dat %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'C')

vizTrends(dat_filter, lines = T, withPerms = T)

dat_filter %>% group_by(neighbor, scale, reference, id) %>% 
  summarise(mean = mean(Z), 
            sd = sd(Z)) %>% 
  mutate(Z = mean)
