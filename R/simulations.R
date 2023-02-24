#' Create a uniform background of points, where each point is a given cell type
#' 
#' @param size number of points to generate
#' @param cts what are the different cell types in the background? Vector of IDs
#' @param prob the probability, or proportion of each cell type, summing to 1. A vector of probabilities
#' @param seed pseudo random seed for the function sample()
#' @param scale amount to expand the coordinates by
#' 
#' @return data.frame of x, y coordinates of each point, and column type indicating its cell type
#'
simulate_background <- function(size = 10000, cts = c("A"), prob = c(1), seed = 1, scale = 1){
  set.seed(seed)
  x <- runif(size, min = 0, max = 1)
  y <- runif(size, min = 0, max = 1)
  p <- data.frame(x = x, y = y, type = sample(cts, size = size, replace = TRUE, prob = prob)
  )
  
  ## scale to 3100 microns for different tile resolutions
  p$x <- p$x * scale
  p$y <- p$y * scale
  
  return(p)
}


#' Create cell type circle patterns in the background of cells
#'
#' @description takes a dataframe of x,y, and type (from simulate_background, for example)
#' and selects cells to create new cell type patterns of circles.
#' THe circles can actually be an outer ring and an inner core
#'
#' @param pos dataframe with columns x, y, and type (typically from simulate_background)
#' @param locs list of 2-d vectors, where each vector is the x and y coordinates of a circle center (0-1)
#' @param radii list of 2-d lists, where each inner list has variables "inner" and "outer" that refer to the
#' radius of the outer ring and the inner cores of each circle. An inner list for each circle in locs
#' @param cts same format as radii, but "inner" and "outer" are vectors of cell types present in core and ring
#' @param probs same format as radii, but "inner" and "outer" are vectors of cell types proportions 
#' 
simulate_circles <- function(pos, locs, radii, cts, probs){
  
  p <- pos
  
  for(i in 1:length(locs)){
    
    a <- locs[[i]][1]
    b <- locs[[i]][2]
    
    ## outer section of circle, mainly for forming the outside ring
    ro <- radii[[i]]$outer
    co <- cts[[i]]$outer
    po <- probs[[i]]$outer
    c1o <- rownames(p[((p$x-a)^2 + (p$y - b)^2 < ro^2),])
    p[c1o,]$type <- sample(co, size = length(c1o), replace = TRUE, prob = po)
    
    ## inner section of the circle
    ri <- radii[[i]]$inner
    ci <- cts[[i]]$inner
    pi <- probs[[i]]$inner
    c1i <- rownames(p[((p$x-a)^2 + (p$y - b)^2 < ri^2),])
    p[c1i,]$type <- sample(ci, size = length(c1i), replace = TRUE, prob = pi)
    
  }
  
  celltypes <- p$type
  names(celltypes) <- rownames(p)
  p$type <- as.factor(celltypes)
  
  return(p)
  
}