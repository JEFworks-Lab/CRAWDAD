
#' get neighbor cells defined as being a distance away from a set of reference cells
#' @description get neighbor cells defined as being a distance away from a set of reference cells.
#'      `reference.ids` can be selected by subsetting rownames from `cells`:
#'       ex: rownames(cells)[which(cells$celltypes == "A")]
#'       or can be an entry in a subset list from `selectSubsets()`
#' 
#' @param cells sf object, with celltypes features and point geometries
#' @param reference.ids vector of cell ids (rownames) in `cells` to be used as the reference cell set
#' @param dist distance to define neighbors (default = 100)
#' @param returnSF boolean to return either an sf object of just the neighbors
#' otherwise returns factor where non neighbor cells are NAs. (default: FALSE)
#'
#' @return sf object of the neighbor cells or factor of neighbor cells
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffle.list <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000, 1500, 2000), ncores = 2)
#' binomMat <- binomialTestMatrix(cells, neigh.dist = 100, ncores = 2)
#' subset.list <- selectSubsets(binomMat, cells$celltypes, sub.type = "near", sub.thresh = 0.05)
#' neighCells <- getNeighbors(cells = cells, reference.ids = subset.list[["C_near_B"]],  dist = 100)
#' }
#' 
#' @export
getNeighbors <- function(cells,
                         reference.ids,
                         removeRef = TRUE,
                         dist = 100,
                         returnSF = FALSE){
  
  ## get the reference cells
  ref.cells <- cells[reference.ids,]
  
  ## define the buffer around the reference cells with distance of
  refs.buffer <- sf::st_buffer(ref.cells, dist)
  
  ## get the neighbor cells
  cells.inbuffer <- sf::st_intersection(cells, refs.buffer$geometry)
  ## remove duplicate neighbors
  neigh.cells <- cells.inbuffer[intersect(rownames(cells.inbuffer), rownames(cells)),]
  
  ## remove the reference cells from the set of neighbor cells
  if(removeRef){
    neigh.cells <- neigh.cells[setdiff(rownames(neigh.cells), reference.ids),]
  }
  
  ## get the new com factor of neighbor cells only
  neigh_coms <- cells$celltypes
  names(neigh_coms) <- rownames(cells)
  ## set the non neighbor cell labels to NA
  non_neighbors <- setdiff(rownames(cells), rownames(neigh.cells))
  neigh_coms[non_neighbors] <- NA
  
  ## if true, the output is sf of just the neighbor cells
  ## otherwise, the output is a new com factor of all the cells in `cells` but
  ## non neighbors are NAs
  if(returnSFE){
    return(neigh.cells)
  } else {
    return(neigh_coms)
  }
}


#' Melt the output list of `findTrends()` into a dataframe
#' 
#' @description idea is that the output of `findTrends()` is a list of dataframes, where each dataframe
#' is for a reference cell type and contains the Z scores at each scale for the neighbor cell types.
#' So melt this list of dataframes into a single dataframe. Idea is to get a single dataframe setup for plotting with
#' ggplot2 and tidyverse functions. `id` parameter allows adding a specific identifier for the given melted results so that
#' one can combine multiple results dataframes and compare downstream.
#' For example, you melt a results list from `findTrends()` that was generated from an analysis using a neighbor distance of 100.
#' Then you generate another one that was done with neighbor distance of 200. 
#' The id column for each dataframe can be set to "100" and "200", respectively.
#' Then both dataframes can be combined into one final dataframe.
#' Now you have identifiers that include: scale, neighbor, reference, and "id" (ie neighbor distance).
#' 
#' @param resultsList list output from `findTrends()`
#' @param id id desired, can add a column that contains an additional identifier for the results. Can use these for plotting and comparing different things
#' @param withPerms if the results list is a list of lists using `returnMeans = FALSE` in `findTrends()`, then column order is different and this flag is needed (default: FALSE)
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffle.list <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000, 1500, 2000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffle.list = shuffle.list, ncores = 2)
#' meltResultsList(results)
#' }
#' 
#' @export
meltResultsList <- function(resultsList, id = NA, withPerms = FALSE){
  
  df <- reshape2::melt(resultsList)
  
  if(withPerms){
    colnames(df) <- c("perm", "neighbor", "Z", "scale", "reference")
  } else {
    colnames(df) <- c("scale", "neighbor", "Z", "reference")
  }
  
  ## add an identifier for the particular results
  df[["id"]] <- id
  
  # scales as numeric:
  df <- df %>%
    dplyr::mutate_at(dplyr::vars(scale), as.numeric)
  
  return(df)
  
}


#' filter for significant cell type association trends that are co-localized
#' 
#' @description filter the results list from `findTrends()` for neighbor cell types that are significantly co-localized with each reference cell type.
#' 
#' @param results list output from `findTrends()`
#' @param alpha significance threshold
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffle.list <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000, 1500, 2000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffle.list = shuffle.list, ncores = 2)
#' filterCoTrends(results = results, alpha = 0.05)
#' }
#' 
filterCoTrends <- function(results, alpha = 0.05) {
  zthresh <- qnorm(1-alpha/2)
  lapply(results, function(x) {
    colIds <- unique(which(x > zthresh, arr.ind=TRUE)[,2])
    celltypes <- colnames(x)[colIds]
    x <- as.matrix(x[,celltypes])
    colnames(x) <- celltypes
    x
  })
}


#' filter for significant cell type association trends that are separated
#' 
#' @description filter the results list from `findTrends()` for neighbor cell types that are significantly separated with each reference cell type.
#' 
#' @param results list output from `findTrends()`
#' @param alpha significance threshold
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffle.list <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000, 1500, 2000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffle.list = shuffle.list, ncores = 2)
#' filterSepTrends(results = results, alpha = 0.05)
#' }
#' 
filterSepTrends <- function(results, alpha = 0.05) {
  zthresh <- qnorm(1-alpha/2)
  lapply(results, function(x) {
    colIds <- unique(which(x < -zthresh, arr.ind=TRUE)[,2])
    celltypes <- colnames(x)[colIds]
    x <- as.matrix(x[,celltypes])
    colnames(x) <- celltypes
    x
  })
}


#' filter for significant cell type association trends that are either co-localized or separated
#' 
#' @description filter the results list from `findTrends()` for neighbor cell types that become significantly co-localized and separated at different scales with each reference cell type.
#' 
#' @param results list output from `findTrends()`
#' @param alpha significance threshold
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffle.list <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000, 1500, 2000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffle.list = shuffle.list, ncores = 2)
#' filterChangeTrends(results = results, alpha = 0.05)
#' }
#' 
filterChangeTrends <- function(results, alpha = 0.05) {
  zthresh <- qnorm(1-alpha/2)
  lapply(results, function(x) {
    co <- unique(which(x > zthresh, arr.ind=TRUE)[,2])
    sep <- unique(which(x < -zthresh, arr.ind=TRUE)[,2])
    colIds <- intersect(co, sep)
    celltypes <- colnames(x)[colIds]
    x <- as.matrix(x[,celltypes])
    colnames(x) <- celltypes
    x
  })
}



#' Perform Bonferroni Multiple-test Correction
#' 
#' @description
#' Calculate the number of tests performed by multiplying the number of cell
#' types used in the analysis. Then, divide the 0.05 pvalue by the number of 
#' tests. Finally, convert the new pvalue to a zscore.
#' 
#' @param dat data.frame; the melted findTrends result.
#' @param pval numeric; the pvalue to be converted. Defaults to 0.05.
#' 
#' @return numeric; the Bonferroni corrected zscore.
#' 
#' @export
#' 
correctZBonferroni <- function(dat, pval = 0.05) {
  ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
  psig <- pval/ntests
  zsig <- round(qnorm(psig/2, lower.tail = F), 2)
  
  return(zsig)
}
