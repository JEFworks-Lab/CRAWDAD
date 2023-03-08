
#' get neighbor cells defined as being a distance away from a set of reference cells
#' @description get neighbor cells defined as being a distance away from a set of reference cells.
#'      `reference.ids` can be selected by subsetting rownames from `cells`:
#'       ex: rownames(cells)[which(cells$celltypes == "A")]
#'       or can be an entry in a subset list from `selectSubsets()`
#' 
#' @param cells sp::SpatialPointsDataFrame object, with celltypes features and point geometries
#' @param reference.ids vector of cell ids (rownames) in `cells` to be used as the reference cell set
#' @param dist distance to define neighbors (default = 100)
#' @param returnSP boolean to return either an sp::SpatialPointsDataFrame object of just the neighbors
#' otherwise returns factor where non neighbor cells are NAs. (default: FALSE)
#'
#' @return sp::SpatialPointsDataFrame object of the neighbor cells or factor of neighbor cells
#' 
#' @export
getNeighbors <- function(cells,
                         reference.ids,
                         removeRef = TRUE,
                         dist = 100,
                         returnSP = FALSE){
  
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
  
  ## if true, the output is sp::SpatialPointsDataFrame of just the neighbor cells
  ## otherwise, the output is a new com factor of all the cells in `cells` but
  ## non neighbors are NAs
  if(returnSP){
    return(neigh.cells)
  } else {
    return(neigh_coms)
  }
}


#' Melt the output list of `findTrends()` into a dataframe
#' 
#' @description idea is that the output of `findTrends()` is a list of dataframes, where each dataframe
#' is for a reference cell type and contains the Z scores at each resolution for the neighbor cell types.
#' So melt this list of dataframes into a single dataframe. Idea is to get a single dataframe setup for plotting with
#' ggplot2 and tidyverse functions. `id` parameter allows adding a specific identifier for the given melted results so that
#' one can combine multiple results dataframes and compare downstream.
#' For example, you melt a results list from `findTrends()` that was generated from an analysis using a neighbor distance of 100.
#' Then you generate another one that was done with neighbor distance of 200. 
#' The id column for each dataframe can be set to "100" and "200", respectively.
#' Then both dataframes can be combined into one final dataframe.
#' Now you have identifiers that include: resolution, neighbor, reference, and "id" (ie neighbor distance).
#' 
#' @param resultsList list output from findTrendsv2
#' @param id id desired, can add a column that contains an additional identifier for the results. Can use these for plotting and comparing different things
#' 
#' @export
meltResultsList <- function(resultsList, id = NA){
  
  df <- reshape2::melt(resultsList)
  colnames(df) <- c("resolution", "neighbor", "Z", "reference")
  ## add an identifier for the particular results
  df[["id"]] <- id
  return(df)
  
}