#' convert an sp::SpatialPointsDataFrame object to a dataframe with x y coords and cell type labels
#' @param cells sp::SpatialPointsDataFrame object, with celltypes features and point geometries
#'
#' @return datafame with columns: x, y, and celltype
#'
spToDF <- function(cells){
  
  pos <- data.frame(sf::st_coordinates(sf::st_cast(cells$geometry,"POINT")))
  pos$type <- cells$celltypes
  colnames(pos) <- c("x", "y", "celltypes")
  return(pos)
  
}


#' get neighbor cells defined as being a distance away from a set of reference cells
#' @description get neighbor cells defined as being a distance away from a set of reference cells.
#'      `reference.ids` can be selected by subsetting rownames from `cells`:
#'       ex: rownames(cells)[which(cells$celltypes == "A")]
#'       or can be an entry in a subset list from `selectSubsets()`
#' 
#' @param cells sp::SpatialPointsDataFrame object, with celltypes features and point geometries
#' @param reference.ids vector of cell ids (rownames) in `cells` to be used as the reference cell set
#' @param dist distance to define neighbors (default = 100)
#'
#' @return sp::SpatialPointsDataFrame object of the neighbor cells
#' 
getNeighbors <- function(cells,
                         reference.ids,
                         dist = 100){
  
  ## get the reference cells
  ref.cells <- cells[reference.ids,]
  
  ## define the buffer around the reference cells with distance of
  refs.buffer <- sf::st_buffer(ref.cells, dist)
  
  ## get the neighbor cells
  cells.inbuffer <- sf::st_intersection(cells, refs.buffer$geometry)
  ## remove duplicate neighbors
  neigh.cells <- cells.inbuffer[intersect(rownames(cells.inbuffer), rownames(cells)),]
  
  return(neigh.cells)
  
}


#' Melt the output list of findTrendsv2 into a dataframe
#' 
#' @description idea is that the output is a list of dataframes, where each dataframe
#' is for a reference cell type and contains the Z scores at each resolution for the neighbor cells.
#' So melt this list of dataframes into a single dataframe. Idea is to get a single dataframe setup for plotting with
#' ggplot2 and tidyverse functions
#' @param resultsList list output from findTrendsv2
#' @param id id desired, can add a column that contains an additional identifier for the results.
#' For example, you melted a resultsList that was generated from a simulation using 2 circles. 
#' Then you generate another one that was done with 1 circle. The id column for each dataframe
#' can be set to "2" and "1". Then both dataframes can be combined into one final dataframe
#' Now you have identifiers that include: resolution, neighbor, reference, and simulation type.
#' Can use these for plotting and comparing different things
meltResultsList <- function(resultsList, id = NA){
  
  df <- reshape2::melt(resultsList)
  colnames(df) <- c("resolution", "neighbor", "Z", "reference")
  ## an identifier for the resultsList analysis
  df[["id"]] <- id
  return(df)
  
}