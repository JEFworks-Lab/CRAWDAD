#' Go from positions and celltype annotations to spatialpointsdataframe object
#' 
#' @param pos data frame; x and y coordinates of each cell
#' @param celltypes character vector; the cell type of each cell provided in pos
#' @param verbose Boolean; verbosity (default TRUE)
#' 
#' @export
toSP <- function(pos, celltypes, verbose=TRUE){
  
  if(length(levels(celltypes)) == 0){
    message("Warning: 'celltypes' does not have levels. Creating levels from values")
    celltypes <- factor(celltypes)
    names(celltypes) <- rownames(pos)
  }
  
  if(verbose){
    message("creating `sp::SpatialPointsDataFrame`")
  }
  
  cells <- sp::SpatialPointsDataFrame(
    coords = as.data.frame(pos),
    data = data.frame(
      celltypes = celltypes
      # name = rownames(pos)
    ))
  cells <- sf::st_as_sf(cells)
  
  ## Change rowname assignments of cells to integers.
  ## Solution to keep rows in same order later on when
  ## randomly shuffling cell labels
  rownames(cells) <- as.character(1:dim(cells)[1])
  
  # make asumption that cell type attribute is constant throughout the geometries of each cell
  ## it removed the warning that keep popping up, which says this assumption is made anyways
  sf::st_agr(cells) <- "constant"
  
  return(cells)
}


#' convert an sp::SpatialPointsDataFrame object to a dataframe with x y coords and cell type labels
#' @param cells sp::SpatialPointsDataFrame object, with celltypes features and point geometries
#'
#' @return datafame with columns: x, y, and celltype
#'
#' @export
spToDF <- function(cells){
  
  pos <- data.frame(sf::st_coordinates(sf::st_cast(cells$geometry,"POINT")))
  pos$type <- cells$celltypes
  colnames(pos) <- c("x", "y", "celltypes")
  return(pos)
  
}