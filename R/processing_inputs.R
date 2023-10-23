#' Convert positions and cell type annotations to sf object
#' 
#' @param pos data frame; x and y coordinates of each cell
#' @param celltypes character vector; the cell type of each cell provided in pos
#' @param verbose Boolean; verbosity (default TRUE)
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(df)
#' }
#' 
#' @export
toSF <- function(pos, celltypes, verbose = TRUE){
  
  if(length(levels(celltypes)) == 0){
    message("Warning: 'celltypes' does not have levels. Creating levels from values")
    celltypes <- factor(celltypes)
    names(celltypes) <- rownames(pos)
  }
  
  df <- pos
  df$celltypes <- celltypes
  
  if(verbose){
    message("creating `sf` object")
  }
  
  cells_sf <- sf::st_as_sf(df, coords = c('x','y'))
  
  ## Change rownames assignments of cells to integers.
  ## Solution to keep rows in same order later on when
  ## randomly shuffling cell labels
  rownames(cells) <- as.character(1:dim(cells)[1])
  
  # make asumption that cell type attribute is constant throughout the geometries of each cell
  ## it removed the warning that keep popping up, which says this assumption is made anyways
  sf::st_agr(cells) <- "constant"
  
  return(cells)
}

#' Go from positions and celltype annotations to spatialpointsdataframe object
#' 
#' @param pos data frame; x and y coordinates of each cell
#' @param celltypes character vector; the cell type of each cell provided in pos
#' @param verbose Boolean; verbosity (default TRUE)
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSP(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' }
#' 
#' @export
toSP <- function(pos, celltypes, verbose=TRUE){
  
  .Deprecated("toSF")
  
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


#' Convert an sf object to a dataframe with x and y coordinates and cell type labels
#' 
#' @param cells sf object, with celltypes features and point geometries
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


#' Convert Seurat object to an sf object
#' 
#' @description Assumes that cell spatial coordinates and cell type annotations/communities are columns in the `meta.data` slot
#' 
#' @param obj the seurat object
#' @param coms name of `meta.dat`a column that contains cell annotations
#' @param posIDs columns in the `meta.data` that contain the spatial coordinates (default: c("x", "y")) 
#' @param verbose Boolean; verbosity (default TRUE)
#' 
#' @export
seuratToSF <- function(obj, coms, posIDs = c("x", "y"), verbose=TRUE){
  
  pos <- obj@meta.data[,posIDs]
  celltypes <- obj@meta.data[,coms]
  
  cells <- toSF(pos, celltypes)
  
  return(cells)
  
}


#' remove celltypes that have few cells in the dataset
#' 
#' @param cells sp::SpatialPointsDataFrame object, with celltypes features column
#' @param thresh threshold to remove cell types that contribute less than this fraction to the dataset (default: 0.0)
#' @param verbose Boolean; verbosity (default TRUE)
#' 
#' @noRd
filterCells <- function(cells, thresh = 0.0, verbose = TRUE){
  
  cellCounts <- table(cells$celltypes)
  typesBelow <- names(which(cellCounts/sum(cellCounts) < thresh))
  
  if(verbose){
    message("Removing celltypes whose fraction make up less than ", thresh, " of dataset:")
    message(cellCounts[typesBelow])
  }
  
  cellsFilt <- cells[!cells$celltypes %in% typesBelow,]
  
  return(cellsFilt)
  
}






