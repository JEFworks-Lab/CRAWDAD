#' Calculate proportions of cell types inside the neighborhood
#' 
#' For a reference cell type and neighborhood distance, calculate the 
#' proportion of each neighbor cell type inside the neighborhood of the reference
#' cell type at that distance.
#' 
#' @param cells sf data.frame; as produced by crawdad::toSF function: cells with 
#' cell types annotated in the celltypes column and point positions in the 
#' geometry column
#' @param ref character; reference cell type to define neighborhood around
#' @param dist numeric vector; distances used to define the neighborhoods
#' 
#' @return named vector; the proportions
#' 
calculateCelltypeProportions <- function(cells, ref, dist) {
  ## create a circle around each reference cell 
  neighborhood <- sf::st_buffer(cells[cells$celltypes == ref,], dist) 
  ## merge the circles into a neighborhood (can take some time to compute)
  # neighborhood <- sf::st_union(buffer)
  ## calculate cells inside the neighborhood
  neighbor_cells <- sf::st_intersection(cells, neighborhood)
  ## remove duplicates
  self_cells <- cells[cells$celltypes == ref, ]
  neighbor_cells <- neighbor_cells[setdiff(rownames(neighbor_cells), 
                                           rownames(self_cells)), ]
  ## hack to accommodate self cells that are neighbors of another self cell
  neighbor_cells <- neighbor_cells[intersect(rownames(neighbor_cells), 
                                             c(rownames(cells), 
                                               paste0(rownames(self_cells), '.1'))), ]
  ## calculate proportions
  proportions <- as.vector(round(100*(table(neighbor_cells$celltypes)) / 
                                   (table(cells$celltypes)), 2))
  names(proportions) <- names(round(100*table(neighbor_cells$celltypes) / 
                                      table(cells$celltypes), 2))
  
  return(proportions)
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
