


# Exported ----------------------------------------------------------------



#' Visualize clusters in the same plot
#' 
#' @description Uses the cells sf object to visualize the clusters together in 
#' the same plot.
#' 
#' @param cells sf object; spatial (x and y) coordinates and celltypes column
#' @param ofInterest character vector; a vector of specific clusters to visualize
#' @param pointSize numeric; size of points
#' @param alpha numeric; transparency of points
#' @param ref character; reference cell type to draw the neighborhood 
#' around. If NULL, it will not create the neighborhood (default: NULL) 
#' @param dist numeric; distance to define neighbor cells with respect to each 
#' reference cell. If NULL, it will not create the neighborhood (default: NULL)
#' @param lineWidth numeric; width of neighborhood line
#' 
#' @return plot
#' 
#' @examples
#' \dontrun{
#' data(slide)
#' cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)
#' vizClusters(cells)
#' }
#' 
#' @export
vizClusters <- function(cells, ofInterest = NULL,
                        pointSize = 1, alpha = 0.5,
                        ref = NULL, dist = NULL, lineWidth = 0.1){
  
  ## if cells are a data.frame with "x" and "y" cell coordinate columns
  if( class(cells)[1] %in% c("data.frame", "matrix") ){
    stop('Use an sf object created by the crawdad::toSF function.')
  }
  
  ## define colors
  cluster_cols <- rainbow(n = length(unique(cells$celltypes)))
  names(cluster_cols) <- unique(cells$celltypes)
  cluster_cols['other'] <- '#E6E6E6'
  
  ## separate cells of interest to plot on top of others
  if(!is.null(ofInterest)){
    cells <- cells %>% 
      dplyr::mutate(celltypes = 
                      dplyr::case_when((!cells$celltypes %in% ofInterest) ~ 'other',
                                       T ~ celltypes))
    # cells$celltypes <- droplevels(cells$celltypes)
  }
  
  ## order cell types based on abundance
  ordered_cts <- c(names(sort(table(cells$celltypes), decreasing = T)))
  ordered_cts <- c('other', ordered_cts[ordered_cts != 'other'])
  cells <- cells %>% 
    dplyr::arrange(match(celltypes, ordered_cts))
  
  ## plot
  plt <- ggplot2::ggplot() +
    ## plot other cells
    ggplot2::geom_sf(data = cells, 
                     ggplot2::aes(color = celltypes), 
                     size = pointSize, alpha = alpha) +
    ## NA to gray
    ggplot2::scale_color_manual(values = cluster_cols, na.value = "#E6E6E6")
  plt
  
  if( (!is.null(ref)) & (!is.null(dist)) ) {
    ## create a circle around each reference cell 
    buffer <- sf::st_buffer(cells[cells$celltypes == ref,], dist) 
    ## merge the circles into a neighborhood (can take some time to compute)
    neighborhood <- sf::st_union(buffer)
    ## add to plot
    plt <- plt +
      ggplot2::geom_sf(data = neighborhood, fill = NA, 
                       color = 'black', linewidth = lineWidth)
  }
  
  ## add labels
  plt <- plt + 
    ggplot2::labs(x = "x",
                  y = "y") +
    ggplot2::theme_minimal() +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), 
                                                   ncol = 2))
  # ggplot2::coord_equal() ## geom_sf seems to be equal already
  
  return(plt)
  
}





#' Visualize each cluster separately
#' 
#' @description Returns a gridExtra of grobs.
#'     A single plot where each panel is a different cluster highlighted on the tissue
#' 
#' @param cells either a data.frame or sf object with cell spatial coordinates
#' @param coms a factor of cell type labels for the cells
#' @param axisAdj how much to increase axis ranges. If tissue, 100 okay, if embedding, 1 ok (default: 100)
#' @param size size of points (default: 0.01)
#' @param a alpha of points (default: 1; no transparency)
#' @param nacol color of the NA values for cells of "other" cluster (default: (transparentCol(color = "gray", percent = 50)))
#' @param clustcol color of the cells of the given cluster (default: red)
#' 
#' @return plot where each panel is a different celltype
#' 
#' @examples 
#' \dontrun{
#' data(slide)
#' vizEachCluster(slide, coms = slide$celltypes)
#' }
#' 
#' @export
vizEachCluster <- function(cells, coms, axisAdj = 1, s = 0.5, a = 1,
                           nacol = transparentCol(color = "gray", percent = 50),
                           clustcol = "red"){
  
  ## if cells are a data.frame or matrix  with "x" and "y" cell coordinate columns
  if( class(cells)[1] %in% c("data.frame", "matrix") ){
    pos <- cells[,c("x", "y")]
    ctemp <- factor(coms)
    names(ctemp) <- rownames(cells)
  }
  
  ## if cells are the sf object
  if( any(class(cells) == "sf") ){
    p <- sfToDF(cells)
    pos <- p[,c("x", "y")]
    ctemp <- factor(coms)
    names(ctemp) <- rownames(cells)
  }
  
  ps <- lapply(levels(ctemp), function(c) {
    
    tempCom <- ctemp
    tempCom[which(!tempCom %in% c(c))] <- NA
    tempCom <- droplevels(tempCom)
    
    cluster_cell_id <- which(tempCom == c)
    other_cells_id <- as.vector(which(is.na(tempCom)))
    
    dat <- data.frame("x" = pos[,"x"],
                      "y" = pos[,"y"])
    
    dat_cluster <- data.frame("x" = pos[cluster_cell_id,"x"],
                              "y" = pos[cluster_cell_id,"y"],
                              "Clusters" = c)
    
    dat_other <- data.frame("x" = pos[other_cells_id,"x"],
                            "y" = pos[other_cells_id,"y"],
                            "Clusters" = NA)
    
    plt <- ggplot2::ggplot() +
      
      ## plot other cells
      scattermore::geom_scattermore(data = dat_other, ggplot2::aes(x = x, y = y,
                                                                   color = Clusters), pointsize = s, alpha = a,
                                    pixels=c(1000,1000)) +
      ## cluster cells on top
      scattermore::geom_scattermore(data = dat_cluster, ggplot2::aes(x = x, y = y,
                                                                     color = Clusters), pointsize = s, alpha = a,
                                    pixels=c(1000,1000)) +
      
      ggplot2::scale_color_manual(values = c("red"), na.value = nacol) +
      
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$y)-axisAdj, max(dat$y)+axisAdj)) +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$x)-axisAdj, max(dat$x)+axisAdj) ) +
      
      ggplot2::labs(title = c, #paste0("c", c, ", volume normalized"),
                    x = "x",
                    y = "y") +
      
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, color = "black"),
                     axis.text.y = ggplot2::element_text(size=10, color = "black"),
                     axis.title.y = ggplot2::element_text(size=10),
                     axis.title.x = ggplot2::element_text(size=10),
                     axis.ticks.x = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size=10),
                     legend.text = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     legend.background = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     panel.grid.major.y =  ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(linewidth = 1, colour = "black"),
                     plot.margin = ggplot2::unit(c(1,1,1,1), "pt"), # change default margins around each plot
                     legend.position="none"
      ) +
      
      ggplot2::coord_equal()
    
    plt
    
  })
  
  numplots <- length(levels(ctemp))
  numcols <- 4
  numrows <- ceiling(numplots/numcols)
  numpanels <- numcols * numrows
  
  panel_layout <- t(matrix(data = 1:numpanels, nrow = numrows, ncol = numcols))
  
  p <- gridExtra::grid.arrange(
    grobs = ps,
    layout_matrix = panel_layout
  )
  
  p
}



#' Visualize grids and clusters
#' 
#' @description Uses the cells sf object and size of grid to visualize the grids 
#' used to create the null background.
#' 
#' @param cells sf object; spatial (x and y) coordinates and celltypes column
#' @param scale numeric; size of the scale to plot
#' @param permutation numeric; the number of the permutation of interest.
#' @param totalPermutations numeric; the total number of permutations used to 
#' shuffle the data.
#' @param square boolean; if true, create a squared grid, if false, make 
#' hexagonal grid (default TRUE)
#' @param ofInterest character vector; a vector of specific clusters to visualize
#' @param pointSize numeric; size of points
#' @param alpha numeric; transparency of points
#' 
#' @return plot
#' 
#' @examples
#' \dontrun{
#' data(slide)
#' cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)
#' vizAllClusters(cells)
#' }
#' 
#' @export
vizGrids <- function(cells, scale,
                     permutation = 1, totalPermutations = NULL, 
                     square = TRUE,
                     ofInterest = NULL, pointSize = 1, alpha = 0.5){
  
  ## the total number of permutations is used to calculate the offset, so if the
  ## user asks for a different permutation, but does not provide it, the code
  ## will not work
  if ((permutation != 1) & (is.null(totalPermutations))){
    stop('provide the total number of permutations')
  }
  
  ## create grids
  if (permutation == 1) {
    ## for no permutations 
    ## create grid with no offset
    grid <- sf::st_make_grid(cells, cellsize = scale, square = square)
  } else {
    ## for permutations 
    ## calculate offset
    offset <- -seq(from = 0, to = scale, by = scale/totalPermutations)[permutation]
    ## get bounding box with min and max coordinates
    bbox <- sf::st_bbox(cells$geometry)
    ## create grid with offset
    grid <- sf::st_make_grid(cells, cellsize = scale, 
                             offset = c(bbox[['xmin']] + offset, 
                                        bbox[['ymin']] + offset),
                             square = square)
  }
  
  
  ## get the coordinates of the centers of the grid tiles to add the tile IDs
  grid_coords_centroids <- as.data.frame(sf::st_coordinates(sf::st_centroid(grid)))
  grid_coords_centroids$name <- as.character(rownames(grid_coords_centroids))
  ## create plot
  crawdad::vizClusters(cells = cells, ofInterest, pointSize, alpha) + 
    ## add in the grid information on top of the plot
    ggplot2::geom_sf(data = grid, fill = NA) +
    ggplot2::geom_text(data = grid_coords_centroids, 
                       ggplot2::aes(X, Y, label = name))
}



#' Visualize grids and clusters after shuffling
#' 
#' @description Uses the cells sf object and size of grid to visualize the grids 
#' and the shuffled null background.
#' 
#' @param cells sf object; spatial (x and y) coordinates and celltypes column
#' @param shuffledList list; cell type labels shuffled at different scales 
#' (output from makeShuffledCells())
#' @param scale numeric; size of the scale to plot
#' @param permutation numeric; the number of the permutation of interest.
#' @param square boolean; if true, create a squared grid, if false, make 
#' hexagonal grid (default TRUE)
#' @param ofInterest character vector; a vector of specific clusters to visualize
#' @param pointSize numeric; size of points
#' @param alpha numeric; transparency of points
#' 
#' @return plot
#' 
#' @examples
#' \dontrun{
#' data(slide)
#' cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)
#' vizAllClusters(cells)
#' }
#' 
#' @export
vizShuffledGrids <- function(cells, shuffledList, scale, 
                             permutation = 1, 
                             square = TRUE,
                             ofInterest = NULL, pointSize = 1, alpha = 0.5){
  
  ## define colors to be consistent with the vizCluster function
  cluster_cols <- rainbow(n = length(unique(cells$celltypes)))
  names(cluster_cols) <- unique(cells$celltypes)
  cluster_cols['other'] <- '#E6E6E6'
  
  ## assing the shuffled labels to the celltypes column
  shuffled_labels <- shuffledList[[as.character(scale)]][[as.character(permutation)]]
  cells$celltypes <- as.factor(shuffled_labels)
  ## determine the total number of permutations to perform the offsetting
  totalPermutations <- length(shuffledList[[1]])
  
  ## use vizGrids to plot
  vizGrids(cells = cells, scale = scale, 
           permutation = permutation, totalPermutations = totalPermutations,
           square = square,
           ofInterest = ofInterest, pointSize = pointSize, alpha = alpha) +
    ## NA to gray
    ggplot2::scale_color_manual(values = cluster_cols, na.value = "#E6E6E6")
  
}



#' Plot Cell-type Proportions
#' 
#' For a chosen a neighborhood distance, plot a histogram of 
#' the proportion of each neighbor cell type inside the neighborhood of the 
#' reference cell type for all reference cell types.
#' 
#' @param cells sf data.frame; as produced by crawdad::toSF function: cells with 
#' cell types annotated in the celltypes column and point positions in the 
#' geometry column
#' @param neighDist numeric; distance used to define the neighborhood
#' 
#' @return ggplot2 plot; the a histogram of the proportions
#' 
#' @export
vizCelltypeProportions <- function(cells, neighDist) {
  
  ## for each cell type
  celltypes <- unique(cells$celltypes)
  props <- lapply(celltypes, calculateCelltypeProportions, 
                  cells = cells, neighDist = neighDist)
  df <- data.frame(proportions = unlist(props))
  
  df %>% ggplot2::ggplot(ggplot2::aes(x = proportions)) + 
    ggplot2::geom_histogram(color='#006437', fill='white', bins = 100) +
    ggplot2::theme_bw()
  
}



#' Plot trends with ggplot2
#' 
#' @description The input data.frame should be the results list from `findTrends()` that has been melted into a data.frame using `meltResultsList()`.
#' 
#' @param dat `findTrends()` results list, or data.frame; the information about the scale, Z-score, reference and the neighbor cell.
#' @param id column name that contains an additional feature to color trend lines (default: "id")
#' @param yaxis column that has significance value across scales (default: "Z")
#' @param zSigThresh threshold for significance, ie Z score significance threshold (default: 1.96).
#' @param nc number of colors to use for labeling features in id column
#' @param colors color assignment for each of the features in id column
#' @param title plot title (default: NULL)
#' @param facet boolean to facet wrap on reference and neighbor cell types. (default: TRUE)
#' @param lines boolean to plot lines (default: TRUE)
#' @param points boolean to plot points (default: TRUE)
#' @param withPerms if the results list is a list of lists using `returnMeans = FALSE` in `findTrends()`, then column order is different and this flag is needed (default: FALSE)
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffleList <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffleList = shuffleList, ncores = 2)
#' vizTrends(dat = results)
#' }
#' 
#' @export
vizTrends <- function(dat, id = "id", yaxis = "Z",
                      zSigThresh = 1.96, # -log10(0.05/nrow(dat)), ## sig thresh for num tests
                      nc = length(unique(dat[[id]])),
                      colors = rainbow(nc),
                      title = NULL,
                      facet = TRUE,
                      lines = TRUE,
                      points = TRUE,
                      withPerms = FALSE){
  
  ## if results is list from `findTrends()` then melt it.
  if (inherits(dat, "list")) {
    message("results detected to be a list. Melting to data.frame.")
    dat <- crawdad::meltResultsList(dat, withPerms = withPerms)
  }
  
  ## adds ref and neigh labels
  if (facet) {
    dat <- dat %>% 
      dplyr::mutate(reference = paste('reference:', reference),
             neighbor = paste('neighbor:', neighbor))
  }
  
  ## creates error bar
  if (withPerms) {
    dat <- dat %>% 
      dplyr::group_by(neighbor, scale, reference, id) %>% 
      dplyr::summarise(mean = mean(Z), 
                sd = sd(Z)) %>% 
      dplyr::mutate(Z = mean)
  }
  
  plt <- ggplot2::ggplot(data = dat)
  if(points){
    plt <- plt + ggplot2::geom_point(ggplot2::aes(x = scale, y = .data[[yaxis]], 
                                                  color = .data[[id]]), 
                                     size = 1.5)
  }
  if(lines){
    plt <- plt + ggplot2::geom_path(ggplot2::aes(x = scale, y = .data[[yaxis]], 
                                                 color = .data[[id]]), 
                                    size = 1)
  }
  plt <- plt + ggplot2::scale_color_manual(values = colors, na.value = "black") +
    ggplot2::geom_hline(yintercept = 0, color = "black", size = 1) +
    ggplot2::geom_hline(yintercept = zSigThresh, color = "black", size = 0.6, 
                        linetype = "dotted") +
    ggplot2::geom_hline(yintercept = -zSigThresh, color = "black", size = 0.6, 
                        linetype = "dotted") +
    # ggplot2::facet_grid(neighbor ~ reference) +
    ggplot2::labs(title = title) +
    # ggplot2::scale_x_log10() +
    # ggplot2::scale_y_continuous(trans = ggallin::pseudolog10_trans) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, color = "black", 
                                                       angle = -90, vjust = 0.5, 
                                                       hjust = 0),
                   axis.text.y = ggplot2::element_text(size=12, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   plot.background = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_line(linewidth = 0.1, 
                                                             colour = "black"),
                   panel.border = ggplot2::element_rect(colour = "black", fill=NA, 
                                                        linewidth=1),
                   axis.line = ggplot2::element_line(linewidth = 0, 
                                                     colour = "black"),
                   panel.spacing = ggplot2::unit(0.1, "lines"),
                   strip.text = ggplot2::element_text(size = 12),
                   legend.title = ggplot2::element_blank(),
                   # legend.position="none"
    )
  
  if (facet){
    plt <- plt + ggplot2::facet_grid(neighbor ~ reference, scales = 'free_y')
  }
  
  ## check if the data has permutations
  if (withPerms){
    # plt <- plt + ggplot2::geom_smooth(data = dat, 
    #                                   ggplot2::aes(x = scale, y = Z),
    #                                   se = TRUE) 
    pd <- position_dodge(0.1) # move them .05 to the left and right
    plt <- plt + ggplot2::geom_errorbar(ggplot2::aes(x=scale, 
                                                     ymin=mean-sd, ymax=mean+sd), 
                                        width=.1, position=pd, color="red")
  }
  
  plt
  
}


#' Plot trends as heatmaps with ggplot2
#' 
#' @description The input data.frame should be the results list from `findTrends()` that has been melted into a data.frame using `meltResultsList()`.
#' 
#' @param dat `findTrends()` results list, or data.frame; the information about the scale, Z-score, reference and the neighbor cell.
#' @param zSigThresh threshold for significance, ie Z score significance threshold (default: 1.96).
#' @param zScoreLimit Z score limits (default +/- 20)
#' @param palette_ color gradient for heatmap (default: grDevices::colorRampPalette(c("blue", "white", "red"))(n = 209))
#' @param title plot title (default: NULL)
#' @param withPerms if the results list is a list of lists using `returnMeans = FALSE` in `findTrends()`, then column order is different and this flag is needed (default: FALSE)
#' @param annotation boolean to show the Z score values in the squares of the heatmap (default: FALSE)
#' @param ncols specify number of columns in facet wrap (default: 4)
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffleList <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffleList = shuffleList, ncores = 2)
#' vizTrends.heatmap(dat = results)
#' }
#' 
#' @export
vizTrends.heatmap <- function(dat,
                      zSigThresh = 1.96, # -log10(0.05/nrow(dat)), ## sig thresh for num tests
                      zScoreLimit = 20,
                      palette_ = grDevices::colorRampPalette(c("blue", "white", "red"))(n = 209),
                      title = NULL,
                      withPerms = FALSE,
                      annotation = FALSE,
                      ncols = 4
                      ){
  
  ## if results is list from `findTrends()` then melt it.
  if(inherits(dat, "list")){
    message("results detected to be a list. Melting to data.frame.")
    dat <- crawdad::meltResultsList(dat, withPerms = withPerms)
  }
  
  ## save original for actual Z scores if annotation
  dat$scale <- factor(as.character(dat$scale), ordered = TRUE,
                           levels = as.character(sort(unique(dat$scale), decreasing = FALSE)))
  
  d <- dat
  ## winsorize high Z scores
  d$Z <- DescTools::Winsorize(d$Z, minval = -zScoreLimit, maxval = zScoreLimit)
  ## all non-significant scores to 0 so they are white on heat map
  d$Z[with(d, Z < zSigThresh & Z > -zSigThresh)] <- 0
  
  plt <- ggplot2::ggplot(data = d) +
    ggplot2::geom_tile(ggplot2::aes(x = scale, y = neighbor, fill=Z)) +
    ggplot2::facet_wrap(~reference, ncol = ncols) +
    ggplot2::scale_fill_gradientn(
      limits = c(-zScoreLimit, zScoreLimit),
      # values = breaks_,
      colors = palette_
    ) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, color = "black", 
                                                       angle = -90, vjust = 0.5, 
                                                       hjust = 0),
                   axis.text.y = ggplot2::element_text(size=10, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   # axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   plot.background = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_line(linewidth = 0.1, 
                                                             colour = "black"),
                   panel.border = ggplot2::element_rect(colour = "black", fill=NA, 
                                                        linewidth=1),
                   axis.line = ggplot2::element_line(linewidth = 0, 
                                                     colour = "black"),
                   panel.spacing = ggplot2::unit(0.1, "lines"),
                   strip.text = ggplot2::element_text(size = 8),
                   # legend.title = ggplot2::element_blank(),
                   # legend.position="none"
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Z score",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 1,
                                                   frame.colour= "black",
                                                   frame.linewidth = 1,
                                                   label.hjust = 0
    ))
  
  if(annotation){
    plt <- plt + ggplot2::geom_text(data = dat, ggplot2::aes(x = scale, y = neighbor, label = format(round(Z, 1), nsmall = 2) ), size = 2)
  }
  
  plt
  
}


#' Update cell type labels to only include specific cell types or subtypes
#' @description make a factor of selected cell type labels for visualizing specific cells with `vizAllClusters()`.
#'      Could append new entries into the `subsetList` to select and label other custom subsets of cells 
#'
#' @param df dataframe or sf object of cells
#' @param com original factor or vector of cell type labels for cells in `df`
#' @param cellIDs vector of cell type labels to include in output (default: NA)
#' @param subsetList list of subsets from `selectSubsets()` (default; NA)
#' @param subsetIDs vector of susbet cell type labels to include in output (names of vectors in `subsetList`) (default: NA)
#'
#' @return factor of specific cell type labels for visualizing with `vizAllClusters()` in the parameter: `clusters`
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffleList <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000), ncores = 2)
#' binomMat <- binomialTestMatrix(cells, neighDist = 100, ncores = 2)
#' subsetList <- selectSubsets(binomMat, cells$celltypes, subType = "near", subThresh = 0.05)
#' annots_temp <- selectLabels(df = cells, com = cells$celltypes, subsetList = subsetList, cellIDs = c("A", "B", "C", "D"), subsetIDs = c("C_near_B"))
#' }
#' 
#' @export
selectLabels <- function(df,
                         com,
                         cellIDs = NA,
                         subsetList = NA,
                         subsetIDs = NA
){
  
  ## get vector to append cell annotations of interest
  annots_temp <- rep(NA, length(rownames(df)))
  names(annots_temp) <- rownames(df)
  
  ## append the neighbor cells
  if(!is.na(cellIDs[1])){
    for(neigh_id in cellIDs){
      annots_temp[com == neigh_id] <- neigh_id
    }
  } else {
    message("`cellIDs` set to NA")
  }
  
  ## get cell ids that are part of subset
  ## these added after the neigbors in case
  ## you want to plot all CD4 T cells first
  # then color the ones that are a subset
  ## note that the order will be important
  ## because labels are overwritten in this way
  if( !is.na(subsetIDs[1]) & !is.na(subsetList[1]) ){
    for(subsetID in subsetIDs){
      cells_temp <- as.numeric(subsetList[[subsetID]])
      ## append the subset label
      annots_temp[cells_temp] <- subsetID
    }
  } else {
    message("`subsetIDs` and `subsetList` set to NA")
  }
  
  annots_temp <- as.factor(annots_temp)
  return(annots_temp)
}


#' Function to make transparent colors
#' 
#' @noRd
transparentCol <- function(color, percent = 50, name = NULL) {
  ## Get RGB values for named color
  rgb.val <- grDevices::col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- grDevices::rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                          maxColorValue = 255,
                          alpha = (100 - percent) * 255 / 100,
                          names = name)
  
  ## Save the color
  invisible(t.col)
}



#' Plot Summary of Relationships
#' 
#' @description This function takes the `findTrends()` melted data frame and 
#' plots the scale and the Z-score in which the trend first crossed the 
#' significance line (Z = 1.96). The Z-score was capped between -3 and 3. Since 
#' the relationships at smaller scales are more important than those at 
#' greater ones, we plotted the inverse of the scale so smaller ones would 
#' correspond to larger dots.
#' 
#' @param dat `findTrends()` data.frame; the information about the scale, 
#' Z-score, reference and the neighbor cell. The input data.frame should be the 
#' results list from `findTrends()` that has been melted into a data.frame 
#' using `meltResultsList()`.
#' @param zSigThresh numeric; the Z score significance threshold (default: 1.96).
#' @param pSigThresh numeric; the two-sided P value significance threshold. It 
#' can be used in place of the zSigThresh parameter. If no value is provided, 
#' the zSigThresh will be used.
#' @param zScoreLimit numeric; limit the Z-score to look better in the graph 
#' scale gradient. Z-score values above zScoreLimit will be represented as 
#' zScoreLimit, scores below -zScoreLimit will be represented as -zScoreLimit
#' (default: NULL).
#' @param reorder boolean; if TRUE, reorder the cell types by clustering on the 
#' z-score. If false, orders in alphabetical order (default: FALSE).
#' @param symmetrical boolean; highligh relationships that are symmetrical between cell 
#' type pairs (default: TRUE).
#' @param onlySignificant boolean; plot only cell types with significant 
#' relationships (default: FALSE).
#' @param colors character vector; colors for the gradient heatmap (low, mid, high) 
#' (default: c("blue", "white", "red")).
#' @param dotSizes numeric vector; minimum and maximum size of the dot 
#' (default: c(6,31)).
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffleList <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffleList = shuffleList, ncores = 2)
#' dat <- meltResultsList(results)
#' vizRelationships(dat)
#' }
#' 
#' @export
vizRelationships <- function(dat, zSigThresh = 1.96, pSigThresh = NULL,
                            zScoreLimit = NULL,  reorder = FALSE,
                            symmetrical = FALSE, onlySignificant = FALSE,
                            colors = c("blue", "white", "red"),
                            dotSizes = c(6,31)){
  
  if (!is.null(pSigThresh)) {
    zSigThresh <- round(qnorm(pSigThresh/2, lower.tail = F), 2)
  }
  
  ## create data.frame with the Z-scores and scales at the first scale
  ## the trend becomes significant
  ## get mean Z
  mean_dat <- dat %>% 
    dplyr::group_by(neighbor, scale, reference) %>% 
    dplyr::summarize(Z = mean(Z))
  ## get values before filtering
  max_scale <- max(dat$scale)
  u_cts <- unique(dat$reference)
  ## calculate sig z scores
  sig_dat <- mean_dat %>%
    dplyr::filter(abs(Z) >= zSigThresh) %>% 
    dplyr::group_by(neighbor, reference) %>% 
    dplyr::filter(scale == min(scale, na.rm = TRUE))
  
  ## limit the z-score for the gradient in the figure to look better
  if (!is.null(zScoreLimit)) {
    sig_dat$Z[sig_dat$Z > zScoreLimit] <- zScoreLimit
    sig_dat$Z[sig_dat$Z < -zScoreLimit] <- -zScoreLimit
  }
  
  ## scale sizes
  lsizes <- sort(unique(sig_dat$scale))
  legend_sizes <- c(lsizes[1],
                    round(mean(c(lsizes[1], lsizes[length(lsizes)]))),
                    lsizes[length(lsizes)])
  ## reorder
  if (reorder) {
    ## merge with all cts
    comb_cts <- tidyr::expand_grid(u_cts, u_cts)
    colnames(comb_cts) <- c('reference', 'neighbor')
    df_all <- dplyr::left_join(comb_cts, sig_dat, by = c('reference', 'neighbor'))
    df_all <- df_all %>% 
      dplyr::mutate(Z = dplyr::coalesce(Z, 0),
                    scale = dplyr::coalesce(scale, max_scale))
    ## create matrix
    sig_mat <- reshape::cast(df_all, neighbor~reference, value='Z')
    rownames(sig_mat) <- sig_mat[,1]
    sig_mat <- sig_mat[,-1]
    sig_mat[is.na(sig_mat)] <- 0
    ## cluster
    hc <- hclust(dist(sig_mat))
    sig_dat$neighbor <- factor(sig_dat$neighbor, 
                               levels=rownames(sig_mat)[hc$order])
    sig_dat$reference <- factor(sig_dat$reference, 
                                levels=colnames(sig_mat)[hc$order])
  }
  
  ## highlight symmetrical
  if (symmetrical) {
    # ## collect all pairs
    # ref_cts <- as.character(sig_dat$reference)
    # ngb_cts <- as.character(sig_dat$neighbor)
    # ## get Z scores and derive type of relationship to compare if it is the same
    # zscores <- sig_dat$Z
    df_pairs <- sig_dat %>% 
      dplyr::select(reference, neighbor, Z) %>% 
      ## calculate the type of relationship
      dplyr::mutate(type = dplyr::case_when(Z > 0 ~ 'enrichment',
                                            Z < 0 ~ 'depletion',
                                            T ~ NA)) %>% 
      dplyr::mutate(pair = paste(sort(c(gsub(" ", "", reference), 
                                        gsub(" ", "", neighbor))), 
                                 collapse = '_'))
    df_same_type <- df_pairs %>% 
      dplyr::group_by(pair) %>% 
      ## check if the type is not different for each ref of the pair and
      ## check if there are two relationships by checking distinct references
      dplyr::summarise(same_type = (dplyr::n_distinct(type) != 2) & 
                         (dplyr::n_distinct(reference) == 2))
    ## merge to reorder
    df_pairs <- df_pairs %>% 
      dplyr::left_join(df_same_type, by = 'pair')
    ## check if pairs are duplicate
    symmetrical_same_relationships <- df_pairs$same_type
    sig_dat$symmetrical <- symmetrical_same_relationships
  }
  
  ## plot figure
  p <- sig_dat %>% 
    ggplot2::ggplot(ggplot2::aes(x=reference, y=neighbor, 
                                 color=Z, size=scale)) +
    ggplot2::geom_point() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                       vjust = 0.5, 
                                                       hjust=1)) +
    {if (symmetrical) ggplot2::geom_point(data = ~dplyr::filter(.x, symmetrical == T),
                                   ggplot2::aes(x=reference, y=neighbor),
                                   shape = 18, color = 'gold', size = 2*dotSizes[1]/3)} + 
    ggplot2::scale_colour_gradient2(
      low = colors[1],
      mid = colors[2],
      high = colors[3],
      na.value = "lightgray"
    ) + 
    ggplot2::scale_radius(trans = 'reverse',
                          breaks = legend_sizes,
                          range = dotSizes) + 
    ggplot2::scale_x_discrete(position = "top") + 
    ggplot2::theme_bw()
  
  if (!onlySignificant) {
    all_cts <- sort(unique(dat$reference))
    alpha_cts <- sort(unique(dat$reference))
    if (reorder) {
      ct_order <- rownames(sig_mat)[hc$order]
      all_cts <- c(ct_order, setdiff(alpha_cts, ct_order))
    }
    p <- p +
      ggplot2::scale_x_discrete(limits = all_cts, position = 'top') +
      ggplot2::scale_y_discrete(limits = all_cts) 
  }
  
  return(p)
}



#' Visualize the samples using their AUC and PCA
#' 
#' @description
#' Uses the AUC values calculated from each cell-type pair in each sample to 
#' represent the samples in the PCA reduced dimension space.
#' 
#' @param aucSamples data.frame; the AUC values for each cell-type pair in each
#' sample as calculated by `calculateAUC()`
#' 
#' @return gglot2 plot; the PCA visualization of the samples
#' 
#' @export
#' 
vizPCASamples <- function(aucSamples) {
  sample_ids <- unique(aucSamples$id)
  
  ## create matrix
  auc_mtx <- aucSamples %>% 
    dplyr::select(c(pair, id, auc)) %>% 
    tidyr::pivot_wider(names_from = id, values_from = auc) %>% 
    dplyr::select(!pair) %>% 
    tidyr::drop_na() %>% ## drop nas
    as.matrix() %>% 
    t()
  
  ## normalize
  auc_mtx <- scale(auc_mtx)
  apply(auc_mtx, 2, mean)
  
  ## calculate pca
  pca <- prcomp(auc_mtx)
  pcs <- pca$x[, 1:2] %>% 
    as.data.frame() %>% 
    dplyr::mutate(id = rownames(pca$x)) %>% 
    dplyr::mutate(patient = 
                    sapply(rownames(pca$x), 
                           FUN = function(x) stringr::str_split(x, '_')[[1]][2]))
  
  ## plot
  pcs %>% 
    ggplot2::ggplot() + 
    ggplot2::geom_point(ggplot2::aes(x = PC1, y = PC2, color = id)) + 
    ggplot2::scale_color_manual(values = rainbow(length(sample_ids))) +
    ggplot2::theme_bw() +
    ggplot2::coord_equal()
}



#' Visualize the variance of each cell-type relationship in the samples
#' 
#' @description
#' Visualize the variance of the AUC values for each cell-type pair in all 
#' samples.
#' 
#' @param aucSamples data.frame; the AUC values for each cell-type pair in each
#' sample as calculated by `calculateAUC()`
#' 
#' @return gglot2 plot; the dot plot visualization of the variance in the 
#' samples
#' 
#' @export
#' 
vizVarianceSamples <- function(aucSamples) {
  
  ## all samples
  aucSamples %>% 
    dplyr::group_by(reference, neighbor) %>%
    dplyr::filter(reference != 'indistinct',
                  neighbor != 'indistinct') %>% 
    dplyr::summarize(variance = (var(auc))) %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = reference, y = neighbor, 
                                     size = variance), color = '#006437') +
    ggplot2::scale_radius(range = c(1, 10)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_discrete(position = 'top') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                       vjust = 0.5, 
                                                       hjust=0)) +
    ggplot2::coord_equal()
}



# To be deprecated --------------------------------------------------------

#' Plot Co-localization Dotplot
#' 
#' @description This function takes the `findTrends()` melted data frame and 
#' plots the scale and the Z-score in which the trend first crossed the 
#' significance line (Z = 1.96). The Z-score was capped between -3 and 3. Since 
#' the co-localization at smaller scales are more important than those at 
#' greater ones, we plotted the inverse of the scale so smaller ones would 
#' correspond to larger dots.
#' 
#' @param dat `findTrends()` data.frame; the information about the scale, 
#' Z-score, reference and the neighbor cell. The input data.frame should be the 
#' results list from `findTrends()` that has been melted into a data.frame 
#' using `meltResultsList()`.
#' @param zSigThresh numeric; the Z score significance threshold (default: 1.96).
#' @param pSigThresh numeric; the two-sided P value significance threshold. It 
#' can be used in place of the zSigThresh parameter. If no value is provided, 
#' the zSigThresh will be used.
#' @param zScoreLimit numeric; limit the Z-score to look better in the graph 
#' scale gradient. Z-score values above zScoreLimit will be represented as 
#' zScoreLimit, scores below -zScoreLimit will be represented as -zScoreLimit
#' (default: NULL).
#' @param reorder boolean; if TRUE, reorder the cell types by clustering on the 
#' z-score. If false, orders in alphabetical order (default: FALSE).
#' @param symmetrical boolean; highligh relationships that are symmetrical between cell 
#' type pairs (default: TRUE).
#' @param onlySignificant boolean; plot only cell types with significant 
#' relationships (default: FALSE).
#' @param colors character vector; colors for the gradient heatmap (low, mid, high) 
#' (default: c("blue", "white", "red")).
#' @param dotSizes numeric vector; minimum and maximum size of the dot 
#' (default: c(6,31)).
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffleList <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffleList = shuffleList, ncores = 2)
#' dat <- meltResultsList(results)
#' vizColocDotplot(dat)
#' }
#' 
#' @export
vizColocDotplot <- function(dat, zSigThresh = 1.96, pSigThresh = NULL,
                            zScoreLimit = NULL,  reorder = FALSE,
                            symmetrical = FALSE, onlySignificant = FALSE,
                            colors = c("blue", "white", "red"),
                            dotSizes = c(6,31)){
  
  ## deprecation warning
  message('this function will be deprecated, use vizRelationships instead')
  
  if (!is.null(pSigThresh)) {
    zSigThresh <- round(qnorm(pSigThresh/2, lower.tail = F), 2)
  }
  
  ## create data.frame with the Z-scores and scales at the first scale
  ## the trend becomes significant
  ## get mean Z
  mean_dat <- dat %>% 
    dplyr::group_by(neighbor, scale, reference) %>% 
    dplyr::summarize(Z = mean(Z))
  ## get values before filtering
  max_scale <- max(dat$scale)
  u_cts <- unique(dat$reference)
  ## calculate sig z scores
  sig_dat <- mean_dat %>%
    dplyr::filter(abs(Z) >= zSigThresh) %>% 
    dplyr::group_by(neighbor, reference) %>% 
    dplyr::filter(scale == min(scale, na.rm = TRUE))
  
  ## limit the z-score for the gradient in the figure to look better
  if (!is.null(zScoreLimit)) {
    sig_dat$Z[sig_dat$Z > zScoreLimit] <- zScoreLimit
    sig_dat$Z[sig_dat$Z < -zScoreLimit] <- -zScoreLimit
  }
  
  ## scale sizes
  lsizes <- sort(unique(sig_dat$scale))
  legend_sizes <- c(lsizes[1],
                    round(mean(c(lsizes[1], lsizes[length(lsizes)]))),
                    lsizes[length(lsizes)])
  ## reorder
  if (reorder) {
    ## merge with all cts
    comb_cts <- tidyr::expand_grid(u_cts, u_cts)
    colnames(comb_cts) <- c('reference', 'neighbor')
    df_all <- dplyr::left_join(comb_cts, sig_dat, by = c('reference', 'neighbor'))
    df_all <- df_all %>% 
      dplyr::mutate(Z = dplyr::coalesce(Z, 0),
                    scale = dplyr::coalesce(scale, max_scale))
    ## create matrix
    sig_mat <- reshape::cast(df_all, neighbor~reference, value='Z')
    rownames(sig_mat) <- sig_mat[,1]
    sig_mat <- sig_mat[,-1]
    sig_mat[is.na(sig_mat)] <- 0
    ## cluster
    hc <- hclust(dist(sig_mat))
    sig_dat$neighbor <- factor(sig_dat$neighbor, 
                               levels=rownames(sig_mat)[hc$order])
    sig_dat$reference <- factor(sig_dat$reference, 
                                levels=colnames(sig_mat)[hc$order])
  }
  
  ## highlight symmetrical
  if (symmetrical) {
    # ## collect all pairs
    # ref_cts <- as.character(sig_dat$reference)
    # ngb_cts <- as.character(sig_dat$neighbor)
    # ## get Z scores and derive type of relationship to compare if it is the same
    # zscores <- sig_dat$Z
    df_pairs <- sig_dat %>% 
      dplyr::select(reference, neighbor, Z) %>% 
      ## calculate the type of relationship
      dplyr::mutate(type = dplyr::case_when(Z > 0 ~ 'enrichment',
                                            Z < 0 ~ 'depletion',
                                            T ~ NA)) %>% 
      dplyr::mutate(pair = paste(sort(c(gsub(" ", "", reference), 
                                        gsub(" ", "", neighbor))), 
                                 collapse = '_'))
    df_same_type <- df_pairs %>% 
      dplyr::group_by(pair) %>% 
      ## check if the type is not different for each ref of the pair and
      ## check if there are two relationships by checking distinct references
      dplyr::summarise(same_type = (dplyr::n_distinct(type) != 2) & 
                         (dplyr::n_distinct(reference) == 2))
    ## merge to reorder
    df_pairs <- df_pairs %>% 
      dplyr::left_join(df_same_type, by = 'pair')
    ## check if pairs are duplicate
    symmetrical_same_relationships <- df_pairs$same_type
    sig_dat$symmetrical <- symmetrical_same_relationships
  }
  
  ## plot figure
  p <- sig_dat %>% 
    ggplot2::ggplot(ggplot2::aes(x=reference, y=neighbor, 
                                 color=Z, size=scale)) +
    ggplot2::geom_point() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                       vjust = 0.5, 
                                                       hjust=1)) +
    {if (symmetrical) ggplot2::geom_point(data = ~dplyr::filter(.x, symmetrical == T),
                                          ggplot2::aes(x=reference, y=neighbor),
                                          shape = 18, color = 'gold', size = 2*dotSizes[1]/3)} + 
    ggplot2::scale_colour_gradient2(
      low = colors[1],
      mid = colors[2],
      high = colors[3],
      na.value = "lightgray"
    ) + 
    ggplot2::scale_radius(trans = 'reverse',
                          breaks = legend_sizes,
                          range = dotSizes) + 
    ggplot2::scale_x_discrete(position = "top") + 
    ggplot2::theme_bw()
  
  if (!onlySignificant) {
    all_cts <- sort(unique(dat$reference))
    alpha_cts <- sort(unique(dat$reference))
    if (reorder) {
      ct_order <- rownames(sig_mat)[hc$order]
      all_cts <- c(ct_order, setdiff(alpha_cts, ct_order))
    }
    p <- p +
      ggplot2::scale_x_discrete(limits = all_cts, position = 'top') +
      ggplot2::scale_y_discrete(limits = all_cts) 
  }
  
  return(p)
}






# Deprecated --------------------------------------------------------------

#' Visualize all clusters on the tissue
#' 
#' @description uses the x and y position information and a chosen set of communities
#' 
#' @param cells either a data.frame or sf object with cell spatial coordinates
#' @param coms a factor of cell type labels for the cells
#' @param ofInterest a vector of specific clusters to visualize (default; NULL)
#' @param title title of plot (default: NULL)
#' @param axisAdj how much to increase axis ranges. If tissue, 100 okay, if embedding, 1 ok (default: 100)
#' @param size size of points (default: 0.01)
#' @param a alpha of points (default: 1; no transparency)
#' @param nacol color of the NA values for cells of "other" cluster (default: (transparentCol(color = "gray", percent = 50)))
#' 
#' @return plot
#' 
#' @examples
#' \dontrun{
#' data(slide)
#' vizAllClusters(slide, coms = slide$celltypes)
#' }
#' 
#' @export
vizAllClusters <- function(cells, coms, ofInterest = NULL,
                           axisAdj = 1, s = 0.5, a = 1, title = NULL,
                           nacol = transparentCol(color = "gray", percent = 50)){
  
  .Deprecated("vizClusters")
  
  ## if cells are a data.frame with "x" and "y" cell coordinate columns
  if( class(cells)[1] %in% c("data.frame", "matrix") ){
    pos <- cells[,c("x", "y")]
    tempCom <- factor(coms)
    names(tempCom) <- rownames(cells)
  }
  
  ## if cells are the sf object
  if( any(class(cells) == "sf") ){
    p <- spToDF(cells)
    pos <- p[,c("x", "y")]
    tempCom <- factor(coms)
    names(tempCom) <- rownames(cells)
  }
  
  if(!is.null(ofInterest)){
    ## goal:
    ## setup so the clusters of interest are plotted on top of everything else
    
    tempCom[which(!tempCom %in% ofInterest)] <- NA
    tempCom <- droplevels(tempCom)
    
    cluster_cell_id <- which(tempCom %in% ofInterest)
    other_cells_id <- as.vector(which(is.na(tempCom)))
    
    cluster_cols <- rainbow(n = length(ofInterest))
    names(cluster_cols) <- ofInterest
    
    dat <- data.frame("x" = pos[,"x"],
                      "y" = pos[,"y"])
    
    ## note: "Clusters" will be a variable id used to assign colors.
    ## for the "other cells" make this NA
    dat_cluster <- data.frame("x" = pos[cluster_cell_id,"x"],
                              "y" = pos[cluster_cell_id,"y"],
                              "Clusters" = as.vector(tempCom[cluster_cell_id]))
    
    dat_other <- data.frame("x" = pos[other_cells_id,"x"],
                            "y" = pos[other_cells_id,"y"],
                            "Clusters" = NA)
    
    plt <- ggplot2::ggplot() +
      
      ## plot other cells
      scattermore::geom_scattermore(data = dat_other, ggplot2::aes(x = x, y = y,
                                                                   color = Clusters), pointsize = s, alpha = a,
                                    pixels=c(1000,1000)) +
      ## cluster cells on top
      scattermore::geom_scattermore(data = dat_cluster, ggplot2::aes(x = x, y = y,
                                                                     color = Clusters), pointsize = s, alpha = a,
                                    pixels=c(1000,1000)) +
      
      ggplot2::scale_color_manual(values = cluster_cols, na.value = nacol)
    
  } else {
    
    tempCom <- droplevels(tempCom)
    dat <- data.frame("x" = pos[,"x"],
                      "y" = pos[,"y"],
                      "Clusters" = tempCom)
    
    plt <- ggplot2::ggplot(data = dat) +
      
      scattermore::geom_scattermore(ggplot2::aes(x = x, y = y,
                                                 color = Clusters), pointsize = s, alpha = a,
                                    pixels=c(1000,1000)) +
      
      ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom))), na.value = nacol)
  }
  
  plt <- plt + ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$y)-axisAdj, max(dat$y)+axisAdj)) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$x)-axisAdj, max(dat$x)+axisAdj) ) +
    
    ggplot2::labs(title = title,
                  x = "x",
                  y = "y") +
    
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                   axis.text.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 12, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(linewidth = 1, colour = "black")
                   # legend.position="none"
    ) +
    
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
    ) +
    
    ggplot2::coord_equal()
  
  plt
  
}



# Not exported ------------------------------------------------------------

#' Plot trends
#'
#' @description Plot panel of Z-score trends for each reference and neighbor cell-type pairs.
#'
#' @param results list or data.frame; the information about the scale, Z-score, reference and the neighbor cell. It can be the result directly obtained by the findTrends function or the melted version created by the `meltResultsList()` function.
#' @param idcol character; if results are a data.frame, this is the column that contains the additional feature to plot multiple trend lines with
#' @param legend boolean to produce legend, if results are a melted data.frame with "idcol" column (default: FALSE)
#' @param ... additional plotting parameters for base R plotting. Fed into "lines()" in script
#'
#' @return nothing
#'
#' @noRd
plotTrends <- function(results,
                       idcol = "id",
                       legend = FALSE,
                       ...){
  
  
  ## setup to check if original list output from `findTrends`, and plot one way
  ## or if a melted dataframe with ids from merging multiple trend analyses, check if dataframe and plot that way
  
  ## if in original list format from `findTrends`:
  if(inherits(results, "list")){
    message("results detected to be a list")
    
    par(mfrow=c(length(names(results)), length(names(results))),
        mar=rep(4,4))
    
    ## for each reference cell type, ie a dataframe in the list..
    invisible(sapply(names(results), function(ct1) {
      # print(ct1)
      results.norm <- results[[ct1]]
      results.norm[is.nan(results.norm)] <- NA
      results.norm[is.infinite(results.norm)] <- NA
      
      ## for each neighbor cell type...
      sapply(colnames(results.norm), function(ct2) {
        rg <- max(abs(results.norm[, ct2]), na.rm = TRUE)
        scales <- rownames(results.norm)
        
        if(is.infinite(rg)){
          rg <- 2.0
        }
        
        ## instantiate a plot and plot trend
        plot(scales, results.norm[,ct2],
             type="l", lwd = 2,
             main=paste0(ct1,' ref \n', ct2, " neighbors"),
             ylim=c(-rg, rg),
             xlab="scale", ylab="Z", ...)
        
        ## threshold lines
        abline(h = -2, col='red')
        abline(h = 2, col='red')
      })
    }))
    
    ## if a melted dataframe,
    ## will have an additional column that can serve to plot
    ## several trend lines on the same plot instance
    ## for example, ref vs neigh at different distances
  } else if(inherits(results, "data.frame")){
    message("results detected to be a data.frame")
    
    refs <- unique(results[,"reference"])
    neighs <- unique(results[,"neighbor"])
    ids <- unique(results[,idcol])
    
    cl <- rainbow(length(ids))
    
    par(mfrow=c(length(neighs), length(refs)),
        mar=c(4,4,4,6)) ## bot, top, left, right
    
    ## for each reference cell type...(rows)
    invisible(sapply(refs, function(ct1) {
      # print(ct1)
      results.norm <- results[results[,"reference"] == ct1,]
      results.norm[is.nan(results.norm[,"Z"]), "Z"] <- NA
      results.norm[is.infinite(results.norm[,"Z"]), "Z"] <- NA
      
      ## for each neighbor cell type...(columns)
      sapply(neighs, function(ct2) {
        results.norm.neigh <- results.norm[results.norm[,"neighbor"] == ct2,]
        
        yl <- max(abs(results.norm.neigh[, "Z"]), na.rm = TRUE)
        xl <- max(as.numeric(results.norm.neigh[,"scale"]))
        
        if(is.infinite(yl)){
          yl <- 2.0
        }
        
        ## instantiate a plot
        plot(0, 0, type = "n",
             main=paste0(ct1,' ref \n', ct2, " neighbors"),
             cex.main=1,
             ylim=c(-yl, yl),
             xlim=c(0, xl),
             xlab="scale", ylab="Z")
        
        ## for each id param, draw a line on plot instance
        for(i in 1:length(ids)){
          id <- ids[i]
          
          results.norm.neigh.id <- results.norm.neigh[results.norm.neigh[,idcol] == id,]
          
          lines(as.numeric(results.norm.neigh.id[,"scale"]), results.norm.neigh.id[,"Z"],
                type="l", lwd=2, col=cl[i], ...)
        }
        
        ## threshold lines
        abline(h = -2, col='red')
        abline(h = 2, col='red')
        
        if(legend){
          legend("topright", inset=c(-0.4,0), xpd=TRUE, legend = id, col=cl, pch=20, cex=0.5, title = idcol)
        }
        
      })
    }))
    
  } else {
    stop("`results` are neither a list from `findTrends()` or a melted data.frame from `meltResultsList()`")
  }
  
}


#' This one overlays each neighbor trend wrt the same reference cell type on the plot
#' 
#' @param results data.frame; the information about the scale, Z-score, reference and the neighbor cell. It can be the result directly obtained by the melted `findTrends` output created by the `meltResultsList()` function.
#' @param ... additional plotting parameters for base R plotting. Fed into "lines()" in script
#' 
#' @return nothing
#' 
#' @noRd
plotTrendsOverlay <- function(results,
                              ...){
  
  
  ## setup to check if original list output from `findTrends`, and plot one way
  ## or if a melted dataframe with ids from merging multiple trend analyses, check if dataframe and plot that way
  
  ## if in original list format from `findTrends`:
  if(inherits(results, "list")){
    message("results detected to be a list")
    
    par(mfrow=c(length(names(results)), length(names(results))),
        mar=rep(4,4))
    
    ## for each reference cell type, ie a dataframe in the list..
    sapply(names(results), function(ct1) {
      # print(ct1)
      results.norm <- results[[ct1]]
      results.norm[is.nan(results.norm)] <- NA
      results.norm[is.infinite(results.norm)] <- NA
      
      ## for each neighbor cell type...
      sapply(colnames(results.norm), function(ct2) {
        rg <- max(abs(results.norm[, ct2]), na.rm = TRUE)
        scales <- rownames(results.norm)
        
        ## instantiate a plot and plot trend
        plot(scales, results.norm[,ct2],
             type="l", lwd = 2,
             main=paste0(ct1,' ref \n', ct2, " neighbors"),
             ylim=c(-rg, rg),
             xlab="scale", ylab="Z", ...)
        
        ## threshold lines
        abline(h = -2, col='red')
        abline(h = 2, col='red')
      })
    })
    
    ## if a melted dataframe,
    ## will have an additional column that can serve to plot
    ## several trend lines on the same plot instance
    ## for example, ref vs neigh at different distances
  } else if(inherits(results, "data.frame")){
    message("results detected to be a data.frame")
    
    results <- results[,c("scale", "neighbor", "reference", "Z")]
    
    refs <- unique(results[,"reference"])
    neighs <- unique(results[,"neighbor"])
    cl <- rainbow(length(neighs))
    
    par(mfrow=c(length(refs),1),
        mar=c(4,4,4,8)) ## bot, top, left, right
    
    ## for each reference cell type...(rows)
    sapply(refs, function(ct1) {
      # print(ct1)
      results.norm <- results[results[,"reference"] == ct1,]
      results.norm[is.nan(results.norm[,"Z"]), "Z"] <- NA
      results.norm[is.infinite(results.norm[,"Z"]), "Z"] <- NA
      
      ## set limits based on trends with other cell types, not self
      results.norm.limits <- results.norm[results.norm[,"neighbor"] != ct1,]
      yl_max <- max(results.norm.limits[, "Z"], na.rm = TRUE)
      yl_min <- min(results.norm.limits[, "Z"], na.rm = TRUE)
      xl <- max(as.numeric(results.norm.limits[,"scale"]))
      if(is.infinite(yl_max)){
        yl_max <- 2.0
      }
      if(is.infinite(yl_min)){
        yl_min <- -2.0
      }
      
      ## instantiate a plot
      plot(0, 0, type = "n",
           main=paste0(ct1," ref"),
           cex.main=1,
           ylim=c(yl_min,yl_max),
           xlim=c(0, xl),
           xlab="scale", ylab="Z")
      
      ## for each neighbor cell type draw a line on plot instance
      for(i in 1:length(neighs)){
        ct2 <- neighs[i]
        # ignore showing trends with self because typically very large and
        # masks relationships with other cell types
        if(ct1 != ct2){
          results.norm.neigh.id <-  results.norm[results.norm[,"neighbor"] == ct2,]
          lines(as.numeric(results.norm.neigh.id[,"scale"]), results.norm.neigh.id[,"Z"],
                type="l", lwd=0.8, col=cl[i], ...)
        }
      }
      
      ## threshold lines
      abline(h = -1, col='black')
      abline(h = 1, col='black')
      
      legend("topright", inset=c(-0.4,0), xpd=TRUE, legend = neighs, col=cl, pch=20, cex=0.5, title = "neighbors")
      
    })
    
  } else {
    stop("`results` are neither a list from `findTrends` or a melted data.frame from `meltResultsList`")
  }
  
}


#' Create color for the "other" cell type
#'
#' @description Given a color and a transparency in percentage, create a color.
#'
#' @param color character; color
#' @param percent numeric (0 to 100); transparency
#' @param name character; name of the color
#' 
#' @return color
#'
#' @noRd
transparentCol <- function(color, percent = 50, name = NULL) {
  ## Get RGB values for named color
  rgb.val <- grDevices::col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- grDevices::rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                          maxColorValue = 255,
                          alpha = (100 - percent) * 255 / 100,
                          names = name)
  
  ## Save the color
  invisible(t.col)
}