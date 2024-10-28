
#' get neighbor cells defined as being a distance away from a set of reference cells
#' @description get neighbor cells defined as being a distance away from a set of reference cells.
#'      `referenceIds` can be selected by subsetting rownames from `cells`:
#'       ex: `rownames(cells)[which(cells$celltypes == "A")]`
#'       or can be an entry in a subset list from `selectSubsets()`
#' 
#' @param cells sf object, with celltypes features and point geometries
#' @param referenceIds vector of cell ids (rownames) in `cells` to be used as the reference cell set
#' @param removeRef boolean; remove reference (default: TRUE)
#' @param dist distance to define neighbors (default = 50)
#' @param returnSF boolean to return either an sf object of just the neighbors
#' otherwise returns factor where non neighbor cells are NAs. (default: FALSE)
#'
#' @return sf object of the neighbor cells or factor of neighbor cells
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffleList <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000), ncores = 2)
#' binomMat <- binomialTestMatrix(cells, neighDist = 100, ncores = 2)
#' subsetList <- selectSubsets(binomMat, cells$celltypes, subType = "near", subThresh = 0.05)
#' neighCells <- getNeighbors(cells = cells, referenceIds = subsetList[["C_near_B"]],  dist = 100)
#' }
#' 
#' @export
getNeighbors <- function(cells,
                         referenceIds,
                         removeRef = TRUE,
                         dist = 50,
                         returnSF = FALSE){
  
  ## get the reference cells
  ref.cells <- cells[referenceIds,]
  
  ## define the buffer around the reference cells with distance of
  refs.buffer <- sf::st_buffer(ref.cells, dist)
  
  ## get the neighbor cells
  cells.inbuffer <- sf::st_intersection(cells, refs.buffer$geometry)
  ## remove duplicate neighbors
  neigh.cells <- cells.inbuffer[intersect(rownames(cells.inbuffer), rownames(cells)),]
  
  ## remove the reference cells from the set of neighbor cells
  if(removeRef){
    neigh.cells <- neigh.cells[setdiff(rownames(neigh.cells), referenceIds),]
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
#' shuffleList <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffleList = shuffleList, ncores = 2)
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
#' shuffleList <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffleList = shuffleList, ncores = 2)
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
#' shuffleList <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffleList = shuffleList, ncores = 2)
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
#' shuffleList <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffleList = shuffleList, ncores = 2)
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



#' Calculate AUC for each cell type pair from a dat file
#' 
#' @description
#' Calculate the AUC values for each cell type pair from the dat variable of a
#' sample.
#' 
#' @param datList list of `findTrends()` data.frames; the information about the 
#' scale, Z-score, reference and the neighbor cell. The input data.frame should 
#' be the results list from `findTrends()` that has been melted into a 
#' data.frame using `meltResultsList()`.
#' @param sharedPairs boolean; if true, will remove the AUC values from the 
#' cell-type pairs that are not available in all samples (default: TRUE).
#' 
#' @return data.frame; the AUC values for the cell-type pairs in the samples
#' 
#' @export
calculateAUC <- function(datList, sharedPairs = TRUE) {
  
  auc_list <- list()
  
  ## calculate
  for (dat in datList) {
    sample_id <- dat$id[1]
    
    ## define pairs
    pairs <- dplyr::distinct(dat[c("neighbor", "reference")])
    
    ## compute AUC for each pair
    auc_sample <- do.call(rbind, lapply(seq(nrow(pairs)), function(i) {
      ## process dat
      summarized_dat <- dat %>%
        dplyr::filter(reference == pairs$reference[i]) %>%
        dplyr::filter(neighbor == pairs$neighbor[i]) %>%
        dplyr::group_by(neighbor, scale, reference) %>%
        dplyr::summarise(mean = mean(Z),
                         sd = sd(Z))
      
      ## calculate auc
      return(data.frame(id = sample_id,
                        reference = pairs$reference[i], 
                        neighbor = pairs$neighbor[i], 
                        auc = pracma::trapz(summarized_dat$scale, summarized_dat$mean)))
    }))
    
    auc_list[[sample_id]] <- auc_sample
  }
  
  ## bind auc from all samples
  auc_df <- dplyr::bind_rows(auc_list) 
  
  ## process auc
  ## create ct pair colunm
  auc_df <- auc_df %>% 
    dplyr::mutate(pair = paste0(paste0(reference, ' - ', neighbor))) 
  
  ## if true, will remove cell type pairs that are not present in all samples
  if (sharedPairs) {
    ## calculate the ct pairs that are shared across all samples
    shared_pairs <- auc_df %>% 
      dplyr::group_by(pair) %>%
      dplyr::summarize(n_unique = dplyr::n_distinct(id)) %>% 
      dplyr::filter(n_unique == dplyr::n_distinct(auc_df$id)) %>% 
      dplyr::pull(pair)
    
    ## filter data based on shared pairs
    auc_df <- auc_df %>% 
      dplyr::filter(pair %in% shared_pairs)
  }
  
  return(auc_df)
}
