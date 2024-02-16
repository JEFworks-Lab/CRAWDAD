#' Shuffle the celltype labels at different shuffling scales
#'
#' @description for each scale, and for i number of permutations, shuffle the celltype labels.
#' Returns a list of lists, where each list is a factor of the shuffled celltype labels.
#' Outer list is for each scale of shuffling, and each inner list is for a given permutation with a given seed.
#' @param cells sf object, with celltypes features and point geometries
#' @param scales numeric vector of the different scales to shuffle at and subsequently compute significance at (default: c(50, 100, 200, 300, 400, 500))
#' @param perms number of permutations to shuffle for each scale (default = 1)
#' @param ncores number of cores for parallelization (default 1)
#' @param seed set seed for shuffling (if more than 1 permutation, each shuffling permutation has a seed equal to the permutation number)
#' @param square if false, make hexagonal grid (default TRUE)
#' @param verbose Boolean for verbosity (default TRUE)
#' 
#' @examples  
#' \dontrun{
#' data(sim)
#' makeShuffledCells(sim, scales = c(50, 100, 200, 300, 400, 500), ncores = 2)
#' }
#'
#' @export
makeShuffledCells <- function(cells,
                              scales = c(50, 100, 200, 300, 400, 500),
                              perms = 1,
                              ncores = 1,
                              seed = 0,
                              square = TRUE,
                              verbose = TRUE){
  
  ## effectively will make a list of lists of randomly shuffled cell labels.
  ## a list for each scale that contains factors of shuffled cell labels for each permutation
  ## use the cell labels to reorder the labels of the `cells`
  
  ## check if cells is an `sf object`
  if( !any(class(cells) == "sf") ){
    stop("`cells` needs to be an `sf` object. You can make this using `toSF()`")
  }
  
  if( !any(grepl("celltypes", colnames(cells))) ){
    stop("`cells` needs a column named `celltypes`. You can make this using `toSF()`")
  }
  
  if(verbose){
    start_time <- Sys.time()
  }
  
  cells_df <- sfToDF(cells)
  
  randomcellslist <- lapply(scales, function(r) {
    
    ## create list of offsets for the permutations
    offsets <- -seq(from = 0, to = r, by = r/perms)
    
    permutations <- lapply(1:perms, function(i){
      
      ## if only 1 permutation, use the set seed. Otherwise seed is equal to permutation number
      if (perms == 1){
        s <- seed
      } else if (perms > 1){
        s <- i
      }
      
      if(verbose){
        message("shuffling permutation ", i, " using seed ", s)
      }
      
      ## create grid after going into permutations
      grid <- sf::st_make_grid(cells, cellsize = r, 
                               offset = c(min(cells_df$x) + offsets[i], min(cells_df$y) + offsets[i]),
                               square = square)
      
      if(verbose){
        message(r, " unit scale")
        message(length(grid), " tiles to shuffle...")
      }
      
      ## disable scientific notation.
      ## apparently, when I sort the names by converting to numeric then back to characters
      ## some numerics are written in scientific notation but then the notation and not
      ## the exact number is converted back to character.
      ## for example: 1e+5 comes back as "1e+5".
      ## this effectively results in NA values and loss of cells
      ## subsequently there is a discrepancy between the number of cells in the shuffled list
      ## vs the cells of the actual data.
      options(scipen = 999)
      
      ## shuffle within grid once
      randcelltype <- unlist(BiocParallel::bplapply(1:length(grid), function(i) {
        ## sometimes can be on boundary, so just pick first one
        int <- sf::st_intersection(cells, grid[[i]])
        
        # randomly grab cell labels for cell in the grid (this is a factor)
        set.seed(s)
        shuffled_cells <- sample(int$celltypes)
        # assign the cell ids to the randomly sampled cell labels
        
        ## convert to named character vector
        ## the factor levels get lost later on
        ## when trying to combine the removed dups
        ## and selected dups
        shuffled_cells <- as.character(shuffled_cells)
        names(shuffled_cells) <- rownames(int)
        return(shuffled_cells)
      }, BPPARAM=BiocParallel::SnowParam(workers=ncores)))
      
      
      
      # reorder so cells are "1", "2", etc. and not in the order based on grids
      randcelltype <- randcelltype[as.character(sort(as.numeric(names(randcelltype))))]
      #print(length(randcelltype))
      
      ## apparently, some cells can end up in multiple grids thus leading to double counting
      ## this is problematic because later on, the celltype names in the original cell dataframe
      ## are replaced with the shuffled values, but if different length, will cause an error
      ## so need to deal with cells that are counted more than once
      
      ## find the duplicates in each case and pick one of them to use
      n_occur <- data.frame(table(names(randcelltype)))
      dups <- as.character(n_occur[n_occur$Freq > 1,]$Var1)
      #print(length(dups))
      removedDups <- randcelltype[!names(randcelltype) %in% dups]
      #print(length(removedDups))
      
      ## deal with duplicates by picking the first one
      selectedDups <- unlist(lapply(dups, function(dup){
        randcelltype[names(randcelltype) == dup][1]
      }))
      #print(length(selectedDups))
      
      ## recombine the unique cells and the selected duplicate cells
      randcelltype2 <- c(removedDups, selectedDups)
      
      ## reorder 1,2,3, etc like done previously
      randcelltype2 <- randcelltype2[as.character(sort(as.numeric(names(randcelltype2))))]
      #print(length(randcelltype2))
      randcelltype2
    })
    
    names(permutations) <- 1:perms
    return(permutations)
  })
  
  if(verbose){
    total_t <- round(difftime(Sys.time(), start_time, units="mins"), 2)
    message(sprintf("Time was %s mins", total_t))
  }
  
  names(randomcellslist) <- scales
  return(randomcellslist)
}


#' compute significant different between real and randomly shuffled cell neighbor proportions
#' 
#' @description for a given reference cell type, computes the Z scores for being colocalized or separated from each query cell type. 
#' If `returnMeans = TRUE`, then the result will be a data.frame where each row is a scale, each column is a query cell type, and each value is the Z score.
#' If `returnMeans = FALSE`, then the result will be a list of data.frames, where each data.frame is a scale,
#' the rows are permutations, the columns are query cell types, and each value is a Z score.
#'
#' @param cells sf object of all the cells
#' @param randomcellslist list of lists of randomly shuffled cell type labels produced from `makeShuffledCells`
#' @param trueNeighCells Simple feature collection of real cells for a given reference cell type, with geometries of a given dist (from sf::st_buffer)
#' @param cellBuffer Simple feature collection of the neighbor cells that are within "dist" of the ref cells (from sf::intersection)
#' @param ncores number of cores for parallelization (default 1)
#' @param removeDups remove duplicate neighbor cells to prevent them from being counted multiple times and inflate the Z scores (default: TRUE)
#' @param returnMeans if multiple permutations, return the mean Z score across the permutations in each scale with respect to each neighbor cell type (default: TRUE) 
#'
#'@export
evaluateSignificance <- function(cells,
                                 randomcellslist,
                                 trueNeighCells,
                                 cellBuffer,
                                 ncores = 1,
                                 removeDups = TRUE,
                                 returnMeans = TRUE){
  
  allcells <- cells
  trueNeighCells <- trueNeighCells
  cellBuffer <- cellBuffer
  
  ## if true, will return a data.frame
  if(returnMeans){
    
    ## for each scale:
    results <- do.call(rbind, BiocParallel::bplapply(randomcellslist, function(cellsAtRes){
      
      ## Iterate through each permutation of a given scale and 
      ## produce the scores for each neighbor cell type. 
      ## Scores for a given permutation added as a row to a dataframe.
      ## Take the column mean of the scores for each neighbor cell type across permutations.
      ## The resulting vector of score means is returned and appended as a row in the `results` data,frame,
      ## where each row is a scale, and contains Z scores for each neighbor cell type (the columns)
      
      scores <- do.call(rbind, lapply(cellsAtRes, function(randomcellslabels){
        
        randomcells <- allcells
        randomcells$celltypes <- as.factor(randomcellslabels)
        sf::st_agr(randomcells) <- "constant"
        
        bufferrandomcells <- sf::st_intersection(randomcells, cellBuffer$geometry)
        
        ## remove duplicate neighbor cells to prevent them from being counted multiple times
        ## and inflate the Z scores
        if(removeDups){
          # message("number of permuted neighbor cells before: ", nrow(bufferrandomcells))
          bufferrandomcells <- bufferrandomcells[intersect(rownames(bufferrandomcells), rownames(randomcells)),]
          # message("number of permuted neighbor cells after removing dups: ", nrow(bufferrandomcells))
        }
        
        ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
        y1 <- table(trueNeighCells$celltypes)
        y2 <- table(bufferrandomcells$celltypes)
        n1 <- length(trueNeighCells$celltypes)
        n2 <- length(bufferrandomcells$celltypes)
        p1 <- y1/n1
        p2 <- y2/n2
        p <- (y1+y2)/(n1+n2)
        Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
        
        rm(bufferrandomcells)
        rm(randomcells)
        gc(verbose = FALSE, reset = TRUE)
        
        return(Z)
      }))
      
      ## returning the mean Z score across permutations for the given scale
      return(colMeans(scores))
      
    }, BPPARAM=BiocParallel::SnowParam(workers=ncores)))
    
    ## otherwise, returns a list
  } else {
    
    ## for each scale:
    results <- BiocParallel::bplapply(randomcellslist, function(cellsAtRes){
      
      ## Iterate through each permutation of a given scale and 
      ## produce the scores for each neighbor cell type. 
      ## Scores for a given permutation added as a row to a dataframe.
      ## Take the column mean of the scores for each neighbor cell type across permutations.
      ## The resulting vector of score means is returned and appended as a row in the `results` data,frame,
      ## where each row is a scale, and contains Z scores for each neighbor cell type (the columns)
      
      scores <- do.call(rbind, lapply(cellsAtRes, function(randomcellslabels){
        
        randomcells <- allcells
        randomcells$celltypes <- as.factor(randomcellslabels)
        sf::st_agr(randomcells) <- "constant"
        
        bufferrandomcells <- sf::st_intersection(randomcells, cellBuffer$geometry)
        
        ## remove duplicate neighbor cells to prevent them from being counted multiple times
        ## and inflate the Z scores
        if(removeDups){
          # message("number of permuted neighbor cells before: ", nrow(bufferrandomcells))
          bufferrandomcells <- bufferrandomcells[intersect(rownames(bufferrandomcells), rownames(randomcells)),]
          # message("number of permuted neighbor cells after removing dups: ", nrow(bufferrandomcells))
        }
        
        ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
        y1 <- table(trueNeighCells$celltypes)
        y2 <- table(bufferrandomcells$celltypes)
        n1 <- length(trueNeighCells$celltypes)
        n2 <- length(bufferrandomcells$celltypes)
        p1 <- y1/n1
        p2 <- y2/n2
        p <- (y1+y2)/(n1+n2)
        Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
        
        rm(bufferrandomcells)
        rm(randomcells)
        gc(verbose = FALSE, reset = TRUE)
        
        return(Z)
      }))
      
      ## returning the data.frame of Z scores for each permutation (row)
      return(scores)
      
    }, BPPARAM=BiocParallel::SnowParam(workers=ncores))
    
    ## will be list of data.frame in this case
    ## each data.frame is a scale
    names(results) <- names(randomcellslist)
    
  }
  
  rm(allcells)
  rm(trueNeighCells)
  rm(cellBuffer)
  gc(verbose = FALSE, reset = TRUE)
  
  return(results)
}


#' Generate matrix of pvalues indicating if a cell is enriched in neighbors of a given cell type
#' @description pvalues are based on a binomial test and neighbors are defined within a given distance from a cell
#' @param cells sf object, with celltypes features and point geometries
#' @param neigh.dist distance to define neighbors (default = 100)
#' @param ncores number of cores for parallelization (default 1)
#' @param verbose Boolean for verbosity (default TRUE)
#' 
#' @return matrix where rows are cells, columns are cell types and values are p-values whether or not a cell is enriched in neighbors of a given cell type based on a binomial test.
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffle.list <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000, 1500, 2000), ncores = 2)
#' binomMat <- binomialTestMatrix(cells, neigh.dist = 100, ncores = 2)
#' }
#' 
#' @export
binomialTestMatrix <- function(cells,
                               neigh.dist = 100,
                               ncores = 1,
                               verbose = TRUE) {
  
  start_time <- Sys.time()
  
  ## check to make sure cells is a sf object
  if( !any(class(cells) == "sf") ){
    stop("`cells` needs to be an `sf` object. You can make this using `toSF()`")
  }
  if( !any(grepl("celltypes", colnames(cells))) ){
    stop("`cells` needs a column named `celltypes`. You can make this using `toSF()`")
  }
  
  if(verbose){
    message("Binomial test for each cell testing if it is enriched in neighbors of a given cell type based on distance of ",
            neigh.dist)
  }
  
  cell.types <- cells$celltypes
  names(cell.types) <- rownames(cells)
  
  if(length(levels(cell.types)) == 0){
    message("Warning: `celltypes` does not have levels. Creating levels from values")
    cell.types <- factor(cell.types)
    names(cell.types) <- rownames(cells)
  }
  
  ## get global fractions of each cell type (hypothesized probability of success)
  p <- table(cell.types)/sum(table(cell.types))
  
  ## get buffer around each cell with diameter of neigh.dist
  refs.buffer <- sf::st_buffer(cells, neigh.dist)
  
  ## define the binomial test to be performed and used in the mapply below
  binom <- function(x, n, p){stats::binom.test(x, n, p, alternative="greater")$p.value}
  
  ## for each cell: (parallelize this)
  ## get number of each cell type that are neighbors (x, ie number of successes)
  ## and total number of neighbors (n, number of trials)
  ## compute pvals for a binomial test for each cell type neighbor
  ## As we loop through each cell, add the pvals for each cell type as a new row of a matrix
  message("Performing tests...")
  results <- do.call(rbind, BiocParallel::bplapply(rownames(cells), function(c){
    
    cells.inbuffer <- sf::st_intersection(cells, refs.buffer[c,]$geometry)
    x <- table(cells.inbuffer$celltypes)
    n <- rep(sum(x), length(x))
    
    ## apply the binomial test across the vectors using mapply
    pvals <- mapply(binom, x, n, p)
    pvals
    
  }, BPPARAM=BiocParallel::SnowParam(workers=ncores)))
  
  rownames(results) <- rownames(cells)
  colnames(results) <- names(p)
  
  if(verbose){
    total_t <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
    message(sprintf("Time to compute was %smins", total_t))
  }
  
  return(results)
}


#' find subsets of cells
#' @description assign cells to subsets based on whether they are enriched in a given neighbor cell type based on the pvalues in the input `binomMatrix`
#'
#' @param binomMatrix matrix where rows are cells, columns are cell types and values are p-values whether or not a cell is enriched in neighbors of a given cell type based on a binomial test. Output from `binomialTestMatrix()`
#' @param celltypes named vector or factor of cell type labels of each cell in `binomMatrix` and in the same order.
#' @param sub.type subset type, either ref cells "near" (ie localized) a neighbor cell type, or "away" (ie separated) from a neighbor cell type.
#' @param sub.thresh significance threshold for the binomial test (default = 0.05)
#' @param ncores number of cores for parallelization (default 1)
#' @param verbose Boolean for verbosity (default TRUE)
#' 
#' @return list where each entry is a subset and the values are the cell ids determined to be in the subset
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffle.list <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000, 1500, 2000), ncores = 2)
#' binomMat <- binomialTestMatrix(cells, neigh.dist = 100, ncores = 2)
#' subset.list <- selectSubsets(binomMat, cells$celltypes, sub.type = "near", sub.thresh = 0.05) 
#' }
#' 
#' @export
selectSubsets <- function(binomMatrix,
                          celltypes,
                          sub.type = c("near", "away"),
                          sub.thresh = 0.05,
                          ncores = 1,
                          verbose = TRUE){
  
  start_time <- Sys.time()
  sub.type <- match.arg(sub.type)
  
  if(!sub.type %in% c("near", "away")){
    stop("`sub.type` must be either 'near' or 'away'")
  }
  
  ## make sure cell types is a factor
  cell.types <- factor(celltypes)
  
  ## setup possible subset combinations:
  combos <- expand.grid(rep(list(1:length(levels(cell.types))),2))
  colnames(combos) <- c("ref", "neighbors")
  combo_ids <- unlist(lapply(rownames(combos), function(i){
    cells.ref.ix <- levels(cell.types)[as.numeric(combos[i,1])]
    cells.neighbors.ix <- levels(cell.types)[as.numeric(combos[i,2])]
    id <- paste0(cells.ref.ix, "_", sub.type, "_", cells.neighbors.ix)
    id
  }))
  
  binomMat <- binomMatrix
  
  subsets <- BiocParallel::bplapply(rownames(combos), function(i){
    
    cells.ref.ix <- levels(cell.types)[as.numeric(combos[i,1])]
    cells.neighbors.ix <- levels(cell.types)[as.numeric(combos[i,2])]
    id <- paste0(cells.ref.ix, "_", sub.type, "_", cells.neighbors.ix)
    message("computing subsets for ", id)
    
    ## get reference cell rows
    ref.cells <- binomMat[which(cell.types == cells.ref.ix),]
    
    ## subset reference cells whose pval is below thresh for given neighbor cell type
    ## return cells that were significant
    if(sub.type == "near"){
      sub.cells <- rownames(ref.cells)[which(ref.cells[,cells.neighbors.ix] < sub.thresh)]
    }
    
    ## return the cells that were not significant
    if (sub.type == "away"){
      ## recommend setting threshold very liberal, like 0.5,
      ## that way, only the cells that couldn't even pass a p-val cutoff of 0.5 would be selected for,
      ## and these would be expected to be very much depleted or separated from the neighbor cell type
      sub.cells <- rownames(ref.cells)[which(ref.cells[,cells.neighbors.ix] < sub.thresh)]
      sub.cells <- rownames(ref.cells)[which(!rownames(ref.cells) %in% sub.cells)]
    }
    
    return(sub.cells)
    
  }, BPPARAM=BiocParallel::SnowParam(workers=ncores))
  
  if(verbose){
    total_t <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
    message(sprintf("Time to compute was %smins", total_t))
  }
  
  names(subsets) <- combo_ids
  return(subsets)
  
}


#' Compute trends of cell type colocalization for each cell type combination across specified scales
#'
#' @description Trends are based on significant differences in cell type proportions between the real and randomly shuffled datasets.
#' Cell type proportions are with respect to the different cell types that are neighboring the cells of a given reference cell type within a certain defined distance.
#' This is done at difference scales, where a scale is whether the cell type labels are shuffled locally or globally.
#' Trends are essentially built from significance values. The significance test basically asks if two cell types are localized or separated by assessing if the proportion of the neighboring cell type is significantly greater, or less than, random chance.
#'
#' @param cells sf object, with celltypes features and point geometries
#' @param dist numeric distance to define neighbor cells with respect to each reference cell (default: 100)
#' @param ncores number of cores for parallelization (default 1)
#' @param shuffle.list a list of cell type labels shuffled at different scales (output from `makeShuffledCells()`)
#' @param subset.list a subset list (output from `selectSubsets()`). Required if computing trends for subsets (default NULL)
#' @param verbose Boolean for verbosity (default TRUE)
#' @param removeDups remove duplicate neighbor cells to prevent them from being counted multiple times and inflate the Z scores (default: TRUE)
#' @param returnMeans if multiple permutations, return the mean Z score across the permutations in each scale with respect to each neighbor cell type (default: TRUE)
#'
#' @return A list that contains a dataframe for each reference cell type, where the dataframe contains the significance values for each neighbor cell type at each scale
#' 
#' @examples
#' \dontrun{
#' data(sim)
#' shuffle.list <- makeShuffledCells(sim, scales = c(50, 100, 200, 300, 400, 500))
#' findTrends(sim, dist = 100, shuffle.list = shuffle.list, ncores = 2)
#' }
#' 
#' @export
findTrends <- function(cells,
                       dist = 100,
                       ncores = 1,
                       shuffle.list,
                       subset.list = NULL,
                       verbose = TRUE,
                       removeDups = TRUE,
                       returnMeans = TRUE){
  
  if(!is.list(shuffle.list)){
    stop("`shuffle.list` is not a list. You can make this using `makeShuffledCells()`")
  }
  
  if( !any(class(cells) == "sf") ){
    stop("`cells` needs to be an `sf` object. You can make this using `toSF()`")
  }
  
  if( !any(grepl("celltypes", colnames(cells))) ){
    stop("`cells` needs a column named `celltypes`. You can make this using `toSF()`")
  }
  
  if(verbose){
    start_time <- Sys.time()
  }
  
  sf::st_agr(cells) <- "constant"
  
  if(verbose){
    message("Evaluating significance for each cell type")
    message("using neighbor distance of ", dist)
  }
  
  ## Evaluate significance (pairwise)
  if(is.null(subset.list)){
    
    if(verbose){
      message("Calculating for pairwise combinations")
    }
    
    celltypes <- factor(cells$celltypes)
    
    d <- dist
    results.all <- lapply(levels(celltypes), function(ct) {
      
      if(verbose){
        message(ct)
      }
      
      # get polygon geometries of reference cells of "celltype" up to defined distance "dist"
      # use this to assess neighbors within "d" um of each cell
      ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
      ## union polygons to avoid memory overflow, too slow
      # ref.buffer_union <- sf::st_union(ref.buffer)
      # get the different types of neighbor cells that are within "d" of the ref cells
      neigh.cells <- sf::st_intersection(cells, ref.buffer$geometry)
      
      ## remove duplicate neighbor cells to prevent them from being counted multiple times
      ## and inflate the Z scores
      if(removeDups){
        ## need to remove self cells else have trivial enrichment when d~0
        self.cells <- cells[cells$celltypes == ct,]
        neigh.cells <- neigh.cells[setdiff(rownames(neigh.cells), rownames(self.cells)),]
        
        # message("number of neighbor cells before: ", nrow(neigh.cells))
        #neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), rownames(cells)),]
        ## hack to accommodate self cells that are neighbors of another self cell
        neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), c(rownames(cells), paste0(rownames(self.cells), '.1'))),]
        # message("number of neighbor cells after removing dups: ", nrow(neigh.cells))
      }
      
      ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
      ## chose to shuffle the scales in parallel, but in each scale, the perms done linearly
      ## I think a bottle neck originally was waiting for certain cell types to finish
      ## so if I split up the scales, might be able to get through each cell type faster and speed up entire process?
      ## I could also split up the permutations, but then each cell type for each scale is done one by one
      ## I could do cell types in parallel, but then for each cell type need to go through each res and each perm one by one
      results <- evaluateSignificance(cells = cells,
                                      randomcellslist = shuffle.list,
                                      trueNeighCells = neigh.cells,
                                      cellBuffer = ref.buffer,
                                      ncores = ncores,
                                      removeDups = removeDups,
                                      returnMeans = returnMeans)
      return(results)
    }) 
    names(results.all) <- levels(celltypes)
    
    ## Evaluate significance (cell type subsets)
  }
  
  if(!is.null(subset.list)){
    
    ## load in the subset file if it exists, or make it and probably be a good
    ## idea to save it, too
    if(!is.list(subset.list)){
      stop(paste0("`subset.list` is not a list. You can build this using: `binomialTestMatrix()` then `selectSubsets()`"))
    }
    
    combo_ids <- names(subset.list)
    d <- dist
    
    if(verbose){
      message("Calculating trends for each subset in `subset.list` with respect to the cell types in `cells$celltypes`")
    }
    
    ## initialize list
    results.all <- list()
    
    ## for each subset of cells, evaluate significance against each reference cell type across scales
    for(i in combo_ids){
      
      if(verbose){
        message(i)
      }
      
      ## get area around the subset cells to identify neighbors
      ref.buffer <- sf::st_buffer(cells[subset.list[[i]], ], d)
      ## union polygons to avoid memory overflow, too slow
      # ref.buffer_union <- sf::st_union(ref.buffer)
      # get the different types of neighbor cells that are within "d" of the ref cells
      neigh.cells <- sf::st_intersection(cells, ref.buffer$geometry)
      
      ## remove duplicate neighbor cells to prevent them from being counted multiple times
      ## and inflate the Z scores
      if(removeDups){
        ## need to remove self cells else have trivial enrichment when d~0
        self.cells <- cells[cells$celltypes == ct,]
        neigh.cells <- neigh.cells[setdiff(rownames(neigh.cells), rownames(self.cells)),]
        
        # message("number of neighbor cells before: ", nrow(neigh.cells))
        #neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), rownames(cells)),]
        ## hack to accommodate self cells that are neighbors of another self cell
        neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), c(rownames(cells), paste0(rownames(self.cells), '.1'))),]
        # message("number of neighbor cells after removing dups: ", nrow(neigh.cells))
      }
      
      ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
      ## chose to shuffle the scales in parallel, but in each scale, the perms done linearly
      ## I think a bottle neck originally was waiting for certain cell types to finish
      ## so if I split up the scales, might be able to get through each cell type faster and speed up entire process?
      ## I could also split up the permutations, but then each cell type for each scale is done one by one
      ## I could do cell types in parallel, but then for each cell type need to go through each res and each perm one by one
      results <- evaluateSignificance(cells = cells,
                                      randomcellslist = shuffle.list,
                                      trueNeighCells = neigh.cells,
                                      cellBuffer = ref.buffer,
                                      ncores = ncores,
                                      removeDups = removeDups,
                                      returnMeans = returnMeans)
      results.all[[i]] <- results
      
      rm(ref.buffer)
      rm(neigh.cells)
      rm(results)
      gc(verbose = FALSE, reset = TRUE)
    }
    
  }
  
  ## return results
  if(verbose){
    total_t <- round(difftime(Sys.time(), start_time, units="mins"), 2)
    message(sprintf("Time was %s mins", total_t))
  }
  
  return(results.all)
  
}

