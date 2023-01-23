# functions for generating simulations

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


# functions for finding trends


#' Shuffle the celltype labels at different shuffling resolutions
#'
#' @description for each resolution, and for i number of permutations, shuffle the celltype labels.
#' Returns a list of lists, where each list is a factor of the shuffled celltype labels.
#' Outer list is for each resolution of shuffling, and each inner list is for a given permutation with a given seed.
#' @param cells sp::SpatialPointsDataFrame object, with celltypes features and point geometries
#' @param resolutions numeric vector of the different resolutions to shuffle at and subsequently compute significance at
#' @param perms number of permutations to shuffle for each resolution (default = 1)
#' @param ncores number of cores for parallelization (default 1)
#' @param seed set seed for shuffling (if more than 1 permutation, then seed equals permutation number)
#' @param verbose Boolean for verbosity (default TRUE)
#'
makeShuffledCells <- function(cells, resolutions, perms = 1, ncores = 1, seed = 0, verbose = TRUE){
  
  ## effectively will make a list of lists of randomly shuffled cell labels.
  ## a list for each resolution that contains factors of shuffled cell labels for each permutation
  ## use the cell labels to reorder the labels of the `cells`
  
  randomcellslist <- lapply(resolutions, function(r) {
    
    grid <- sf::st_make_grid(cells, cellsize = r)
    
    if(verbose){
      message(r, " micron resolution")
      message(length(grid), " tiles to shuffle...")
    }
    
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
      randcelltype <- unlist(parallel::mclapply(1:length(grid), function(i) {
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
      }, mc.cores=ncores))
      
      
      
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
  
  names(randomcellslist) <- resolutions
  return(randomcellslist)
}


#' compute significant different between real and randomly shuffled cell neighbor proportions
#'
#' @param cells sp::SpatialPointsDataFrame of all the cells
#' @param randomcellslist list of lists of randomly shuffled cell type labels produced from `makeShuffledCells`
#' @param trueNeighCells Simple feature collection of real cells for a given reference cell type, with geometries of a given dist (from sf::st_buffer)
#' @param cellBuffer Simple feature collection of the neighbor cells that are within "dist" of the ref cells (from sf::intersection)
#' @param ncores number of cores for parallelization (default 1)
#'
evaluateSignificance <- function(cells, randomcellslist, trueNeighCells, cellBuffer, ncores = 1){
  
  allcells <- cells
  trueNeighCells <- trueNeighCells
  cellBuffer <- cellBuffer
  
  results <- do.call(rbind, parallel::mclapply(randomcellslist, function(cellsAtRes){
    
    ## iterate through each permutation of a given resolution
    ## produce the scores for each neighbor cell type
    ## combine the tables into a dataframe, take the mean of the scores for each cell type
    ## score means for each resolution returned as a column where rows are the cell types using sapply
    scores <- do.call(rbind, lapply(cellsAtRes, function(randomcellslabels){
      
      randomcells <- allcells
      randomcells$celltypes <- as.factor(randomcellslabels)
      sf::st_agr(randomcells) <- "constant"
      
      bufferrandomcells <- sf::st_intersection(randomcells, cellBuffer$geometry)
      
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
    })
    )
    return(colMeans(scores))
  }, mc.cores = ncores)
  )
  
  rm(allcells)
  rm(trueNeighCells)
  rm(cellBuffer)
  gc(verbose = FALSE, reset = TRUE)
  
  return(results)
}


#' find subsets of cells
#' @description find the subset cells of a reference cell type defined by being either significantly "near" or "away" with respect to a given neighbor cell type.
#' Significantly "near" or "away" defined by whether the proportion of adjacent cells within a given distance are either enriched or depleted in a given neighbor cell type compared to the global proportion of the given neighbor cell type.
#' Note that for near, or localized, just test if significantly enriched.
#' For "away", do the same, but then take the cells that were not significant. However, would recommend setting the p-value threshold to be very liberal, like 0.5, that way, only the cells that couldn't even pass a p-val cutoff of 0.5 would be selected for, and these would be expected to be very much depleted or separated from the neighbor cell type.
#'
#' @param cells sp::SpatialPointsDataFrame object, with celltypes features and point geometries
#' @param sub.dist distance to define subsets relative to a neighbor cell type (default = 100)
#' @param sub.type subset type, either ref cells "near" (ie localized) a neighbor cell type, or "away" (ie separated) from a neighbor cell type.
#' @param sub.thresh significance threshold for the binomial test (default = 0.05)
#' @param ncores number of cores for parallelization (default 1)
#' @param verbose Boolean for verbosity (default TRUE)
#' 
#' @return list where each entry is a subset and the values are the cell ids determined to be in the subset
#' 
getSubsets <- function(cells,
                       sub.dist = 100,
                       sub.type = c("near", "away"),
                       sub.thresh = 0.05,
                       ncores = 1,
                       verbose = TRUE) {
  
  start_time <- Sys.time()
  sub.type <- match.arg(sub.type)
  
  if(!sub.type %in% c("near", "away")){
    stop("`sub.type` must be either 'near' or 'away'")
  }
  
  if(verbose){
    message("Identifying subsets of cells that are ",
            sub.type,
            " a given neighbor cell type based on distance of ",
            sub.dist)
  }
  
  cell.types <- cells$celltypes
  names(cell.types) <- rownames(cells)
  
  ## setup possible subset combinations:
  combos <- expand.grid(rep(list(1:length(levels(cell.types))),2))
  colnames(combos) <- c("ref", "neighbors")
  combo_ids <- unlist(lapply(rownames(combos), function(i){
    cells.ref.ix <- levels(cell.types)[as.numeric(combos[i,1])]
    cells.neighbors.ix <- levels(cell.types)[as.numeric(combos[i,2])]
    id <- paste0(cells.ref.ix, "_", sub.type, "_", cells.neighbors.ix)
    id
  }))
  
  ## loop through possible combinations and determine cells in each subset
  subsets <- lapply(rownames(combos), function(i){
    
    cells.ref.ix <- levels(cell.types)[as.numeric(combos[i,1])]
    cells.neighbors.ix <- levels(cell.types)[as.numeric(combos[i,2])]
    id <- paste0(cells.ref.ix, "_", sub.type, "_", cells.neighbors.ix)
    message("computing subsets for ", id)
    
    cells.ref <- names(cell.types[cell.types == cells.ref.ix])
    cells.neighbors <- names(cell.types[cell.types == cells.neighbors.ix])
    
    ## get global fraction of the neighbor cell type
    globalFrac <- length(cells.neighbors)/nrow(cells)
    
    ## get area around the reference cells to test if they are enriched with the neighbor cell type
    refs.buffer <- sf::st_buffer(cells[cells$celltypes == cells.ref.ix,], sub.dist)
    
    start_time_subset <- Sys.time()
    # -------------------------------------------------------------------------
    ## for each ref cell, look for neighbor cells that are within given distance
    ## and compare to all cells within this distance to get local fraction
    tests <- unlist(parallel::mclapply(rownames(refs.buffer), function(i){
      neighs.inbuffer <- sf::st_intersection(cells[cells$celltypes == cells.neighbors.ix,], refs.buffer[i,]$geometry)
      cells.inbuffer <- sf::st_intersection(cells, refs.buffer[i,]$geometry)
      t <- stats::binom.test(x = nrow(neighs.inbuffer),
                             n = nrow(cells.inbuffer),
                             p = globalFrac,
                             alternative = "greater")
      # return the p-value for the reference cell that was tested
      t$p.value
    }, mc.cores = ncores
    ))
    # -------------------------------------------------------------------------
    
    if(sub.type == "near"){
      ## return cells that were significant
      sub.cells <- rownames(refs.buffer)[which(tests < sub.thresh)]
    }
    if (sub.type == "away"){
      ## return the cells that were not significant
      ## however, would recommend setting threshold very liberal, like 0.5,
      ## that way, only the cells that couldn't even pass a p-val cutoff of 0.5 would be selected for,
      ## and these would be expected to be very much depleted or separated from the neighbor cell type
      sub.cells <- rownames(refs.buffer)[which(tests < sub.thresh)]
      sub.cells <- rownames(refs.buffer)[which(!rownames(refs.buffer) %in% sub.cells)]
    }
    
    total_t_subset <- round(difftime(Sys.time(), start_time_subset, units = "mins"), 2)
    if(verbose){
      message(sprintf("Time to compute was %smins", total_t_subset))
    }
    
    return(sub.cells)
  })
  
  names(subsets) <- combo_ids
  
  if(verbose){
    total_t <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
    message(sprintf("Time to compute was %smins", total_t))
  }
  
  return(subsets)
  
}


# use this one, more developed:

#' Compute trends of cell type colocalization for each cell type combination across specified resolutions
#'
#' @description Trends are based on significant differences in cell type proportions between the real and randomly shuffled datasets.
#' Cell type proportions are with respect to the different cell types that are neighboring the cells of a given reference cell type within a certain defined distance.
#' This is done at difference resolutions, where a resolution is whether the cell type labels are shuffled locally or globally.
#' Trends are essentially built from significance values. The significance test basically asks if two cell types are localized or separated by assessing if the proportion of the neighboring cell type is significantly greater, or less than, random chance.
#'
#' @param pos matrix of x and y coordinates of each cell
#' @param resolutions numeric vector of the different resolutions to shuffle at and subsequently compute significance at
#' @param dist numeric distance to define neighbor cells with respect to each reference cell (default 50)
#' @param sub.dist distance to define subsets relative to a neighbor cell type (default = 100)
#' @param sub.type subset type, either "pairwise", to not use subsets, or susbets of ref cells "near" (ie localized) a neighbor cell type, or "away" (ie separated) from a neighbor cell type.
#' @param sub.thresh significance threshold for the binomial test (default = 0.05)
#' @param perms number of permutations to shuffle for each resolution (default = 1)
#' @param seed set seed for shuffling (if more than 1 permutation, then seed equals permutation number)
#' @param ncores number of cores for parallelization (default 1)
#' @param loadShuffleFile path to a preshuffled randomcellslist rds object (default NA)
#' @param saveShuffleFilePath can save the shuffled cell labels to speed things up later (default NA)
#' @param loadSubsetFile path to a premade subset list rds object (default NA)
#' @param saveSubsetFile can save the subset cell list to speed things up later, or use for subsequent plotting (default NA)
#' @param plot Boolean to return plots (default TRUE)
#' @param verbose Boolean for verbosity (default TRUE)
#'
#' @return A list that contains a dataframe for each reference cell type, where the dataframe contains the significance values for each neighbor cell type at each resolution
#' 
#' @examples 
#' 
#' export
findTrendsv2 <- function(pos,
                         celltypes,
                         resolutions,
                         dist = 50,
                         sub.dist = 100,
                         sub.type = c("pairwise", "near", "away"),
                         sub.thresh = 0.05,
                         perms = 1,
                         seed = 0,
                         ncores = 1,
                         loadShuffleFile = NA,
                         saveShuffleFilePath = NA,
                         loadSubsetFile = NA,
                         saveSubsetFile = NA,
                         plot = FALSE,
                         verbose = TRUE){
  
  sub.type <- match.arg(sub.type)
  
  if(!sub.type %in% c("pairwise", "near", "away")){
    stop("`sub.type` must be either 'pairwise', 'near' or 'away'")
  }
  
  if(length(levels(celltypes)) == 0){
    message("Warning: `celltypes` does not have levels. Creating levels from values")
    celltypes <- factor(celltypes)
    names(celltypes) <- rownames(pos)
  }
  
  seed <- seed
  
  if(verbose){
    start_time <- Sys.time()
  }
  
  ## load the real data
  # ============================================================
  if(verbose){
    message("creating `sp::SpatialPointsDataFrame`")
  }
  
  cells <- sp::SpatialPointsDataFrame(
    coords = as.data.frame(pos),
    data=data.frame(
      celltypes=celltypes
      #name=rownames(pos)
    ))
  cells <- sf::st_as_sf(cells)
  
  ## Change rowname assignments of cells to integers.
  ## Solution to keep rows in same order later on when
  ## randomly shuffling cell labels
  rownames(cells) <- as.character(1:dim(cells)[1])
  
  # make asumption that cell type attribute is constant throughout the geometries of each cell
  ## it removed the warning that keep popping up, which says this assumption is made anyways
  sf::st_agr(cells) <- "constant"
  
  print(cells)
  
  if(verbose){
    message("Generate randomly permuted background at each resolution")
  }
  
  # visualize real data
  # if(plot){
  #   plot.new()
  #   par(mfrow=c(1,1), mar=rep(1,4))
  #   plot(cells, pch=".", main = "original dataset")
  #   # dev.off()
  # }
  
  ## Generate randomly permuted background at each resolution
  # ============================================================
  
  if(assertthat::is.string(loadShuffleFile)){
    randomcellslist <- readRDS(file = loadShuffleFile)
  } else {
    
    ## parallel shuffling of the grids in each resolution
    randomcellslist <- makeShuffledCells(cells,
                                         resolutions,
                                         perms = perms,
                                         ncores = ncores,
                                         seed = seed,
                                         verbose = verbose)
    
    if(assertthat::is.string(saveShuffleFilePath)){
      saveRDS(object = randomcellslist, file = saveShuffleFilePath)
    }
  }
  
  if(verbose){
    message("Evaluating significance for each cell type")
    message("using neighbor distance of ", dist)
  }
  
  ## Evaluate significance (pairwise)
  # ============================================================
  if(sub.type == "pairwise"){
    
    if(verbose){
      message("Calculating for pairwise combinations")
    }
    
    d <- dist
    results.all <- lapply(levels(celltypes), function(ct) {
      
      if(verbose){
        message(ct)
      }
      
      # get polygon geometries of reference cells of "celltype" up to defined distance "dist"
      # use this to assess neighbors within "d" um of each cell
      ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
      # get the different types of neighbor cells that are within "d" of the ref cells
      neigh.cells <- sf::st_intersection(cells, ref.buffer$geometry)
      
      ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
      ## chose to shuffle the resolutions in parallel, but in each resolution, the perms done linearly
      ## I think a bottle neck originally was waiting for certain cell types to finish
      ## so if I split up the resolutions, might be able to get through each cell type faster and speed up entire process?
      ## I could also split up the permutations, but then each cell type for each resolution is done one by one
      ## I could do cell types in parallel, but then for each cell type need to go through each res and each perm one by one
      results <- evaluateSignificance(cells = cells,
                                      randomcellslist = randomcellslist,
                                      trueNeighCells = neigh.cells,
                                      cellBuffer = ref.buffer,
                                      ncores = ncores)
      return(results)
    }) 
    names(results.all) <- levels(celltypes)
    
    ## Evaluate significance (cell type subsets)
    # ============================================================
  } else if (sub.type %in% c("near", "away")){
    
    ## load in the subset file if it exists, or make it and probably be a good
    ## idea to save it, too
    if(assertthat::is.string(loadSubsetFile)){
      subset.list <- readRDS(file = loadSubsetFile)
    } else {
      
      ## parallel shuffling of the grids in each resolution
      subset.list <- getSubsets(cells = cells,
                                sub.dist = sub.dist,
                                sub.type = sub.type,
                                sub.thresh = sub.thresh,
                                ncores = ncores,
                                verbose = verbose)
      if(assertthat::is.string(saveSubsetFile)){
        saveRDS(object = subset.list, file = saveSubsetFile)
      }
    }
    
    combo_ids <- names(subset.list)
    d <- dist
    
    if(verbose){
      message("Calculating for subset type: ", sub.type)
      message("within subset distance of: ", sub.dist)
    }
    
    ## initialize list
    results.all <- list()
    
    ## for each subset of cells, evaluate significance against each reference cell type across resolutions
    for(i in combo_ids){
      
      if(verbose){
        message(i)
      }
      
      ## get area around the subset cells to identify neighbors
      ref.buffer <- sf::st_buffer(cells[subset.list[[i]], ], d)
      # get the different types of neighbor cells that are within "d" of the ref cells
      neigh.cells <- sf::st_intersection(cells, ref.buffer$geometry)
      
      ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
      ## chose to shuffle the resolutions in parallel, but in each resolution, the perms done linearly
      ## I think a bottle neck originally was waiting for certain cell types to finish
      ## so if I split up the resolutions, might be able to get through each cell type faster and speed up entire process?
      ## I could also split up the permutations, but then each cell type for each resolution is done one by one
      ## I could do cell types in parallel, but then for each cell type need to go through each res and each perm one by one
      results <- evaluateSignificance(cells = cells,
                                      randomcellslist = randomcellslist,
                                      trueNeighCells = neigh.cells,
                                      cellBuffer = ref.buffer,
                                      ncores = ncores)
      results.all[[i]] <- results
      
      rm(ref.buffer)
      rm(neigh.cells)
      rm(results)
      gc(verbose = FALSE, reset = TRUE)
    }
    
  } else {
    stop("`sub.type` must be either 'pairwise', 'near' or 'away'")
  }
  
  if(verbose){
    total_t <- round(difftime(Sys.time(), start_time, units="mins"), 2)
    message(sprintf("Time was %s mins", total_t))
  }
  
  return(results.all)
  
}

# ignore this one:

#' Compute trends of cell type colocalization for each cell type combination across specified resolutions
#'
#' @description Trends are based on significant differences in cell type proportions between the real and randomly shuffled datasets.
#' Cell type proportions are with respect to the different cell types that are neighboring the cells of a given reference cell type within a certain defined distance.
#' This is done at difference resolutions, where a resolution is whether the cell type labels are shuffled locally or globally.
#' Trends are essentially built from significance values. The significance test basically asks if two cell types are localized or separated by assessing if the proportion of the neighboring cell type is significantly greater, or less than, random chance.
#'
#' @param pos matrix of x and y coordinates of each cell
#' @param resolutions numeric vector of the different resolutions to shuffle at and subsequently compute significance at
#' @param dist numeric distance to define neighbor cells with respect to each reference cell (default 30)
#' @param subsetdist if subsetting, numeric distance to define subset of reference cells within distance to another cell type (default NaN)
#' @param ncores number of cores for parellelization (default 1)
#' @param plots Boolean to return plots (default TRUE)
#' @param verbose Boolean for verbosity (default TRUE)
#' @param seed set the seed for shuffling (default = 0)
#'
#' @return A list that contains a dataframe for each reference cell type, where the dataframe contains the significance values for each neighbor cell type at each resolution
#' 
#' @examples 
#' 
#' export
findTrends <- function(pos, celltypes, resolutions, dist = 30, subsetdist = NaN, ncores = 1, plot = FALSE, verbose = TRUE, seed = 0){
  
  if(verbose){
    start_time <- Sys.time()
  }
  
  ## load the real data
  # ============================================================
  if(verbose){
    message("creating `sp::SpatialPointsDataFrame`")
  }
  
  cells <- sp::SpatialPointsDataFrame(
    coords = as.data.frame(pos),
    data=data.frame(
      celltypes=celltypes
      #name=rownames(pos)
    ))
  cells <- sf::st_as_sf(cells)
  
  # make asumption that cell type attribute is constant throughout the geometries of each cell
  ## it removed the warning that keep popping up, which says this assumption is made anyways
  sf::st_agr(cells) <- "constant"
  
  if(verbose){
    message("Generate randomly permuted background at each resolution")
  }
  
  # visualize real data
  if(plot){
    plot.new()
    #   par(mfrow=c(1,1), mar=rep(1,4))
    #   plot(cells, pch=".", main = "original dataset")
    #   # dev.off()
  }
  
  ## Generate randomly permuted background at each resolution
  # ============================================================
  randomcellslist <- lapply(resolutions, function(r) {
    
    grid = sf::st_make_grid(cells, cellsize = r)
    
    if(verbose){
      message(r, " micron resolution")
      message(length(grid), " tiles to shuffle...")
    }
    
    ## shuffle within grid once
    randcelltype <- unlist(parallel::mclapply(1:length(grid), function(i) {
      ## sometimes can be on boundary, so just pick first one
      int <- sf::st_intersection(cells, grid[[i]])
      set.seed(seed)
      # randomly grab cell labels for cell in the grid
      shuffled_cells <- sample(int$celltypes)
      # assign the cell ids to the randomly sampled cell labels
      names(shuffled_cells) <- rownames(int)
      shuffled_cells
    }, mc.cores = ncores))
    # reorder so cells are "1", "2", etc. and not in the order based on grids
    randcelltype <- randcelltype[as.character(sort(as.numeric(names(randcelltype))))]
    
    # generate the "cells" SpatialPointsDataframe for the shuffled labels
    randomcells <- cells
    randomcells$celltypes <- randcelltype
    
    # visualize shuffled cells
    if(plot){
      par(mfrow=c(1,1), mar=rep(2,4))
      # plot(grid)
      plot(randomcells, pch=".", add=TRUE, main = paste0(r, " resolution shuffled"))
      # dev.off()
    }
    
    sf::st_agr(randomcells) <- "constant"
    
    return(randomcells)
  })
  names(randomcellslist) <- resolutions
  
  if(verbose){
    message("Evaluating significance for each cell type")
    message("using neighbor distance of ", dist)
  }
  
  ## Evaluate significance (pairwise)
  # ============================================================
  if(is.na(subsetdist)){
    
    if(verbose){
      message("Calculating for pairwise combinations")
    }
    
    d <- dist
    results.all <- parallel::mclapply(levels(celltypes), function(ct) {
      
      if(verbose){
        message(ct)
      }
      
      # get polygon geometries of reference cells of "celltype" up to defined distance "dist"
      buffer1 <- sf::st_buffer(cells[cells$celltypes == ct,], d) # assessing neighbors within 30 um of each cell for example; think of this like the K for knn
      # the neighbor cells that are within "dist" of the ref cells
      buffer1cells <- sf::st_intersection(cells, buffer1$geometry)
      
      # compare the real neighbor cell results above to the random shuffles at each resolution
      results <- sapply(randomcellslist, function(randomcells){
        buffer1randomcells <- sf::st_intersection(randomcells, buffer1$geometry)
        
        ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
        y1 <- table(buffer1cells$celltypes)
        y2 <- table(buffer1randomcells$celltypes)
        n1 <- length(buffer1cells$celltypes)
        n2 <- length(buffer1randomcells$celltypes)
        p1 <- y1/n1
        p2 <- y2/n2
        p <- (y1+y2)/(n1+n2)
        Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
        return(Z)
      })
      return(t(results))
    }, mc.cores=ncores) # im getting an error is this is higher, at least for the 3 ct sims with circles...
    names(results.all) <- levels(celltypes)
    
    ## Evaluate significance (triplet subsetting)
    # ============================================================
  } else if (isTRUE(as.double(subsetdist) > 0)){
    
    if(verbose){
      message("Calculating for reference subsets defined by subset distance of ", subsetdist)
    }
    
    ## the subset combinations
    combos <- expand.grid(rep(list(1:length(levels(celltypes))),2))
    colnames(combos) <- c("ref", "neighbors")
    
    combo_ids <- unlist(lapply(rownames(combos), function(i){
      cells.ref.ix <- levels(celltypes)[as.numeric(combos[i,1])]
      cells.neighbors.ix <- levels(celltypes)[as.numeric(combos[i,2])]
      id <- paste0(cells.ref.ix, "_near_", cells.neighbors.ix)}))
    
    d <- dist
    sd <- subsetdist
    results.all <- parallel::mclapply(1:nrow(combos), function(i) {
      
      ## get reference subset that is within `subsetdist` of another given cell type
      ct1 <- levels(celltypes)[combos[i,1]]
      ct2 <- levels(celltypes)[combos[i,2]]
      
      if(verbose){
        message(ct1, " cells within ", subsetdist, " microns of ", ct2, " cells")
      }
      
      ## get area around the cells for which the subset is defined by
      bufferct2 <- sf::st_buffer(cells[cells$celltypes == ct2,], sd)
      ## look for ref cells within this area of ct2
      subsetnearbufferct2 <- sf::st_intersection(cells[cells$celltypes == ct1,], bufferct2$geometry)
      
      ## now for this subset, test for association with all cell-types, using the defined `dist`
      buffersubsetcells <- sf::st_buffer(subsetnearbufferct2, d)
      cellsnearsubset <- sf::st_intersection(cells, buffersubsetcells$geometry)
      
      # compare the real neighbor cell results above to the random shuffles at each resolution
      results <- sapply(randomcellslist, function(randomcells){
        buffer1randomcells <- sf::st_intersection(randomcells, buffersubsetcells$geometry)
        
        ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
        y1 <- table(cellsnearsubset$celltypes)
        y2 <- table(buffer1randomcells$celltypes)
        n1 <- length(cellsnearsubset$celltypes)
        n2 <- length(buffer1randomcells$celltypes)
        p1 <- y1/n1
        p2 <- y2/n2
        p <- (y1+y2)/(n1+n2)
        Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
        return(Z)
      })
      return(t(results))
    }, mc.cores=ncores) # im getting an error is this is higher, at least for the 3 ct sims with circles...
    names(results.all) <- combo_ids
    
  } else {
    stop("`subsetdist` needs to be either NaN for pairwise comparisons or a positive distance.")
  }
  
  if(verbose){
    total_t <- round(difftime(Sys.time(), start_time, units="mins"), 2)
    message(sprintf("Time was %s mins", total_t))
  }
  
  # if(plot){
  #   pdf(trendPlotName, width=8, height=8)
  #   par(mfrow=c(length(levels(celltypes)),length(levels(celltypes))), mar=rep(2,4))
  #   sapply(levels(celltypes), function(ct1) {
  #     # print(ct1)
  #     results.norm <- results.all[[ct1]]
  #     results.norm[is.nan(results.norm)] <- NA
  #     results.norm[is.infinite(results.norm)] <- NA
  #     sapply(colnames(results.norm), function(ct2) {
  #       rg <- max(abs(results.norm[, ct2]))
  #       plot(resolutions, results.norm[,ct2], type="l", main=paste0(ct1,'\n', ct2), ylim=c(-rg, rg))
  #       abline(h = -2, col='red')
  #       abline(h = 2, col='red')
  #     })
  #   })
  #   dev.off()
  # }
  
  return(results.all)
  
}


# functions for processing output


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


# functions for visualization

#' @param idcol if results are a data.frame, this is the column that contains the additional feature to plot multiple trend lines with
#' by default it is 'id', which would be the name given to the column describing each list after melting together into a single dataframe
#' using `meltResultsList()`
#' @param ... additional plotting parameters for base R plotting. Fed into "lines()" in script
plotTrends <- function(results, idcol = "id", figPath = "results.pdf", width = 8, height = 8, legend = TRUE, ...){
  
  
  ## setup to check if original list output from `findTrends`, and plot one way
  ## or if a melted dataframe with ids from merging multiple trend analyses, check if dataframe and plot that way
  
  ## if in original list format from `findTrends`:
  if(inherits(results, "list")){
    message("results detected to be a list")
    
    pdf(figPath, width=width, height=height)
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
        resolutions <- rownames(results.norm)
        
        ## instantiate a plot and plot trend
        plot(resolutions, results.norm[,ct2],
             type="l", lwd = 2,
             main=paste0(ct1,' ref \n', ct2, " neighbors"),
             ylim=c(-rg, rg),
             xlab="resolution", ylab="Z", ...)
        
        ## threshold lines
        abline(h = -2, col='red')
        abline(h = 2, col='red')
      })
    })
    dev.off()
    
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
    
    pdf(figPath, width=width, height=height)
    # par(mfrow=c(length(refs), length(neighs)),
    #     mar=rep(4,4))
    par(mfrow=c(length(neighs), length(refs)),
        mar=rep(4,4))
    # par(mfrow=c(2, 6),
    #     mar=rep(4,4))
    
    ## for each reference cell type...(rows)
    sapply(refs, function(ct1) {
      # print(ct1)
      results.norm <- results[results[,"reference"] == ct1,]
      results.norm[is.nan(results.norm[,"Z"]), "Z"] <- NA
      results.norm[is.infinite(results.norm[,"Z"]), "Z"] <- NA
      
      ## for each neighbor cell type...(columns)
      sapply(neighs, function(ct2) {
        results.norm.neigh <- results.norm[results.norm[,"neighbor"] == ct2,]
        
        yl <- max(abs(results.norm.neigh[, "Z"]), na.rm = TRUE)
        xl <- max(as.numeric(results.norm.neigh[,"resolution"]))
        
        if(is.infinite(yl)){
          yl <- 2.0
        }
        
        ## instantiate a plot
        plot(0, 0, type = "n",
             main=paste0(ct1,' ref \n', ct2, " neighbors"),
             cex.main=1,
             ylim=c(-yl, yl),
             xlim=c(0, xl),
             xlab="resolution", ylab="Z")
        
        ## for each id param, draw a line on plot instance
        for(i in 1:length(ids)){
          id <- ids[i]
          
          results.norm.neigh.id <- results.norm.neigh[results.norm.neigh[,idcol] == id,]
          
          lines(as.numeric(results.norm.neigh.id[,"resolution"]), results.norm.neigh.id[,"Z"],
                type="l", lwd=2, col=cl[i], ...)
        }
        
        ## threshold lines
        abline(h = -2, col='red')
        abline(h = 2, col='red')
        if(legend){
          legend("topright", inset=c(-0.4,0), xpd=TRUE, legend = ids, col=cl, pch=20, cex=0.5, title = "ids")
        }
      })
    })
    dev.off()
    
  } else {
    stop("`results` are neither a list from `findTrends` or a melted data.frame from `meltResultsList`")
  }
  
}


#' This one overlays each neighbor trend wrt the same reference cell type on the plot
#' @param ... additional plotting parameters for base R plotting. Fed into "lines()" in script
plotTrendsOverlay <- function(results, figPath = "results.pdf", width = 4, height = 10, legend = TRUE, ...){
  
  
  ## setup to check if original list output from `findTrends`, and plot one way
  ## or if a melted dataframe with ids from merging multiple trend analyses, check if dataframe and plot that way
  
  ## if in original list format from `findTrends`:
  if(inherits(results, "list")){
    message("results detected to be a list")
    
    pdf(figPath, width=width, height=height)
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
        resolutions <- rownames(results.norm)
        
        ## instantiate a plot and plot trend
        plot(resolutions, results.norm[,ct2],
             type="l", lwd = 2,
             main=paste0(ct1,' ref \n', ct2, " neighbors"),
             ylim=c(-rg, rg),
             xlab="resolution", ylab="Z", ...)
        
        ## threshold lines
        abline(h = -2, col='red')
        abline(h = 2, col='red')
      })
    })
    dev.off()
    
    ## if a melted dataframe,
    ## will have an additional column that can serve to plot
    ## several trend lines on the same plot instance
    ## for example, ref vs neigh at different distances
  } else if(inherits(results, "data.frame")){
    message("results detected to be a data.frame")
    
    results <- results[,c("resolution", "neighbor", "reference", "Z")]
    
    refs <- unique(results[,"reference"])
    neighs <- unique(results[,"neighbor"])
    cl <- rainbow(length(neighs))
    
    pdf(figPath, width=width, height=height)
    # par(mfrow=c(length(refs), length(neighs)),
    #     mar=rep(4,4))
    par(mfrow=c(length(refs),1),
        mar=c(4,4,4,8)) ## bot, top, left, right
    # par(mfrow=c(2, 6),
    #     mar=rep(4,4))
    
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
      xl <- max(as.numeric(results.norm.limits[,"resolution"]))
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
           xlab="resolution", ylab="Z")
      
      ## for each neighbor cell type draw a line on plot instance
      for(i in 1:length(neighs)){
        ct2 <- neighs[i]
        # ignore showing trends with self because typically very large and
        # masks relationships with other cell types
        if(ct1 != ct2){
          results.norm.neigh.id <-  results.norm[results.norm[,"neighbor"] == ct2,]
          lines(as.numeric(results.norm.neigh.id[,"resolution"]), results.norm.neigh.id[,"Z"],
                type="l", lwd=0.8, col=cl[i], ...)
        }
      }
      
      ## threshold lines
      abline(h = -1, col='black')
      abline(h = 1, col='black')
      if(legend){
        legend("topright", inset=c(-0.4,0), xpd=TRUE, legend = neighs, col=cl, pch=20, cex=0.5, title = "neighbors")
      }
        
    })
    dev.off()
  } else {
    stop("`results` are neither a list from `findTrends` or a melted data.frame from `meltResultsList`")
  }
  
}


#ignore this, the plotTrends function above has been developed more

# plotTrends2 <- function(results, figPath = "results.pdf", width = 8, height = 8, ...){
#   
#   refs <- unique(results[,"reference"])
#   neighs <- unique(results[,"neighbor"])
#   ids <- levels(results[,"id"])
#   
#   cl <- rainbow(length(ids))
#   
#   pdf(figPath, width=width, height=height)
#   par(mfrow=c(length(refs), length(neighs)),
#       mar=rep(4,4))
#   
#   ## for each reference cell type...
#   sapply(refs, function(ct1) {
#     # print(ct1)
#     results.norm <- results[results[,"reference"] == ct1,]
#     results.norm[is.nan(results.norm[,"Z"]), "Z"] <- NA
#     results.norm[is.infinite(results.norm[,"Z"]), "Z"] <- NA
#     
#     ## for each neighbor cell type...
#     sapply(neighs, function(ct2) {
#       results.norm.neigh <- results.norm[results.norm[,"neighbor"] == ct2,]
#       
#       yl <- max(abs(results.norm.neigh[, "Z"]), na.rm = TRUE)
#       xl <- max(as.numeric(results.norm.neigh[,"resolution"]))
#       
#       ## instantiate a plot
#       plot(0, 0, type = "n",
#            main=paste0(ct1,' ref \n', ct2, " neighbors"),
#            ylim=c(-yl, yl),
#            xlim=c(0, xl),
#            xlab="resolution", ylab="Z")
#       
#       ## for each id param, draw a line on plot instance
#       for(i in 1:length(ids)){
#         id <- ids[i]
#         results.norm.neigh.id <- results.norm.neigh[results.norm.neigh[,"id"] == id,]
#         
#         lines(as.numeric(results.norm.neigh.id[,"resolution"]), results.norm.neigh.id[,"Z"],
#               type="l", lwd=2, col=cl[i], ...)
#       }
#       
#       ## threshold lines
#       abline(h = -2, col='red')
#       abline(h = 2, col='red')
#       
#       legend("topleft", legend = ids, col=cl, pch=20, cex=0.5, title = "ids")
#       
#     })
#   })
#   
#   dev.off()
#   
# }


#originally set up to visualize trends and color by different trend "clusters", ie similar trends.
#Still useful but used as much anymore


## dat = the data.frame of pvals, trend cluster assignments, etc for each pairwise combo. i.e. PKHL
## could subset this to get specific interactions
## clusters = the name of the column that has the clusters to color by
## nc = the number of unique clusters to color by.
## default should be number of unique entries in the clusters column
## but its possible that one dataset might not have a cluster because all 
## datasets clustered together so a cluster could be specific to a particular
## datasets
## colors = by default used rainbow(nc), but can specifically change to vector 
## to be used in scale_color_manual

vizTrends <- function(dat, clusters, yaxis = "zscore",
                      sig = -log10(0.05/nrow(dat)), ## sig thresh for num tests
                      nc = length(unique(dat[[clusters]])),
                      colors = rainbow(nc),
                      title = NA){
  
  plt <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_point(ggplot2::aes(x = resolution, y = .data[[yaxis]], color = .data[[clusters]]), size = 1.5) +
    ggplot2::geom_path(ggplot2::aes(x = resolution, y = .data[[yaxis]], color = .data[[clusters]]), size = 1.5) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_hline(yintercept = 0, color = "black", size = 1) +
    ggplot2::geom_hline(yintercept = sig, color = "red", size = 0.6) +
    ggplot2::geom_hline(yintercept = -sig, color = "red", size = 0.6) +
    ggplot2::facet_grid(neighbor ~ reference) +
    ggplot2::scale_x_log10() +
    # ggplot2::scale_y_continuous(trans = ggallin::pseudolog10_trans) +
    # ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$y)-axisAdj, max(dat$y)+axisAdj)) +
    # ggplot2::scale_y_log10() +
    ggplot2::ggtitle(title) +
    
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, color = "black", angle = -90, vjust = 0.5, hjust = 0),
                   axis.text.y = ggplot2::element_text(size=12, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_line(size = 0.1, colour = "black"),
                   panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                   axis.line = ggplot2::element_line(size = 0, colour = "black"),
                   panel.spacing = ggplot2::unit(0.1, "lines"),
                   strip.text = ggplot2::element_text(size = 12),
                   legend.title = ggplot2::element_blank()
                   # legend.position="none"
    )
  plt
  
}


# functions for trend comparison


#where each row is a cell type combination with a tag for a given sample
#and each column is the zscore of this for a given resolution


#' Test if intra sample vs inter samples trends are significantly different
#' 
#' @description Computes summed distance between z scores at each resolution for two given trends.
#'     Gets summed distances for combinations of intra patient trends, and inter patient trends.
#'     T-test between intra and inter distances to test for significant difference.
#'     Also can return heatmap matrix of distances between each dataset and distributions.
#'     The idea is that if inter trends are similar, then their distances will be small like the intra trends.
#'     But if the inter trends are different, their distances should be larger than the intra trends.
#'     Likewise, if the intra trends are highly variable, then might be expected to also have large differences between inter as well.
#' 
#' @param samples list of sample set names for each patient ex: list(c("PKHL", "XXCD"),c("KSFB", "NGPL"),c("PBVN", "FSLD"))
#' @param refID name of reference cell type for trend of interest
#' @param neighID name of neighbor cell type for trend of interest
#' @param zscores table of zscores for each cell type combo (rows) for each sample vs each resolution tested (columns)
#' @param heatmap return heatmap of distances between sample trends (boolean, default: TRUE)
#' @param distplot return plot of distirbutions of intra and inter trend distances (boolean, default: TRUE)
#'
#' @return pvalues for each patient
#' 
#' @examples 
#' diffTrendTesting(list(c("PKHL", "XXCD"),c("KSFB", "NGPL"),c("PBVN", "FSLD")),
#'                  refID = "Ki67 proliferating",
#'                  neighID = "CD4 Memory T cells",
#'                  zscores = combined_zscores)
#' 
#'
diffTrendTesting <- function(samples, refID, neighID, zscores, heatmap = TRUE, distplot = TRUE){
  
  r <- refID
  n <- neighID
  id <- paste0(r, " vs ", n)
  m <- expand.grid(rep(list(unlist(samples)),2))
  
  ## get table of summed distances for each trend combination
  distances <- unlist(lapply(1:nrow(m), function(i){
    t1 <- paste0(id, "_", as.vector(m[i,1]))
    t2 <- paste0(id, "_", as.vector(m[i,2]))
    t1_zscores <- zscores[t1,]
    t2_zscores <- zscores[t2,]
    sum(abs(t1_zscores - t2_zscores))
  }))
  m$dist <- distances
  
  ## get list of inter trend comparisons for each set of patient datasets
  s <- 1:length(samples)
  inter.pairs <- lapply(s, function(ix){
    patient <- samples[[which(s %in% ix)]]
    others <- unlist(samples[which(!s %in% ix)])
    var1 <- c()
    for(sample in patient){
      var1 <- c(var1, rep(sample, length(others)))
    }
    var2 <- rep(others, length(patient))
    inter_pairs <- cbind(var1, var2) 
  })
  
  ## vector of distances for the intra trends for each patient
  intra.vals <- unlist(lapply(samples, function(i){
    m[m$Var1 == i[1] & m$Var2 == i[2], "dist"]
  }))
  
  ## list of distances between inter trend comparisons for each set of patient datasets
  inter.vals <- lapply(1:length(inter.pairs), function(ix){
    pairs <- inter.pairs[[ix]]
    vals <- unlist(lapply(1:nrow(pairs), function(i){
      m[m$Var1 == pairs[i,1] & m$Var2 == pairs[i,2], "dist"]
    }))
    vals
  })
  
  ## get distributions of the distances
  intra.dists <- MASS::fitdistr(intra.vals, "normal")
  inter.dists <- lapply(1:length(inter.vals), function(ix){
    MASS::fitdistr(inter.vals[[ix]], "normal")
  })
  
  ## perform T-tests
  inter.tests <- lapply(1:length(inter.vals), function(ix){
    t.test(intra.vals, inter.vals[[ix]])
  })
  
  ## plotting
  if(heatmap){
    dat <- m
    plt <- ggplot2::ggplot(data = dat) + ggplot2::geom_tile(ggplot2::aes(x = Var1, y = Var2, fill = dist)) + 
      ggplot2::scale_y_discrete(breaks = as.character(dat$Var2), 
                                labels = as.character(dat$Var2)) +
      ggplot2::labs(title = id,
                    x = "samples",
                    y = "samples")
    
    plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, 
                                                                    color = "black",
                                                                    hjust = 0.5,
                                                                    vjust = 0.5),
                                axis.text.y = ggplot2::element_text(size = 12,
                                                                    color = "black"),
                                axis.title.y = ggplot2::element_text(size = 13),
                                axis.title.x = ggplot2::element_text(size = 13),
                                plot.title = ggplot2::element_text(size = 15),
                                legend.text = ggplot2::element_text(size = 15,
                                                                    colour = "black"),
                                legend.title = ggplot2::element_text(size = 15,
                                                                     colour = "black",
                                                                     angle = 90),
                                panel.background = ggplot2::element_blank(),
                                panel.border = ggplot2::element_rect(fill = NA,
                                                                     color = "black",
                                                                     size = 2),
                                plot.background = ggplot2::element_blank()) +
      ggplot2::geom_text(ggplot2::aes(x = as.character(Var1),
                                      y = as.character(Var2),
                                      label = format(round(dist,1),nsmall = 2))) +
      ggplot2::scale_fill_gradientn(limits = c(0,max(m$dist)),
                                    breaks = c(0,max(m$dist)),
                                    colors = (grDevices::colorRampPalette(c("white","red")))(n = 209)) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Distance",
                                                     title.position = "left",
                                                     title.hjust = 0.5,
                                                     ticks.colour = "black",
                                                     ticks.linewidth = 2,
                                                     frame.colour = "black",
                                                     frame.linewidth = 2,
                                                     label.hjust = 0)) +
      ggplot2::coord_fixed()
    print(plt)
  }
  
  ## return plot of distributions of intra and inter patient trend distances
  if(distplot){
    xlim <- max(unlist(inter.vals)) * 1.1
    xvals <- seq(0,xlim,length = 500)
    
    intra <- dnorm(xvals, mean = intra.dists$estimate[[1]], sd = intra.dists$estimate[[2]])
    
    df <- do.call(cbind.data.frame, lapply(1:length(inter.dists), function(ix){
      dist <- inter.dists[[ix]]
      dnorm(xvals, mean = dist$estimate[[1]], sd = dist$estimate[[2]])
    }))
    colnames(df) <- paste0("inter_", 1:length(inter.dists))
    df[["intra"]] <- intra
    
    dat <- reshape2::melt(df)
    dat[["xvals"]] <- rep(xvals, length(inter.dists)+1)
    
    cols <- rainbow(n = length(inter.dists))
    plt <- ggplot2::ggplot(data = dat) +
      ggplot2::geom_line(ggplot2::aes(x = xvals, y = value, color = variable)) +
      ggplot2::scale_color_manual(values = c(cols, "black")) +
      ggplot2::labs(title = id,
                    x = "distance",
                    y = "density") +
      
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
                     panel.grid.major.y =  ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(size = 1, colour = "black")
                     # legend.position="none"
      )
    print(plt)
    
  }
  
  return(inter.tests)
  
}
