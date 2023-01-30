#' Shuffle the celltype labels at different resolutions
#'
#' @description For each resolution, and for i number of permutations, shuffle the celltype labels.
#' Returns a list of lists, where each list is a factor of the shuffled celltype labels.
#' Outer list is for each resolution of shuffling, and each inner list is for a given permutation with a given seed.
#'
#' @param cells sp::SpatialPointsDataFrame object; celltypes features and point geometries
#' @param resolutions numeric vector; the different resolutions to shuffle at and subsequently compute significance at
#' @param perms numeric; number of permutations to shuffle for each resolution (default 1)
#' @param ncores numeric; number of cores for parallelization (default 1)
#' @param seed numeric; set seed for shuffling (if more than 1 permutation, then seed equals permutation number)
#' @param verbose Boolean; verbosity (default TRUE)
#'
#' @return list;
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



#' Find subsets of cells
#'
#' @description find the subset cells of a reference cell type defined by being either significantly "near" or "away" with respect to a given neighbor cell type.
#' Significantly "near" or "away" defined by whether the proportion of adjacent cells within a given distance are either enriched or depleted in a given neighbor cell type compared to the global proportion of the given neighbor cell type.
#' Note that for near, or localized, just test if significantly enriched.
#' For "away", do the same, but then take the cells that were not significant. However, would recommend setting the p-value threshold to be very liberal, like 0.5, that way, only the cells that couldn't even pass a p-val cutoff of 0.5 would be selected for, and these would be expected to be very much depleted or separated from the neighbor cell type.
#'
#' @param cells sp::SpatialPointsDataFrame object;  with celltypes features and point geometries
#' @param sub.dist distance to define subsets relative to a neighbor cell type (default = 50)
#' @param sub.type subset type, either ref cells "near" (ie localized) a neighbor cell type, or "away" (ie separated) from a neighbor cell type.
#' @param sub.thresh significance threshold for the binomial test (default = 0.05)
#' @param ncores number of cores for parallelization (default 1)
#' @param verbose Boolean for verbosity (default TRUE)
#'
#' @return list; each entry is a subset and the values are the cell ids determined to be in the subset
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




#' Evaluate significance
#'
#' @description Compute significant different between real and randomly shuffled cell neighbor proportions.
#'
#' @param cells sp::SpatialPointsDataFrame; all the cells
#' @param randomcellslist nested list; of randomly shuffled cell type labels produced from `makeShuffledCells`
#' @param trueNeighCells Simple feature collection of real cells for a given reference cell type, with geometries of a given dist (from sf::st_buffer)
#' @param cellBuffer Simple feature collection of the neighbor cells that are within "dist" of the ref cells (from sf::intersection)
#' @param ncores numeric; number of cores for parallelization (default 1)
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




#' Compute trends of cell type colocalization for each cell type combination across specified resolutions
#'
#' @description Trends are based on significant differences in cell type proportions between the real and randomly shuffled datasets.
#' Cell type proportions are with respect to the different cell types that are neighboring the cells of a given reference cell type within a certain defined distance.
#' This is done at difference resolutions, where a resolution is whether the cell type labels are shuffled locally or globally.
#' Trends are essentially built from significance values. The significance test basically asks if two cell types are localized or separated by assessing if the proportion of the neighboring cell type is significantly greater, or less than, random chance.
#'
#' @param pos data frame; x and y coordinates of each cell
#' @param celltypes character vector; the cell type of each cell provided in pos
#' @param resolutions numeric vector; the different resolutions to shuffle at and subsequently compute significance at
#' @param dist numeric; distance to define neighbor cells with respect to each reference cell (default = 50)
#' @param sub.dist numeric; distance to define subsets relative to a neighbor cell type (default = 100)
#' @param sub.type character; subset type, either "pairwise", to not use subsets, or subsets of ref cells "near" (ie localized) a neighbor cell type, or "away" (ie separated) from a neighbor cell type
#' @param sub.thresh numeric; significance threshold for the binomial test (default = 0.05)
#' @param perms numeric; number of permutations to shuffle for each resolution (default = 1)
#' @param seed numeric; set seed for shuffling (if more than 1 permutation, then seed equals permutation number)
#' @param ncores numeric; number of cores for parallelization (default 1)
#' @param loadShuffleFile character; path to a preshuffled randomcellslist rds object (default NA)
#' @param saveShuffleFilePath character; can save the shuffled cell labels to speed things up later (default NA)
#' @param loadSubsetFile character; path to a premade subset list rds object (default NA)
#' @param saveSubsetFile character; can save the subset cell list to speed things up later, or use for subsequent plotting (default NA)
#' @param plot Boolean; return plots (default TRUE)
#' @param verbose Boolean; verbosity (default TRUE)
#'
#' @return A list that contains a dataframe for each reference cell type, where the dataframe contains the significance values for each neighbor cell type at each resolution
#'
#' @export
findTrends <- function(pos,
                         celltypes,
                         resolutions = c(50, 100, 200, 300),
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

    if(length(levels(celltypes)) == 0){
        message("Warning: 'celltypes' does not have levels. Creating levels from values")
        celltypes <- factor(celltypes)
        names(celltypes) <- rownames(pos)
    }

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
        data = data.frame(
            celltypes = celltypes
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
