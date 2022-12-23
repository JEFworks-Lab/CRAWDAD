## --------------------------------------
## libraries and functions

library(parallel)
library(sf)
library(sp)
library(ggplot2)
library(assertthat)
library(reshape2)
library(MASS)
library(stats)

source("/home/bmille79/code/multiscale_trend_analysis/multiscale_trend_analysis_functions.r")

## commands loaded in from command line, but really from the slurm script. A different job per dataset, and in its own node with x number of cores for parallelization
args <- commandArgs(TRUE)

## --------------------------------------
## load arguments from the command line

## path to the metadata table of cell coordinates and celltype names
meta <- read.csv2(file = args[1], row.names = 1)

## name of column with cell type annotations to use
## for example, "celltypes", or "celltypes_folBcombined"
meta <- meta[,c("x", "y", as.character(args[2]))]

subsetType <- as.character(args[3])

subdist <- as.numeric(args[4]) ## distance to define subsets

if(subsetType == "near"){
  subThresh <- 0.05
} else if(subsetType == "away"){
  subThresh <- 0.5
} else (
  stop("subsetType must be either 'near' or 'away'")
)

## full path to the subset dataset list rds
subsetDataPath <- args[5]

ncs <- args[6] ## number of cores


message(args[1])
message(args[2])
message(args[3])
message(args[4])
message(args[5])
message(args[6])

## --------------------------------------
## run the function

pos <- meta[,c("x", "y")]
celltypes <- meta[,3]

if(length(levels(celltypes)) == 0){
  message("Warning: `celltypes` does not have levels. Creating levels from values")
  celltypes <- factor(celltypes)
  names(celltypes) <- rownames(pos)
}

cells <- sp::SpatialPointsDataFrame(
    coords = as.data.frame(pos),
    data = data.frame(
        celltypes = celltypes
        #name=rownames(pos)
))
cells <- sf::st_as_sf(cells)
# make asumption that cell type attribute is constant throughout the geometries of each cell
## it removed the warning that keep popping up, which says this assumption is made anyways
sf::st_agr(cells) <- "constant"



## if subset data exists, set up parameters to load it in, otherwise set them up to save it once made
if(file.exists(subsetDataPath)){
  loadSubset <- subsetDataPath
  saveSubset <- NA
} else {
  loadSubset <- NA
  saveSubset <- subsetDataPath
}

if(assertthat::is.string(loadSubset)){
    subset.list <- readRDS(file = loadSubset)
} else {
    
    ## parallel shuffling of the grids in each resolution
    subset.list <- getSubsets(cells = cells,
                            sub.dist = subdist,
                            sub.type = subsetType,
                            sub.thresh = subThresh,
                            ncores = ncs,
                            verbose = TRUE)
    if(assertthat::is.string(saveSubset)){
        saveRDS(object = subset.list, file = saveSubset)
    }
}