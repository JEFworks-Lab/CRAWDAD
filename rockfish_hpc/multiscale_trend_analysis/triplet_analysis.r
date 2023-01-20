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
library(dplyr)

source("/home/bmille79/code/multiscale_trend_analysis/multiscale_trend_analysis_functions.r")

## commands loaded in from command line, but really from the slurm script. A different job per dataset, and in its own node with x number of cores for parallelization
args <- commandArgs(TRUE)

## --------------------------------------
## load arguments from the command line

## path to the metadata table of cell coordinates and celltype names
meta <- read.csv2(file = args[1], row.names = 1, sep=',')

## name of column with cell type annotations to use
## for example, "celltypes", or "celltypes_folBcombined"
meta <- meta[,c("x", "y", as.character(args[2]))]

## make sure the coordinates are numeric
meta <- meta %>% 
  dplyr::mutate_at(vars(x, y), as.numeric)

## full path to the shuffled data rds
shuffledDataPath <- args[3]

subsetType <- args[4]

if(subsetType == "near"){
  subThresh <- 0.05
} else if(subsetType == "away"){
  subThresh <- 0.5
} else (
  stop("subsetType must be either 'near' or 'away'")
)

subdist <- as.numeric(args[5]) ## distance to define subsets

## full path to the subset dataset list rds
subsetDataPath <- args[6]

## full path and name of the output rds file to save (which is the list format output of the function `findTrendsv2`)
outfile <- args[7]

ncs <- args[8] ## number of cores


message(args[1])
message(args[2])
message(args[3])
message(args[4])
message(args[5])
message(args[6])
message(args[7])
message(args[8])

## --------------------------------------
## run the function

# distances <- c(30, 50, 100, 200)
distances <- c(200)
resolutions <- c(50, 100, 200, 300, 500, 600, 700, 800, 900, 1000, 1200, 1500, 3000)

# resolutions <- c(50, 100, 200, 300, 500, 600, 700, 800)


# resolutions <- c(3000)
prms <- 1 ## number of permutations
sd <- 1 ## seed



## if subset data exists, set up parameters to load it in, otherwise set them up to save it once made
if(file.exists(subsetDataPath)){
  loadSubset <- subsetDataPath
  saveSubset <- NA
} else {
  loadSubset <- NA
  saveSubset <- subsetDataPath
}

triplet_results <- lapply(distances, function(d){
  
  # shuffledData <- paste0("/home/brendan/projects/HuBMAP/data/datasets/spleen/pkhl_shuffled_res_50-3000.rds")
  if(file.exists(shuffledDataPath)){
    
    results <- findTrendsv2(pos = meta[,c("x", "y")],
                             celltypes = meta[,3],
                             resolutions,
                             dist = d,
                             sub.dist = subdist,
                             sub.type = subsetType,
                             sub.thresh = subThresh,
                             perms = prms,
                             ncores = ncs,
                             verbose = TRUE,
                             loadShuffleFile = shuffledDataPath,
                             loadSubsetFile = loadSubset,
                             saveSubsetFile = saveSubset,
                             seed = sd)
    # gc()
  } else {
    results <- findTrendsv2(pos = meta[,c("x", "y")],
                             celltypes = meta[,3],
                             resolutions,
                             dist = d,
                             sub.dist = subdist,
                             sub.type = subsetType,
                             sub.thresh = subThresh,
                             perms = prms,
                             ncores = ncs,
                             verbose = TRUE,
                             saveShuffleFilePath = shuffledDataPath,
                             loadSubsetFile = loadSubset,
                             saveSubsetFile = saveSubset,
                             seed = sd)
    # gc()
  }
  # print(results)
  return(results)
  
})

names(triplet_results) <- distances

saveRDS(object = triplet_results, file = outfile)

# saveRDS(object = pkhl_triplet_results, file = "/home/brendan/projects/HuBMAP/vignettes/20221025-pkhl.triplet.trends.subdist50.dist100.res50-200.rds")