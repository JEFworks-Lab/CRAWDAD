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

## full path to the shuffled data rds
shuffledDataPath <- args[3]

## full path and name of the output rds file to save (which is the list format output of the function `findTrendsv2`)
outfile <- args[4]

ncs <- args[5] ## number of cores


message(args[1])
message(args[2])
message(args[3])
message(args[4])
message(args[5])

## --------------------------------------
## run the function

distances <- c(30, 50, 100, 200)
resolutions <- c(50, 100, 200, 300, 500, 600, 700, 800, 900, 1000, 1200, 1500, 3000)
# resolutions <- c(3000)
prms <- 1 ## number of permutations
sd <- 1 ## seed


pairwise_results <- lapply(distances, function(d){
  
  # shuffledData <- paste0("/home/brendan/projects/HuBMAP/data/datasets/spleen/pkhl_shuffled_res_50-3000.rds")
  if(file.exists(shuffledDataPath)){
    
    results <- findTrendsv2(pos = meta[,c("x", "y")],
                             celltypes = meta[,3],
                             resolutions,
                             dist = d,
                             sub.type = "pairwise",
                             perms = prms,
                             ncores = ncs,
                             verbose = TRUE,
                             loadShuffleFile = shuffledDataPath,
                             seed = sd)
  } else {
    results <- findTrendsv2(pos = meta[,c("x", "y")],
                             celltypes = meta[,3],
                             resolutions,
                             dist = d,
                             sub.type = "pairwise",
                             perms = prms,
                             ncores = ncs,
                             verbose = TRUE,
                             saveShuffleFilePath = shuffledDataPath,
                             seed = sd)
  }
  # print(results)
  return(results)
  
})

names(pairwise_results) <- distances

saveRDS(object = pairwise_results, file = outfile)

# saveRDS(object = pairwise_results, file = "/home/brendan/projects/HuBMAP/vignettes/20221025-pkhl.pairwise.trends.rds")