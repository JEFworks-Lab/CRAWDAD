#!/bin/bash
#SBATCH --job-name=multiscaleMakeSubsets
#SBATCH --error /home/bmille79/logfiles/log-%x.%j.stderr
#SBATCH --output /home/bmille79/logfiles/log-%x.%j.stdout
#SBATCH --account=jfan9
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0:4000:0

## note that the above #SBATCH lines are very specific to set flags for the sbattch command
## they do not take variables, and only have certain placeholders like %x (jobname) %j (jobid)
## check out other options and flags to  set.

## the triplets fail a lot due to memory issues.
## in some cases I am seeing memory hit +70G.
## not sure why but probably why it was failing on Easley
## Each node of Rockfish has up to 184G
## Let's alocate out to 112G for now (14 cores x 8G)
## and only doing one distance at a time
## I think it's 3.91G / core by default?

## where the meta data objects are located
dataDir="/home/bmille79/data_jfan9/brendan/multiscale_trend_analysis/"

## name of the metadata dataset, ie "PKHL"
dataset=$dataDir$1".meta.csv.gz"

## column of the celltype annotations to use, ie "celltypes", or "celltypes_folBcombined"
annots=$2

subsetType=$3

subDist=$4

subsetDataPath=$dataDir"subsetDataSets/"$5

## taken from the number of tasks/cores selected above in SBATCH options
cores=$SLURM_NTASKS_PER_NODE 


## load the appropriate modules needed to run code
## apparently need the following to get sp loaded in rockfish
## otherwise the manual install is a pain
ml r/4.0.2
ml udunits/2.2.28
ml gdal/3.4.1
##module load r-mass
##module load r-matrix
ml ## checks which modules are loaded. Printed to stderr. Note that all messages are to stderr and all print statements and output to stdout

##echo $SLURM_NTASKS_PER_NODE # this env variable should return the value of --ntasks-per-node
##echo $cores

echo Rscript --vanilla /home/bmille79/code/multiscale_trend_analysis/makeSubsets.r $dataset $annots $subsetType $subDist $subsetDataPath $cores

Rscript --vanilla /home/bmille79/code/multiscale_trend_analysis/makeSubsets.r $dataset $annots $subsetType $subDist $subsetDataPath $cores


## example running this script:

# sbatch ./run_makeSubsets.sh PKHL celltypes_folBcombined near 100 pkhl.folBcombined.subsets.near.subdist100.rds
