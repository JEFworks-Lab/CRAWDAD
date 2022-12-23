#!/bin/bash
#SBATCH --job-name=multiscalePairwise
#SBATCH --error /home/bmille79/logfiles/log-%x.%j.stderr
#SBATCH --output /home/bmille79/logfiles/log-%x.%j.stdout
#SBATCH --account=jfan9
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --time=0:200:0

## note that the above #SBATCH lines are very specific to set flags for the sbattch command
## they do not take variables, and only have certain placeholders like %x (jobname) %j (jobid)
## check out other options and flags to  set.

## where the meta data objects are located
dataDir="/home/bmille79/data_jfan9/brendan/multiscale_trend_analysis/"

## name of the metadata dataset, ie "PKHL"
dataset=$dataDir$1".meta.csv.gz"

## column of the celltype annotations to use, ie "celltypes", or "celltypes_folBcombined"
annots=$2

## name and path to the shuffled data rds, ie "pkhl_shuffled_res_50-3000.rds"
## can use the same shuffled for pairwise, triplet, etc
## its just the cell labels shuffled at different resolutions
shuffledDataPath=$dataDir"shuffledDataSets/"$3

## name and path to the output rds file
outfile=$dataDir"results/pairwise/"$4

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

echo Rscript --vanilla /home/bmille79/code/multiscale_trend_analysis/pairwise_analysis.r $dataset $annots $shuffledDataPath $outfile $cores

Rscript --vanilla /home/bmille79/code/multiscale_trend_analysis/pairwise_analysis.r $dataset $annots $shuffledDataPath $outfile $cores


## example running this script:

# sbatch ./run_pairwise.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res50-3000.rds pkhl.pairwise.folBcombined.results.rds