#!/bin/bash

## commands to run multiscale trend analysis sbatch slurm scripts


## from /home/bmille79/code/multiscale_trend_analysis/slurm_scripts/

## pairwise:

# sbatch ./run_pairwise.sh PKHL celltypes pkhl.shuffled_res50-3000.rds pkhl.pairwise.results.rds

# sbatch ./run_pairwise.sh XXCD celltypes xxcd.shuffled_res50-3000.rds xxcd.pairwise.results.rds

# sbatch ./run_pairwise.sh NGPL celltypes ngpl.shuffled_res50-3000.rds ngpl.pairwise.results.rds

# sbatch ./run_pairwise.sh KSFB celltypes ksfb.shuffled_res50-3000.rds ksfb.pairwise.results.rds

# sbatch ./run_pairwise.sh PBVN celltypes pbvn.shuffled_res50-3000.rds pbvn.pairwise.results.rds

# sbatch ./run_pairwise.sh FSLD celltypes fsld.shuffled_res50-3000.rds fsld.pairwise.results.rds

## pairwise with Fol B cells combined:

# sbatch ./run_pairwise.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res50-3000.rds pkhl.pairwise.folBcombined.results.rds

# sbatch ./run_pairwise.sh XXCD celltypes_folBcombined xxcd.folBcombined.shuffled_res50-3000.rds xxcd.pairwise.folBcombined.results.rds

# sbatch ./run_pairwise.sh NGPL celltypes_folBcombined ngpl.folBcombined.shuffled_res50-3000.rds ngpl.pairwise.folBcombined.results.rds

# sbatch ./run_pairwise.sh KSFB celltypes_folBcombined ksfb.folBcombined.shuffled_res50-3000.rds ksfb.pairwise.folBcombined.results.rds

# sbatch ./run_pairwise.sh PBVN celltypes_folBcombined pbvn.folBcombined.shuffled_res50-3000.rds pbvn.pairwise.folBcombined.results.rds

# sbatch ./run_pairwise.sh FSLD celltypes_folBcombined fsld.folBcombined.shuffled_res50-3000.rds fsld.pairwise.folBcombined.results.rds

# ----------------------------------------------------
## triplets

### make the subsets with binomial test first, because takes a long time for each one:

# sbatch ./run_makeSubsets.sh PKHL celltypes_folBcombined near 100 pkhl.folBcombined.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh XXCD celltypes_folBcombined near 100 xxcd.folBcombined.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh NGPL celltypes_folBcombined near 100 ngpl.folBcombined.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh KSFB celltypes_folBcombined near 100 ksfb.folBcombined.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh PBVN celltypes_folBcombined near 100 pbvn.folBcombined.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh FSLD celltypes_folBcombined near 100 fsld.folBcombined.subsets.near.subdist100.rds

# ---------

# sbatch ./run_makeSubsets.sh PKHL celltypes_folBcombined near 200 pkhl.folBcombined.subsets.near.subdist200.rds

# sbatch ./run_makeSubsets.sh XXCD celltypes_folBcombined near 200 xxcd.folBcombined.subsets.near.subdist200.rds

# sbatch ./run_makeSubsets.sh NGPL celltypes_folBcombined near 200 ngpl.folBcombined.subsets.near.subdist200.rds

# sbatch ./run_makeSubsets.sh KSFB celltypes_folBcombined near 200 ksfb.folBcombined.subsets.near.subdist200.rds

# sbatch ./run_makeSubsets.sh PBVN celltypes_folBcombined near 200 pbvn.folBcombined.subsets.near.subdist200.rds

# sbatch ./run_makeSubsets.sh FSLD celltypes_folBcombined near 200 fsld.folBcombined.subsets.near.subdist200.rds

# ---------

### with binomial subsetting:
# Rscript --vanilla /home/bmille79/code/multiscale_trend_analysis/triplet_analysis.r $dataset $annots $shuffledDataPath $subsetType $subsetDist $subsetDataPath $outfile $cores

### folBcombined, near, binom test, subdist 100, dist 100

# sbatch ./run_triplet.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res50-3000.rds near 100 pkhl.folBcombined.subsets.near.subdist100.rds pkhl.triplet.near.binom.subdist100.dist100.folBcombined.results.rds

# sbatch ./run_triplet.sh XXCD celltypes_folBcombined xxcd.folBcombined.shuffled_res50-3000.rds near 100 xxcd.folBcombined.subsets.near.subdist100.rds xxcd.triplet.near.binom.subdist100.dist100.folBcombined.results.rds

# sbatch ./run_triplet.sh NGPL celltypes_folBcombined ngpl.folBcombined.shuffled_res50-3000.rds near 100 ngpl.folBcombined.subsets.near.subdist100.rds ngpl.triplet.near.binom.subdist100.dist100.folBcombined.results.rds

# sbatch ./run_triplet.sh KSFB celltypes_folBcombined ksfb.folBcombined.shuffled_res50-3000.rds near 100 ksfb.folBcombined.subsets.near.subdist100.rds ksfb.triplet.near.binom.subdist100.dist100.folBcombined.results.rds

# sbatch ./run_triplet.sh PBVN celltypes_folBcombined pbvn.folBcombined.shuffled_res50-3000.rds near 100 pbvn.folBcombined.subsets.near.subdist100.rds pbvn.triplet.near.binom.subdist100.dist100.folBcombined.results.rds

# sbatch ./run_triplet.sh FSLD celltypes_folBcombined fsld.folBcombined.shuffled_res50-3000.rds near 100 fsld.folBcombined.subsets.near.subdist100.rds fsld.triplet.near.binom.subdist100.dist100.folBcombined.results.rds

# ---------

### folBcombined, near, binom test, subdist 200, dist 100

# sbatch ./run_triplet.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res50-3000.rds near 200 pkhl.folBcombined.subsets.near.subdist200.rds pkhl.triplet.near.binom.subdist200.dist100.folBcombined.results.rds

# sbatch ./run_triplet.sh XXCD celltypes_folBcombined xxcd.folBcombined.shuffled_res50-3000.rds near 200 xxcd.folBcombined.subsets.near.subdist200.rds xxcd.triplet.near.binom.subdist200.dist100.folBcombined.results.rds

# sbatch ./run_triplet.sh NGPL celltypes_folBcombined ngpl.folBcombined.shuffled_res50-3000.rds near 200 ngpl.folBcombined.subsets.near.subdist200.rds ngpl.triplet.near.binom.subdist200.dist100.folBcombined.results.rds

# sbatch ./run_triplet.sh KSFB celltypes_folBcombined ksfb.folBcombined.shuffled_res50-3000.rds near 200 ksfb.folBcombined.subsets.near.subdist200.rds ksfb.triplet.near.binom.subdist200.dist100.folBcombined.results.rds

# sbatch ./run_triplet.sh PBVN celltypes_folBcombined pbvn.folBcombined.shuffled_res50-3000.rds near 200 pbvn.folBcombined.subsets.near.subdist200.rds pbvn.triplet.near.binom.subdist200.dist100.folBcombined.results.rds

# sbatch ./run_triplet.sh FSLD celltypes_folBcombined fsld.folBcombined.shuffled_res50-3000.rds near 200 fsld.folBcombined.subsets.near.subdist200.rds fsld.triplet.near.binom.subdist200.dist100.folBcombined.results.rds

# ---------

### folBcombined, near, binom test, subdist 200, dist 200

# sbatch ./run_triplet.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res50-3000.rds near 200 pkhl.folBcombined.subsets.near.subdist200.rds pkhl.triplet.near.binom.subdist200.dist200.folBcombined.results.rds

# sbatch ./run_triplet.sh XXCD celltypes_folBcombined xxcd.folBcombined.shuffled_res50-3000.rds near 200 xxcd.folBcombined.subsets.near.subdist200.rds xxcd.triplet.near.binom.subdist200.dist200.folBcombined.results.rds

# sbatch ./run_triplet.sh NGPL celltypes_folBcombined ngpl.folBcombined.shuffled_res50-3000.rds near 200 ngpl.folBcombined.subsets.near.subdist200.rds ngpl.triplet.near.binom.subdist200.dist200.folBcombined.results.rds

# sbatch ./run_triplet.sh KSFB celltypes_folBcombined ksfb.folBcombined.shuffled_res50-3000.rds near 200 ksfb.folBcombined.subsets.near.subdist200.rds ksfb.triplet.near.binom.subdist200.dist200.folBcombined.results.rds

# sbatch ./run_triplet.sh PBVN celltypes_folBcombined pbvn.folBcombined.shuffled_res50-3000.rds near 200 pbvn.folBcombined.subsets.near.subdist200.rds pbvn.triplet.near.binom.subdist200.dist200.folBcombined.results.rds

# sbatch ./run_triplet.sh FSLD celltypes_folBcombined fsld.folBcombined.shuffled_res50-3000.rds near 200 fsld.folBcombined.subsets.near.subdist200.rds fsld.triplet.near.binom.subdist200.dist200.folBcombined.results.rds

