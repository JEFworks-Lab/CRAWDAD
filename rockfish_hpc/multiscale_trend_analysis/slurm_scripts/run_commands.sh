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


# ----------------------------------------------------
## slide seq cerebellum data
# ----------------------------------------------------

# sbatch ./run_pairwise.sh slideseqPuck.190926_08.rctd class slideseqPuck.190926_08.rctd.shuffled_res30-6000.rds slideseqPuck.190926_08.rctd.pairwise.30-300.results.rds


# make subsets, subsetdist 100
# sbatch ./run_makeSubsets.sh slideseqPuck.190926_08.rctd class near 100 slideseqPuck.190926_08.rctd.subsets.near.subdist100.rds


## remove dups
# dist 100, subsetdist 100
# sbatch ./run_triplet.sh slideseqPuck.190926_08.rctd class slideseqPuck.190926_08.rctd.shuffled_res50-3000.rds near 100 slideseqPuck.190926_08.rctd.subsets.near.subdist100.rds slideseqPuck.190926_08.rctd.near.binom.subdist100.dist100.results.removeDups.rds

# dist 100, subsetdist 200
# sbatch ./run_triplet.sh slideseqPuck.190926_08.rctd class slideseqPuck.190926_08.rctd.shuffled_res50-3000.rds near 200 slideseqPuck.190926_08.rctd.subsets.near.subdist200.rds slideseqPuck.190926_08.rctd.near.binom.subdist200.dist100.results.removeDups.rds

# dist 200, subsetdist 200
# sbatch ./run_triplet.sh slideseqPuck.190926_08.rctd class slideseqPuck.190926_08.rctd.shuffled_res50-3000.rds near 200 slideseqPuck.190926_08.rctd.subsets.near.subdist200.rds slideseqPuck.190926_08.rctd.near.binom.subdist200.dist200.results.removeDups.rds


## with dups
# dist 100, subsetdist 100
# sbatch ./run_triplet.sh slideseqPuck.190926_08.rctd class slideseqPuck.190926_08.rctd.shuffled_res50-3000.rds near 100 slideseqPuck.190926_08.rctd.subsets.near.subdist100.rds slideseqPuck.190926_08.rctd.near.binom.subdist100.dist100.results.withDups.rds

# dist 100, subsetdist 200
# sbatch ./run_triplet.sh slideseqPuck.190926_08.rctd class slideseqPuck.190926_08.rctd.shuffled_res50-3000.rds near 200 slideseqPuck.190926_08.rctd.subsets.near.subdist200.rds slideseqPuck.190926_08.rctd.near.binom.subdist200.dist100.results.withDups.rds

# dist 200, subsetdist 200
# sbatch ./run_triplet.sh slideseqPuck.190926_08.rctd class slideseqPuck.190926_08.rctd.shuffled_res50-3000.rds near 200 slideseqPuck.190926_08.rctd.subsets.near.subdist200.rds slideseqPuck.190926_08.rctd.near.binom.subdist200.dist200.results.withDups.rds




# sbatch ./run_pairwise.sh slideseqPuck.190926_08.rctd class slideseqPuck.190926_08.rctd.shuffled_res30-6000.rds slideseqPuck.190926_08.rctd.pairwise.30-300.results.removeDups.rds

# sbatch ./run_pairwise.sh slideseqPuck.190926_08.rctd class slideseqPuck.190926_08.rctd.shuffled_res30-6000.rds slideseqPuck.190926_08.rctd.pairwise.30-300.results.withDups.rds



# make subsets, subsetdist 100
# sbatch ./run_makeSubsets.sh slideseqPuck.190926_08.rctd class near 100 slideseqPuck.190926_08.rctd.subsets.near.subdist100.rds

## remove dups
# dist 100, subsetdist 100
# sbatch ./run_triplet.sh slideseqPuck.190926_08.rctd class slideseqPuck.190926_08.rctd.shuffled_res30-6000.rds near 100 slideseqPuck.190926_08.rctd.subsets.near.subdist100.rds slideseqPuck.190926_08.rctd.near.binom.subdist100.dist100.results.removeDups.rds




# ----------------------------------------------------
## squidpy datasets
# ----------------------------------------------------

# IMC, for co_occurance notebook vignette

# sbatch ./run_pairwise.sh compute_co_occurance_data cluster sp.cooccuranceData.shuffled_res30-800.rds sp.cooccuranceData.pairwise.results.10-300.rds

# sbatch ./run_triplet.sh compute_co_occurance_data cluster sp.cooccuranceData.shuffled_res30-800.rds near 100 sp.cooccuranceData.subsets.near.subdist100.rds sp.cooccuranceData.near.binom.subdist100.dist20-100.results.rds


## dists are 10,30,50,100,200,300
# sbatch ./run_pairwise.sh compute_co_occurance_data cluster sp.cooccuranceData.shuffled_res30-800.rds sp.cooccuranceData.pairwise.results.10-300.removeDups.rds

# sbatch ./run_pairwise.sh compute_co_occurance_data cluster sp.cooccuranceData.shuffled_res30-800.rds sp.cooccuranceData.pairwise.results.10-300.withDups.rds






## triplets

### make the subsets with binomial test first, because takes a long time for each one:

# sbatch ./run_makeSubsets.sh compute_co_occurance_data cluster near 100 sp.cooccuranceData.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh compute_ripley_data cluster near 100 sp.ripleyData.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh compute_co_occurance_data cluster near 200 sp.cooccuranceData.subsets.near.subdist200.rds


# sbatch ./run_makeSubsets.sh compute_ripley_data cluster near 200 sp.ripleyData.subsets.near.subdist200.rds

# sbatch ./run_makeSubsets.sh compute_ripley_data cluster near 300 sp.ripleyData.subsets.near.subdist300.rds


# dist 100, subsetdist 100; co_occurance
# sbatch ./run_triplet.sh compute_co_occurance_data cluster sp.cooccuranceData.shuffled_res50-800.rds near 100 sp.cooccuranceData.subsets.near.subdist100.rds sp.cooccuranceData.near.binom.subdist100.dist100.results.rds

# dist 100, subsetdist 200; co_occurance
# sbatch ./run_triplet.sh compute_co_occurance_data cluster sp.cooccuranceData.shuffled_res50-800.rds near 200 sp.cooccuranceData.subsets.near.subdist200.rds sp.cooccuranceData.near.binom.subdist200.dist100.results.rds

# dist 200, subsetdist 200; co_occurance
# sbatch ./run_triplet.sh compute_co_occurance_data cluster sp.cooccuranceData.shuffled_res50-800.rds near 200 sp.cooccuranceData.subsets.near.subdist200.rds sp.cooccuranceData.near.binom.subdist200.dist200.results.rds


# dist 100, subsetdist 100; ripley
# sbatch ./run_triplet.sh compute_ripley_data cluster sp.ripleyData.shuffled_res50-3000.rds near 100 sp.ripleyData.subsets.near.subdist100.rds sp.ripleyData.near.binom.subdist100.dist100.results.rds

# dist 100, subsetdist 200; ripley
# sbatch ./run_triplet.sh compute_ripley_data cluster sp.ripleyData.shuffled_res50-3000.rds near 200 sp.ripleyData.subsets.near.subdist200.rds sp.ripleyData.near.binom.subdist200.dist100.results.rds

# dist 200, subsetdist 200; ripley
# sbatch ./run_triplet.sh compute_ripley_data cluster sp.ripleyData.shuffled_res50-3000.rds near 200 sp.ripleyData.subsets.near.subdist200.rds sp.ripleyData.near.binom.subdist200.dist200.results.rds

# dist 300, subsetdist 300; ripley
# sbatch ./run_triplet.sh compute_ripley_data cluster sp.ripleyData.shuffled_res50-3000.rds near 300 sp.ripleyData.subsets.near.subdist300.rds sp.ripleyData.near.binom.subdist300.dist300.results.rds


# ripley's pairwise
# sbatch ./run_pairwise.sh compute_ripley_data cluster sp.ripleyData.shuffled_res50-3000.rds sp.ripleyData.pairwise.results.rds


# sbatch ./run_pairwise.sh compute_ripley_data cluster sp.slideseqv2.shuffled_res30-6000.rds sp.slideseqv2.pairwise.30-300.results.rds


# triplet, subset 100, multiple distances

# sbatch ./run_triplet.sh compute_ripley_data cluster sp.slideseqv2.shuffled_res30-6000.rds near 100 sp.ripleyData.subsets.near.subdist100.rds sp.slideseqv2.near.binom.subdist100.dist30-200.results.rds

# sbatch ./run_triplet.sh compute_ripley_data cluster sp.slideseqv2.shuffled_res30-6000.rds near 200 sp.ripleyData.subsets.near.subdist200.rds sp.slideseqv2.near.binom.subdist200.dist30-200.results.rds

# sbatch ./run_triplet.sh compute_ripley_data cluster sp.slideseqv2.shuffled_res30-6000.rds near 300 sp.ripleyData.subsets.near.subdist300.rds sp.slideseqv2.near.binom.subdist300.dist30-200.results.rds


# ----------------------------------------------------
## squidpy datasets
# ----------------------------------------------------
# four_i

# pairwise
# sbatch ./run_pairwise.sh four_i_data cluster sp.fouri.shuffled_res50-800.rds sp.fouri.pairwise.results.rds

# triplet

# sbatch ./run_makeSubsets.sh four_i_data cluster near 50 sp.fouri.subsets.near.subdist50.rds

# sbatch ./run_makeSubsets.sh four_i_data cluster near 50 sp.fouri.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh four_i_data cluster near 50 sp.fouri.subsets.near.subdist200.rds


# ----------------------------------------------------
## sim datasets
# ----------------------------------------------------
# 3 circs

# pairwise

# sbatch ./run_pairwise.sh sim3circ type sim3circ.shuffled_res30-3000.rds sim3circ.pairwise.results.30-300.rds

# triplets

# sbatch ./run_makeSubsets.sh sim3circ type near 100 sim3circ.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh sim3circ type near 200 sim3circ.subsets.near.subdist200.rds

# sbatch ./run_makeSubsets.sh sim3circ type near 300 sim3circ.subsets.near.subdist300.rds


# sbatch ./run_triplet.sh sim3circ type sim3circ.shuffled_res30-3000.rds near 100 sim3circ.subsets.near.subdist100.rds sim3circ.near.binom.subdist100.dist30-300.results.rds

# sbatch ./run_triplet.sh sim3circ type sim3circ.shuffled_res30-3000.rds near 200 sim3circ.subsets.near.subdist200.rds sim3circ.near.binom.subdist200.dist30-300.results.rds

# sbatch ./run_triplet.sh sim3circ type sim3circ.shuffled_res30-3000.rds near 300 sim3circ.subsets.near.subdist300.rds sim3circ.near.binom.subdist300.dist30-300.results.rds



# ----------------------------------------------------
## check remove neigh dups
# ----------------------------------------------------

# sbatch ./run_pairwise.sh slideseqPuck.190926_08.rctd class slideseqPuck.190926_08.rctd.shuffled_res30-6000.rds slideseqPuck.190926_08.rctd.pairwise.30-300.results.removeDupNeighs.rds

# sbatch ./run_pairwise.sh compute_ripley_data cluster sp.slideseqv2.shuffled_res30-6000.rds sp.slideseqv2.pairwise.30-300.results.removeDupNeighs.rds

sbatch ./run_pairwise.sh oneCircle type oneCircle.shuffled_res30-6000.rds oneCircle.pairwise.30-300.results.removeDupNeighs.rds

sbatch ./run_pairwise.sh oneCircle type oneCircle.shuffled_res30-6000.rds oneCircle.pairwise.30-300.results.keepDupNeighs.rds


# sbatch ./run_pairwise.sh PKHL celltypes pkhl.shuffled_res100-6000.rds pkhl.pairwise.100-200.results.removeDupNeighs.rds

# sbatch ./run_pairwise.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res100-6000.rds pkhl.pairwise.100-200.folBcombined.results.removeDupNeighs.rds





sbatch ./run_makeSubsets.sh oneCircle.subsets type near 100 oneCircle.subsets.near.subdist100.rds

sbatch ./run_triplet.sh oneCircle.subsets type oneCircle.shuffled_res30-6000.rds near 100 oneCircle.subsets.near.subdist100.rds

oneCircle.triplet.30-300.results.removeDupNeighs.rds


# ----------------------------------------------------
## check remove neigh dups
# ----------------------------------------------------

##pkhl

## pairwise

## dist 100 and 200
## celltypes
## res100-6000
## remove Dups
# pkhl.pairwise.50-200.results.removeDupNeighs.rds
# rename:
# pkhl.pairwise.100-200.results.res100-6000.removeDups.rds

## dist 100 and 200
## celltypes_folBcombined
## res100-6000
## remove Dups
# pkhl.pairwise.100-200.folBcombined.results.removeDupNeighs.rds
# rename:
# pkhl.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds


### redo these pairwise:

# ## dist 30,50,100 and 200
# ## celltypes_folBcombined
# ## res100,200,300,500,600,700,800,1000,1200,1500,3000
# ## Dups
# pkhl.pairwise.folBcombined.results.rds


# ## dist 30,50,100 and 200
# ## celltypes
# ## res200,300,500,600,700,800,900,1000,1200,3000
# ## Dups
# pkhl.pairwise.results.rds


## dist to 100,200
## removeDups = FALSE
# sbatch ./run_pairwise.sh PKHL celltypes pkhl.shuffled_res100-6000.rds pkhl.pairwise.100-200.results.res100-6000.withDups.rds

## dist to 100,200
## removeDups = FALSE
# sbatch ./run_pairwise.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res100-6000.rds pkhl.pairwise.100-200.folBcombined.results.res100-6000.withDups.rds



## triplets


## subsets don't need to have dups removed because evaluates one cell at a time

# sbatch ./run_makeSubsets.sh PKHL celltypes_folBcombined near 100 pkhl.folBcombined.subsets.near.subdist100.removeDups.rds

# sbatch ./run_makeSubsets.sh PKHL celltypes near 100 pkhl.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh PKHL celltypes near 200 pkhl.subsets.near.subdist200.rds



### celltypes_folBcombined
## dist = 100, subsetdist = 100
# sbatch ./run_triplet.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res100-6000.rds near 100 pkhl.folBcombined.subsets.near.subdist100.rds pkhl.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds

# change the remove dups in the function to FALSE
# sbatch ./run_triplet.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res100-6000.rds near 100 pkhl.folBcombined.subsets.near.subdist100.rds pkhl.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.withDups.rds

## dist = 200,  subsetdist = 100
# change the dist in triplet_analysis.r to 200
# change the remove dups in the function to TRUE
# sbatch ./run_triplet.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res100-6000.rds near 100 pkhl.folBcombined.subsets.near.subdist100.rds pkhl.triplet.near.binom.subdist100.dist200.folBcombined.results.res100-6000.removeDups.rds

# change the remove dups in the function to FALSE
# sbatch ./run_triplet.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res100-6000.rds near 100 pkhl.folBcombined.subsets.near.subdist100.rds pkhl.triplet.near.binom.subdist100.dist200.folBcombined.results.res100-6000.withDups.rds

## dist = 100,  subsetdist = 200
# change the dist in triplet_analysis.r to 100
# change the remove dups in the function to TRUE
# sbatch ./run_triplet.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res100-6000.rds near 200 pkhl.folBcombined.subsets.near.subdist200.rds pkhl.triplet.near.binom.subdist200.dist100.folBcombined.results.res100-6000.removeDups.rds

# change the remove dups in the function to FALSE
# sbatch ./run_triplet.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res100-6000.rds near 200 pkhl.folBcombined.subsets.near.subdist200.rds pkhl.triplet.near.binom.subdist200.dist100.folBcombined.results.res100-6000.withDups.rds

## dist = 200,  subsetdist = 200
# change the dist in triplet_analysis.r to 200
# sbatch ./run_triplet.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res100-6000.rds near 200 pkhl.folBcombined.subsets.near.subdist200.rds pkhl.triplet.near.binom.subdist200.dist200.folBcombined.results.res100-6000.removeDups.rds

# change the remove dups in the function to FALSE
# sbatch ./run_triplet.sh PKHL celltypes_folBcombined pkhl.folBcombined.shuffled_res100-6000.rds near 200 pkhl.folBcombined.subsets.near.subdist200.rds pkhl.triplet.near.binom.subdist200.dist200.folBcombined.results.res100-6000.withDups.rds


### celltypes
## dist = 100, subsetdist = 100
# sbatch ./run_triplet.sh PKHL celltypes pkhl.shuffled_res100-6000.rds near 100 pkhl.subsets.near.subdist100.rds pkhl.triplet.near.binom.subdist100.dist100.results.res100-6000.removeDups.rds

# change the remove dups in the function to FALSE
# sbatch ./run_triplet.sh PKHL celltypes pkhl.shuffled_res100-6000.rds near 100 pkhl.subsets.near.subdist100.rds pkhl.triplet.near.binom.subdist100.dist100.results.res100-6000.withDups.rds


## dist = 200,  subsetdist = 100
# change the dist in triplet_analysis.r to 200
# sbatch ./run_triplet.sh PKHL celltypes pkhl.shuffled_res100-6000.rds near 100 pkhl.subsets.near.subdist100.rds pkhl.triplet.near.binom.subdist100.dist200.results.res100-6000.removeDups.rds

# change the remove dups in the function to FALSE
# sbatch ./run_triplet.sh PKHL celltypes pkhl.shuffled_res100-6000.rds near 100 pkhl.subsets.near.subdist100.rds pkhl.triplet.near.binom.subdist100.dist200.results.res100-6000.withDups.rds


## dist = 100,  subsetdist = 200
# change the dist in triplet_analysis.r to 100
# sbatch ./run_triplet.sh PKHL celltypes pkhl.shuffled_res100-6000.rds near 200 pkhl.subsets.near.subdist200.rds pkhl.triplet.near.binom.subdist200.dist100.results.res100-6000.removeDups.rds

# change the remove dups in the function to FALSE
# sbatch ./run_triplet.sh PKHL celltypes pkhl.shuffled_res100-6000.rds near 200 pkhl.subsets.near.subdist200.rds pkhl.triplet.near.binom.subdist200.dist100.results.res100-6000.withDups.rds


## dist = 200,  subsetdist = 200
# change the dist in triplet_analysis.r to 200
# sbatch ./run_triplet.sh PKHL celltypes pkhl.shuffled_res100-6000.rds near 200 pkhl.subsets.near.subdist200.rds pkhl.triplet.near.binom.subdist200.dist200.results.res100-6000.removeDups.rds

# change the remove dups in the function to FALSE
# sbatch ./run_triplet.sh PKHL celltypes pkhl.shuffled_res100-6000.rds near 200 pkhl.subsets.near.subdist200.rds pkhl.triplet.near.binom.subdist200.dist200.results.res100-6000.withDups.rds


## other spleen datasets, pairwise and triplets
## just removal of duplicates
## dist of 100 looks good

## memory issues with this one, 150Gb, at 300dist
## with 100 and 200 -> 5 min and 19 min, but still 70Gb memory with 14 cores
# sbatch ./run_pairwise.sh XXCD celltypes_folBcombined xxcd.folBcombined.shuffled_res100-6000.rds xxcd.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds

## memory issues with this one, 150Gb, at 300dist
# sbatch ./run_pairwise.sh NGPL celltypes_folBcombined ngpl.folBcombined.shuffled_res100-6000.rds ngpl.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds

# sbatch ./run_pairwise.sh KSFB celltypes_folBcombined ksfb.folBcombined.shuffled_res100-6000.rds ksfb.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds

# sbatch ./run_pairwise.sh PBVN celltypes_folBcombined pbvn.folBcombined.shuffled_res100-6000.rds pbvn.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds

## memory issues with this one, 150Gb, at 300dist
# sbatch ./run_pairwise.sh FSLD celltypes_folBcombined fsld.folBcombined.shuffled_res100-6000.rds fsld.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds



## triplets actually faster and more memory efficient (like 14gb). I think this is because with the subsetting the neighbors and reference sets are a lot smaller so a lot less memory for each iteration. With the full, some cell types have tens of thousands of cells, so makes the data a lot bigger. Also above for pairwise, I accidently ran using distances of 50,100,200,300. It seems to crash for NGPL and FSLD at 300, which makes sense. Tens of thousands of cells and all of their neighbors out to 300 microns are huge connectivity matrices


# sbatch ./run_triplet.sh XXCD celltypes_folBcombined xxcd.folBcombined.shuffled_res100-6000.rds near 100 xxcd.folBcombined.subsets.near.subdist100.rds xxcd.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds

# sbatch ./run_triplet.sh NGPL celltypes_folBcombined ngpl.folBcombined.shuffled_res100-6000.rds near 100 ngpl.folBcombined.subsets.near.subdist100.rds ngpl.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds

# sbatch ./run_triplet.sh KSFB celltypes_folBcombined ksfb.folBcombined.shuffled_res100-6000.rds near 100 ksfb.folBcombined.subsets.near.subdist100.rds ksfb.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds

# sbatch ./run_triplet.sh PBVN celltypes_folBcombined pbvn.folBcombined.shuffled_res100-6000.rds near 100 pbvn.folBcombined.subsets.near.subdist100.rds pbvn.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds

# sbatch ./run_triplet.sh FSLD celltypes_folBcombined fsld.folBcombined.shuffled_res100-6000.rds near 100 fsld.folBcombined.subsets.near.subdist100.rds fsld.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds



# ----------------------------------------------------
## seqfish
# ----------------------------------------------------

## pairwise

# sbatch ./run_pairwise.sh seqfish cluster sp.seqfish.shuffled_res100-6000.rds sp.seqfish.pairwise.50-300.results.res100-6000.removeDups.rds


## triplet

# sbatch ./run_makeSubsets.sh seqfish cluster near 100 sp.seqfish.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh seqfish cluster near 200 sp.seqfish.subsets.near.subdist200.rds


# sbatch ./run_triplet.sh seqfish cluster sp.seqfish.shuffled_res100-6000.rds near 100 sp.seqfish.subsets.near.subdist100.rds sp.seqfish.triplet.near.binom.subdist100.dist100.results.res100-6000.removeDups.rds

# sbatch ./run_triplet.sh seqfish cluster sp.seqfish.shuffled_res100-6000.rds near 200 sp.seqfish.subsets.near.subdist200.rds sp.seqfish.triplet.near.binom.subdist200.dist100.results.res100-6000.removeDups.rds



# ----------------------------------------------------
## simulations
# ----------------------------------------------------

# sbatch ./run_pairwise.sh sim.pairwise type sim.pairwise.shuffled_res100-6000.rds sim.pairwise.50-200.results.res100-6000.removeDups.rds



# sbatch ./run_pairwise.sh sim.subsets type sim.subsets.shuffled_res100-6000.rds sim.subsets.pairwise.50-200.results.res100-6000.removeDups.rds

# sbatch ./run_makeSubsets.sh sim.subsets type near 50 sim.subsets.near.subdist50.rds

# sbatch ./run_makeSubsets.sh sim.subsets type near 100 sim.subsets.near.subdist100.rds

# sbatch ./run_makeSubsets.sh sim.subsets type near 200 sim.subsets.near.subdist200.rds

# sbatch ./run_makeSubsets.sh sim.subsets type near 300 sim.subsets.near.subdist300.rds

## dist 50, subdist50
# sbatch ./run_triplet.sh sim.subsets type sim.subsets.shuffled_res100-6000.rds near 50 sim.subsets.near.subdist50.rds sim.subsets.triplet.near.binom.subdist50.dist50.results.res100-6000.removeDups.rds

## dist 100, subdist100
# sbatch ./run_triplet.sh sim.subsets type sim.subsets.shuffled_res100-6000.rds near 100 sim.subsets.near.subdist100.rds sim.subsets.triplet.near.binom.subdist100.dist100.results.res100-6000.removeDups.rds

## dist 200, subdist200
sbatch ./run_triplet.sh sim.subsets type sim.subsets.shuffled_res100-6000.rds near 200 sim.subsets.near.subdist200.rds sim.subsets.triplet.near.binom.subdist200.dist200.results.res100-6000.removeDups.rds

## dist 300, subdist300
# sbatch ./run_triplet.sh sim.subsets type sim.subsets.shuffled_res100-6000.rds near 300 sim.subsets.near.subdist300.rds sim.subsets.triplet.near.binom.subdist300.dist300.results.res100-6000.removeDups.rds