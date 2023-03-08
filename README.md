# Multi-scale Cell-type Co-localization Analysis

## Project goals
- formalize computational approach to characterize cell-type co-localizations
- identify cell-type pairs, triplets,etc that are significantly co-localized 
- characterize molecular or other feature differences between co-localized vs non-co-localized cells of a particular cell-type
- enable user-friendly application to new datasets


## GitHub branch structure and files

branches:
- `brendan`
- `devel`
- `master` (default)

All branches have the following directories:
- data
- man
- R
- vignettes

Idea is that `master` (default) will eventually be submitted to `bioconductor`. They are strict in terms of the files and folders that can be included. the `data/` will only contain the  `.rda` files to allow users to load selected sample datasets. As per `bioconductor` submission rules, there will be one primary tutorial.Rmd stored in the `vignettes/` folder.

`devel` will contain `docs/`, which contains vignettes to reproduce the analyses on the different datasets. These will also be knitted to generate the pages and figures to host these analyses on the lab website. Links to these will be provided in the `README` so users can still access them on the default branch page of the GitHub repo. `devel` will also contain additional data outputs which are the results of the different analyses performed on the different datasets. Additionally, `devel` will also have `plots/`, which can be used to store the figures for the paper, which are also generated when reproducing the analyses on the different datasets

`brendan` contains records and code for slightly more exploratory analyses. For example, this branch also contains `squidpy/`, which contains `jupyter notebooks` used to compare and assess `squidpy` to `crawdad`. `rockfish_hpc/` are the scripts that were used to run some of these analyses on larger spleen datasets on the RockFish HPC. These scripts use functions that are now depricated, however, I checked the outputs using these scripts and the new functions and they are identical. Finally, `misc/` contains code and files that are related to the project but have not been incorporated.

Note that the `.gitignore` files for each branch reflect the differences in the directory structure, but for some folders, like `data/` some files were specifically removed in `devel` or `master`.

## File structures

`master` (default)

```
├── DESCRIPTION
├── NAMESPACE
├── R
│   ├── comparing_trends.R
│   ├── data.R
│   ├── finding_trends.R
│   ├── processing_inputs.R
│   ├── processing_outputs.R
│   ├── simulations.R
│   └── visualization.R
├── README.md
├── data
│   ├── pkhl.rda
│   ├── seq.rda
│   └── slide.rda
├── man
│   ├── binomialTestMatrix.Rd
│   ├── diffTrendTesting.Rd
│   ├── evaluateSignificance.Rd
│   ├── findTrends.Rd
│   ├── getNeighbors.Rd
│   ├── makeShuffledCells.Rd
│   ├── meltResultsList.Rd
│   ├── pkhl.Rd
│   ├── plotTrends.Rd
│   ├── plotTrendsOverlay.Rd
│   ├── selectLabels.Rd
│   ├── selectSubsets.Rd
│   ├── seq.Rd
│   ├── simulate_background.Rd
│   ├── simulate_circles.Rd
│   ├── slide.Rd
│   ├── spToDF.Rd
│   ├── toSP.Rd
│   ├── vizAllClusters.Rd
│   └── vizEachCluster.Rd
│   └── vizTrends.Rd
└── vignettes
    └── tutorial.Rmd
```

`devel`

```
├── DESCRIPTION
├── NAMESPACE
├── R
│   ├── comparing_trends.R
│   ├── data.R
│   ├── finding_trends.R
│   ├── processing_inputs.R
│   ├── processing_outputs.R
│   ├── simulations.R
│   └── visualization.R
├── README.md
├── data
│   ├── pkhl.rda
│   ├── seq.rda
│   ├── seqfish
│   │   ├── seqfish.meta.csv.gz
│   │   ├── sp.seqfish.binomMat.near.subdist100.rds
│   │   ├── sp.seqfish.pairwise.results.dist100.res100-6000.rds
│   │   ├── sp.seqfish.shuffled_res100-6000.rds
│   │   ├── sp.seqfish.subset.results.dist100.res100-6000.rds
│   │   └── sp.seqfish.subsets.near.subdist100.rds
│   ├── sim
│   │   └── sim.pairwise.meta.csv.gz
│   ├── slide.rda
│   ├── slideseq
│   │   ├── slideseq.puck190926_08.rctd.meta.csv.gz
│   │   ├── slideseqPuck.190926_08.binomMat.near.subdist100.rds
│   │   ├── slideseqPuck.190926_08.pairwise.results.dist100.res30-6000.rds
│   │   ├── slideseqPuck.190926_08.shuffled_res30-6000.rds
│   │   ├── slideseqPuck.190926_08.subset.results.dist100.res30-6000.rds
│   │   └── slideseqPuck.190926_08.subsets.near.subdist100.rds
│   └── spleen
│       ├── FSLD.meta.csv.gz
│       ├── KSFB
│       │   ├── ksfb.folBcombined.shuffled_res100-6000.rds
│       │   ├── ksfb.folBcombined.subsets.near.subdist100.rds
│       │   ├── ksfb.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds
│       │   └── ksfb.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds
│       ├── KSFB.meta.csv.gz
│       ├── NGPL.meta.csv.gz
│       ├── PBVN.meta.csv.gz
│       ├── PKHL
│       │   ├── pkhl.folBcombined.binomMat.near.subdist100.rds
│       │   ├── pkhl.folBcombined.pairwise.results.dist100.res100-6000.rds
│       │   ├── pkhl.folBcombined.shuffled_res100-6000.rds
│       │   ├── pkhl.folBcombined.subset.results.dist100.res100-6000.rds
│       │   └── pkhl.folBcombined.subsets.near.subdist100.rds
│       ├── PKHL.meta.csv.gz
│       └── XXCD.meta.csv.gz
├── docs
│   ├── 1_simulations.Rmd
│   ├── 2_seqfish.Rmd
│   ├── 3_slideseq.Rmd
│   ├── 4_spleen.Rmd
│   ├── _config.yml
│   └── index.md
├── man
│   ├── binomialTestMatrix.Rd
│   ├── diffTrendTesting.Rd
│   ├── evaluateSignificance.Rd
│   ├── findTrends.Rd
│   ├── getNeighbors.Rd
│   ├── makeShuffledCells.Rd
│   ├── meltResultsList.Rd
│   ├── pkhl.Rd
│   ├── plotTrends.Rd
│   ├── plotTrendsOverlay.Rd
│   ├── selectLabels.Rd
│   ├── selectSubsets.Rd
│   ├── seq.Rd
│   ├── simulate_background.Rd
│   ├── simulate_circles.Rd
│   ├── slide.Rd
│   ├── spToDF.Rd
│   ├── toSP.Rd
│   ├── vizAllClusters.Rd
│   └── vizEachCluster.Rd
│   └── vizTrends.Rd
├── plots
│   ├── seqfish
│   └── slideseq
└── vignettes
    └── tutorial.Rmd
```

`brendan`

```
├── DESCRIPTION
├── NAMESPACE
├── R
│   ├── comparing_trends.R
│   ├── data.R
│   ├── finding_trends.R
│   ├── processing_inputs.R
│   ├── processing_outputs.R
│   ├── simulations.R
│   └── visualization.R
├── README.md
├── data
│   ├── pkhl.rda
│   ├── seq.rda
│   ├── seqfish
│   │   ├── seqfish.meta.csv.gz
│   │   ├── sp.seqfish.binomMat.near.subdist100.rds
│   │   ├── sp.seqfish.pairwise.50-300.results.res100-6000.removeDups.rds
│   │   ├── sp.seqfish.pairwise.results.dist100.res100-6000.rds
│   │   ├── sp.seqfish.shuffled_res100-6000.rds
│   │   ├── sp.seqfish.subset.results.dist100.res100-6000.rds
│   │   ├── sp.seqfish.subsets.near.subdist100.rds
│   │   ├── sp.seqfish.subsets.near.subdist200.rds
│   │   ├── sp.seqfish.triplet.near.binom.subdist100.dist100.results.res100-6000.removeDups.rds
│   │   └── sp.seqfish.triplet.near.binom.subdist200.dist100.results.res100-6000.removeDups.rds
│   ├── sim
│   │   ├── sim.pairwise.50-200.results.res100-6000.removeDups.rds
│   │   ├── sim.pairwise.meta.csv.gz
│   │   ├── sim.pairwise.shuffled_res100-6000.rds
│   │   ├── sim.subsets.meta.csv.gz
│   │   ├── sim.subsets.near.subdist100.rds
│   │   ├── sim.subsets.near.subdist200.rds
│   │   ├── sim.subsets.near.subdist300.rds
│   │   ├── sim.subsets.near.subdist50.rds
│   │   ├── sim.subsets.pairwise.50-200.results.res100-6000.removeDups.rds
│   │   ├── sim.subsets.shuffled_res100-6000.rds
│   │   ├── sim.subsets.triplet.near.binom.subdist100.dist100.results.res100-6000.removeDups.rds
│   │   ├── sim.subsets.triplet.near.binom.subdist300.dist300.results.res100-6000.removeDups.rds
│   │   └── sim.subsets.triplet.near.binom.subdist50.dist50.results.res100-6000.removeDups.rds
│   ├── slide.rda
│   ├── slideseq
│   │   ├── slideseq.puck190926_08.rctd.meta.csv.gz
│   │   ├── slideseqPuck.190926_08.binomMat.near.subdist100.rds
│   │   ├── slideseqPuck.190926_08.pairwise.results.dist100.res30-6000.rds
│   │   ├── slideseqPuck.190926_08.rctd.near.binom.subdist100.dist100.results.removeDups.rds
│   │   ├── slideseqPuck.190926_08.rctd.pairwise.30-300.results.removeDups.rds
│   │   ├── slideseqPuck.190926_08.rctd.shuffled_res30-6000.rds
│   │   ├── slideseqPuck.190926_08.rctd.subsets.near.subdist100.rds
│   │   ├── slideseqPuck.190926_08.rctd.subsets.near.subdist200.rds
│   │   ├── slideseqPuck.190926_08.shuffled_res30-6000.rds
│   │   ├── slideseqPuck.190926_08.subset.results.dist100.res30-6000.rds
│   │   └── slideseqPuck.190926_08.subsets.near.subdist100.rds
│   └── spleen
│       ├── FSLD
│       │   ├── fsld.folBcombined.shuffled_res100-6000.rds
│       │   ├── fsld.folBcombined.subsets.near.subdist100.rds
│       │   ├── fsld.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds
│       │   └── fsld.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds
│       ├── FSLD.meta.csv.gz
│       ├── KSFB
│       │   ├── ksfb.folBcombined.shuffled_res100-6000.rds
│       │   ├── ksfb.folBcombined.subsets.near.subdist100.rds
│       │   ├── ksfb.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds
│       │   └── ksfb.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds
│       ├── KSFB.meta.csv.gz
│       ├── NGPL
│       │   ├── ngpl.folBcombined.shuffled_res100-6000.rds
│       │   ├── ngpl.folBcombined.subsets.near.subdist100.rds
│       │   ├── ngpl.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds
│       │   └── ngpl.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds
│       ├── NGPL.meta.csv.gz
│       ├── PBVN
│       │   ├── pbvn.folBcombined.shuffled_res100-6000.rds
│       │   ├── pbvn.folBcombined.subsets.near.subdist100.rds
│       │   ├── pbvn.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds
│       │   └── pbvn.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds
│       ├── PBVN.meta.csv.gz
│       ├── PKHL
│       │   ├── pkhl.binomMat.near.subdist100.rds
│       │   ├── pkhl.folBcombined.binomMat.near.subdist100.rds
│       │   ├── pkhl.folBcombined.pairwise.results.dist100.res100-6000.rds
│       │   ├── pkhl.folBcombined.shuffled_res100-6000.rds
│       │   ├── pkhl.folBcombined.subset.results.dist100.res100-6000.rds
│       │   ├── pkhl.folBcombined.subsets.near.subdist100.rds
│       │   ├── pkhl.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds
│       │   ├── pkhl.pairwise.100-200.results.res100-6000.removeDups.rds
│       │   ├── pkhl.subsets.near.subdist100.rds
│       │   ├── pkhl.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds
│       │   └── pkhl.triplet.near.binom.subdist100.dist100.results.res100-6000.removeDups.rds
│       ├── PKHL.meta.csv.gz
│       ├── XXCD
│       │   ├── xxcd.folBcombined.shuffled_res100-6000.rds
│       │   ├── xxcd.folBcombined.subsets.near.subdist100.rds
│       │   ├── xxcd.pairwise.100-200.folBcombined.results.res100-6000.removeDups.rds
│       │   └── xxcd.triplet.near.binom.subdist100.dist100.folBcombined.results.res100-6000.removeDups.rds
│       └── XXCD.meta.csv.gz
├── docs
│   ├── 1_simulations.Rmd
│   ├── 2_seqfish.Rmd
│   ├── 3_slideseq.Rmd
│   ├── 4_spleen.Rmd
│   ├── _config.yml
│   └── index.md
├── man
│   ├── binomialTestMatrix.Rd
│   ├── diffTrendTesting.Rd
│   ├── evaluateSignificance.Rd
│   ├── findTrends.Rd
│   ├── getNeighbors.Rd
│   ├── makeShuffledCells.Rd
│   ├── meltResultsList.Rd
│   ├── pkhl.Rd
│   ├── plotTrends.Rd
│   ├── plotTrendsOverlay.Rd
│   ├── selectLabels.Rd
│   ├── selectSubsets.Rd
│   ├── seq.Rd
│   ├── simulate_background.Rd
│   ├── simulate_circles.Rd
│   ├── slide.Rd
│   ├── spToDF.Rd
│   ├── toSP.Rd
│   ├── vizAllClusters.Rd
│   ├── vizEachCluster.Rd
│   └── vizTrends.Rd
├── misc
│   ├── misc.R
│   └── xenium_init.R
├── plots
│   ├── seqfish
│   │   ├── 2A_endothelium_tissue_legend.pdf
│   │   ├── 2A_endothelium_tissue_zoom.pdf
│   │   ├── 2A_mesoderm_tissue.pdf
│   │   ├── 2A_mesoderm_tissue_legend.pdf
│   │   ├── 2A_mesoderm_tissue_zoom.pdf
│   │   ├── 2B_endothelium_trends.pdf
│   │   └── 2B_mesoderm_trends.pdf
│   ├── sim
│   │   ├── sim.pairwise.circle.png
│   │   ├── sim.pairwise.d100.trends.pdf
│   │   ├── sim.subset.circle.png
│   │   ├── sim.subsets.pairwise.d100.trends.pdf
│   │   ├── sim.subsets.pairwise.d200.trends.pdf
│   │   ├── sim.subsets.pairwise.d50.trends.pdf
│   │   ├── sim.subsets.triplet.d100.sd100.trends.pdf
│   │   ├── sim.subsets.triplet.d300.sd300.trends.pdf
│   │   └── sim.subsets.triplet.d50.sd50.trends.pdf
│   ├── slideseq
│   │   ├── 3A_UBCs_tissue.pdf
│   │   ├── 3A_UBCs_tissue_zoom.pdf
│   │   ├── 3A_purk_berg_tissue.pdf
│   │   ├── 3A_purk_berg_tissue_zoom.pdf
│   │   ├── 3B_berg_trends.pdf
│   │   ├── 3B_granule_trends.pdf
│   │   ├── 3B_purk_trends.pdf
│   │   ├── 3B_ubc_trends.pdf
│   │   ├── slideseq.pairwise.d100.trends.pdf
│   │   ├── slideseq.rctd.pairwise.dist30-300.removeDups.pdf
│   │   ├── slideseq.rctd.pairwise.dist30-300.withDups.pdf
│   │   ├── slideseq.rctd.pairwise.dist30.dupsVsNoDups.pdf
│   │   ├── slideseq.rctd.pairwise.dist30.removeDups.overlay.pdf
│   │   ├── slideseqPuck.190926_08.rctd.pairwise.d100.neighDupVsRemoveDup.pdf
│   │   ├── slideseqPuck.190926_08.rctd.pairwise.d200.neighDupVsRemoveDup.pdf
│   │   ├── slideseqPuck.190926_08.rctd.pairwise.d50.neighDupVsRemoveDup.pdf
│   │   ├── slideseqPuck.190926_08.rctd.pairwise.pdf
│   │   ├── slideseqPuck.190926_08_dist100_subDist100.subsets.berg.pdf
│   │   ├── slideseqPuck.190926_08_dist100_subDist100.subsets.gran.overlay.pdf
│   │   ├── slideseqPuck.190926_08_dist100_subDist100.subsets.gran.pdf
│   │   ├── slideseqPuck.190926_08_dist100_subDist100.subsets.purk.overlay.pdf
│   │   └── slideseqPuck.190926_08_dist100_subDist100.subsets.purk.pdf
│   └── spleen
│       ├── 4A_ksfb_tissue.pdf
│       ├── 4A_pkhl_grid_38_39_tissue_zoom.pdf
│       ├── 4A_pkhl_grid_66_77_tissue_zoom.pdf
│       ├── 4A_pkhl_tissue.pdf
│       ├── 4B_cd4_subsets_trends.pdf
│       ├── 4B_cd4_trends.pdf
│       ├── 4B_folB_subsets_trends.pdf
│       ├── 4B_folB_trends.pdf
│       ├── 4B_ksfb_cd4_trends.pdf
│       ├── 4B_ksfb_folB_trends.pdf
│       ├── 4B_ksfb_podo_trends.pdf
│       ├── 4B_podo_subsets_trends.pdf
│       ├── 4B_podo_trends.pdf
│       ├── pkhl.cd4_near.d100.subd100.results.res100-6000.removeDups.pdf
│       ├── pkhl.cd4_subs.d100.subd100.removeDups.cd4_podo_folb.overlay.pdf
│       ├── pkhl.folb_near.d100.subd100.results.res100-6000.removeDups.pdf
│       ├── pkhl.folb_subs.d100.subd100.removeDups.cd4_podo_folb.overlay.pdf
│       ├── pkhl.pairwise.100-200.results.res100-6000.removeDups.pdf
│       ├── pkhl.pairwise.100-200.results.res100-6000.withDups.pdf
│       ├── pkhl.pairwise.d100.removeDups.overlay.pdf
│       ├── pkhl.pairwise.d100.results.res100-6000.dupsVsNoDups.pdf
│       ├── pkhl.pairwise.d100.trends.pdf
│       ├── pkhl.podo_near.d100.subd100.results.res100-6000.removeDups.pdf
│       ├── pkhl.podo_subs.d100.sd100.trends.pdf
│       └── pkhl.podo_subs.d100.subd100.removeDups.cd4_podo_folb.overlay.pdf
├── rockfish_hpc
│   └── multiscale_trend_analysis
│       ├── makeSubsets.r
│       ├── pairwise_analysis.r
│       ├── slurm_scripts
│       │   ├── run_commands.sh
│       │   ├── run_makeSubsets.sh
│       │   ├── run_pairwise.sh
│       │   └── run_triplet.sh
│       └── triplet_analysis.r
├── squidpy
│   ├── compute_co_occurrence.ipynb
│   ├── compute_ripley.ipynb
│   ├── seqfish.ipynb
│   ├── slideseq_rctd.ipynb
│   ├── squidpy_imc.ipynb
│   ├── squidpy_pkhl.ipynb
│   └── squidpy_slideseq.ipynb
└── vignettes
    └── tutorial.Rmd
```


