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
├── NEWS.md
├── R
│   ├── comparing_trends.R
│   ├── data.R
│   ├── finding_trends.R
│   ├── processing_inputs.R
│   ├── processing_outputs.R
│   ├── simulations.R
│   └── visualization.R
├── README.md
├── codex
│   ├── codex_akoya_data_functions.R
│   ├── codex_annotation.Rmd
│   └── codex_data_to_csv.Rmd
├── data
│   ├── CODEX
│   │   ├── ProposedAnnotations.xlsx
│   │   ├── README.md
│   │   ├── akoya_output_paths.csv
│   │   └── metadata.xlsx
│   ├── pkhl.rda
│   ├── seq.rda
│   ├── seqfish
│   │   ├── seqfish.meta.csv
│   │   ├── seqfish.meta.csv.gz
│   │   ├── sp.seqfish.binomMat.near.subdist100.rds
│   │   ├── sp.seqfish.co_occurance_array.rds
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
│   ├── sim.rda
│   ├── slide.rda
│   ├── slideseq
│   │   ├── slideseq.puck190926_08.rctd.co_occurance_array.rds
│   │   ├── slideseq.puck190926_08.rctd.co_occurance_intervals.rds
│   │   ├── slideseq.puck190926_08.rctd.meta.csv.gz
│   │   ├── slideseqPuck.190926_08.binomMat.near.subdist100.rds
│   │   ├── slideseqPuck.190926_08.pairwise.results.dist100.res100-6000.rds
│   │   ├── slideseqPuck.190926_08.pairwise.results.dist100.res30-6000.rds
│   │   ├── slideseqPuck.190926_08.rctd.near.binom.subdist100.dist100.results.removeDups.rds
│   │   ├── slideseqPuck.190926_08.rctd.pairwise.30-300.results.removeDups.rds
│   │   ├── slideseqPuck.190926_08.rctd.shuffled_res30-6000.rds
│   │   ├── slideseqPuck.190926_08.rctd.subsets.near.subdist100.rds
│   │   ├── slideseqPuck.190926_08.rctd.subsets.near.subdist200.rds
│   │   ├── slideseqPuck.190926_08.shuffled_res100-6000.rds
│   │   ├── slideseqPuck.190926_08.shuffled_res30-6000.rds
│   │   ├── slideseqPuck.190926_08.subset.results.dist100.res100-6000.rds
│   │   ├── slideseqPuck.190926_08.subset.results.dist100.res30-6000.rds
│   │   ├── slideseqPuck.190926_08.subsets.near.subdist100.rds
│   │   └── sp.seqfish.co_occurance_intervals.rds
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
│   ├── 2_slideseq.Rmd
│   ├── 3_spleen.Rmd
│   ├── _config.yml
│   ├── index.md
│   └── tutorial_full.Rmd
├── inst
│   └── CITATION
├── man
│   ├── binomialTestMatrix.Rd
│   ├── diffTrendTesting.Rd
│   ├── evaluateSignificance.Rd
│   ├── filterCells.Rd
│   ├── filterChangeTrends.Rd
│   ├── filterCoTrends.Rd
│   ├── filterSepTrends.Rd
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
│   ├── seuratToSP.Rd
│   ├── sim.Rd
│   ├── simulate_background.Rd
│   ├── simulate_circles.Rd
│   ├── slide.Rd
│   ├── spToDF.Rd
│   ├── toSP.Rd
│   ├── vizAllClusters.Rd
│   ├── vizEachCluster.Rd
│   ├── vizTrends.Rd
│   └── vizTrends.heatmap.Rd
├── misc
│   ├── compare_ripleys.R
│   ├── compare_ripleys.Rmd
│   ├── misc.R
│   └── xenium_init.R
├── multiscale_celltype_colocalization_analysis.Rproj
├── plots
│   ├── seqfish
│   │   ├── 2B_endothelium_tissue_grid.pdf
│   │   ├── 2B_mesoderm_tissue_grid.pdf
│   │   ├── 2_seqfish_endothelium_tissue.pdf
│   │   ├── 2_seqfish_endothelium_tissue_legend.pdf
│   │   ├── 2_seqfish_endothelium_tissue_zoom.pdf
│   │   ├── 2_seqfish_mesoderm_tissue.pdf
│   │   ├── 2_seqfish_mesoderm_tissue_legend.pdf
│   │   ├── 2_seqfish_mesoderm_tissue_zoom.pdf
│   │   ├── 2_seqfish_trends_permutations.pdf
│   │   ├── 2_seqfish_trends_permutations_legend.pdf
│   │   ├── 2_slideseq_trends_heatmap.pdf
│   │   ├── S3_seqfish_endo_cooccurance.pdf
│   │   ├── S3_seqfish_endo_ripleys.pdf
│   │   ├── S3_seqfish_haem_cooccurance.pdf
│   │   ├── S3_seqfish_haem_ripleys.pdf
│   │   ├── S3_seqfish_intMeso_cooccurance.pdf
│   │   ├── S3_seqfish_intMeso_ripleys.pdf
│   │   ├── S3_seqfish_latMeso_cooccurance.pdf
│   │   ├── S3_seqfish_latplate_ripleys.pdf
│   │   └── S3_seqfish_merge_ripleys.pdf
│   ├── sim
│   │   ├── 1A_sim_tissue_legend.pdf
│   │   ├── 1B_sim_tissue_trends.pdf
│   │   └── 1B_sim_tissue_trends_heatmap.pdf
│   ├── slideseq
│   │   ├── 2A_slideseq_purk_berg_tissue_grid.pdf
│   │   ├── 2_slideseq_bergmann_trends.pdf
│   │   ├── 2_slideseq_mli1_trends.pdf
│   │   ├── 2_slideseq_mli2_trends.pdf
│   │   ├── 2_slideseq_purk_berg_mli_trends_heatmap.pdf
│   │   ├── 2_slideseq_purk_berg_tissue.pdf
│   │   ├── 2_slideseq_purk_berg_tissue_legend.pdf
│   │   ├── 2_slideseq_purk_berg_tissue_zoom.pdf
│   │   ├── 2_slideseq_purkinje_trends.pdf
│   │   ├── S1_slideseq_bergmann_trends_all.pdf
│   │   ├── S1_slideseq_mli1_trends_all.pdf
│   │   ├── S1_slideseq_mli2_trends_all.pdf
│   │   ├── S1_slideseq_purk_berg_mli_trends_heatmap_all.pdf
│   │   ├── S1_slideseq_purk_berg_tissue_sep.pdf
│   │   ├── S1_slideseq_purk_berg_tissue_sep_legend.pdf
│   │   ├── S1_slideseq_purkinje_trends_all.pdf
│   │   ├── S2_slideseq_berg_cooccurance.pdf
│   │   ├── S2_slideseq_berg_ripleys.pdf
│   │   ├── S2_slideseq_mli1_cooccurance.pdf
│   │   ├── S2_slideseq_mli1_ripleys.pdf
│   │   ├── S2_slideseq_mli2_cooccurance.pdf
│   │   ├── S2_slideseq_mli2_ripleys.pdf
│   │   ├── S2_slideseq_purk_cooccurance.pdf
│   │   └── S2_slideseq_purk_ripleys.pdf
│   ├── spleen
│   │   ├── 3B_spleen_pkhl_all_trends_heatmap.pdf
│   │   ├── 3B_spleen_pkhl_grid_tissue.pdf
│   │   ├── 3_spleen_pkhl_cd4_near_folB_trends.pdf
│   │   ├── 3_spleen_pkhl_cd4_near_neutro_trends.pdf
│   │   ├── 3_spleen_pkhl_cd4_subset_trends_heatmap.pdf
│   │   ├── 3_spleen_pkhl_cd4_subsets_tissue.pdf
│   │   ├── 3_spleen_pkhl_cd4_subsets_tissue_legend.pdf
│   │   ├── 3_spleen_pkhl_cd4_trends.pdf
│   │   ├── 3_spleen_pkhl_cd4_trends_heatmap.pdf
│   │   ├── 3_spleen_pkhl_grid_38_39_tissue_zoom.pdf
│   │   ├── 3_spleen_pkhl_grid_66_77_tissue_zoom.pdf
│   │   ├── 3_spleen_pkhl_tissue.pdf
│   │   ├── 3_spleen_pkhl_tissue_legend.pdf
│   │   ├── S11_spleen_fsld_tissue.pdf
│   │   ├── S11_spleen_ksfb_tissue.pdf
│   │   ├── S11_spleen_ngpl_tissue.pdf
│   │   ├── S11_spleen_pbvn_tissue.pdf
│   │   ├── S11_spleen_xxcd_tissue.pdf
│   │   ├── S12_spleen_cd4_near_folB_trends_fsld.pdf
│   │   ├── S12_spleen_cd4_near_folB_trends_heatmap.pdf
│   │   ├── S12_spleen_cd4_near_folB_trends_ksfb.pdf
│   │   ├── S12_spleen_cd4_near_folB_trends_ngpl.pdf
│   │   ├── S12_spleen_cd4_near_folB_trends_pbvn.pdf
│   │   ├── S12_spleen_cd4_near_folB_trends_xxcd.pdf
│   │   ├── S12_spleen_cd4_near_neutro_trends_fsld.pdf
│   │   ├── S12_spleen_cd4_near_neutro_trends_heatmap.pdf
│   │   ├── S12_spleen_cd4_near_neutro_trends_ksfb.pdf
│   │   ├── S12_spleen_cd4_near_neutro_trends_ngpl.pdf
│   │   ├── S12_spleen_cd4_near_neutro_trends_pbvn.pdf
│   │   ├── S12_spleen_cd4_near_neutro_trends_xxcd.pdf
│   │   ├── S12_spleen_cd4_trends_heatmap.pdf
│   │   ├── S12_spleen_fsld_cd4_folB_podo_trends.pdf
│   │   ├── S12_spleen_ksfb_cd4_folB_podo_trends.pdf
│   │   ├── S12_spleen_ngpl_cd4_folB_podo_trends.pdf
│   │   ├── S12_spleen_pbvn_cd4_folB_podo_trends.pdf
│   │   ├── S12_spleen_xxcd_cd4_folB_podo_trends.pdf
│   │   ├── S4_spleen_pkhl_tissue_all.pdf
│   │   └── S4_spleen_pkhl_tissue_all_legend.pdf
│   ├── squidpy
│   │   ├── seqfish.nhood.r100.pdf
│   │   └── slideseq.nhood.r100.png
│   └── supplementary
│       ├── cluster_corr_heatmap_KSFB_vs_NGPL.pdf
│       ├── cluster_corr_heatmap_PBVN_vs_FSLD.pdf
│       ├── cluster_corr_heatmap_PKHL_vs_KSFB.pdf
│       ├── cluster_corr_heatmap_PKHL_vs_PBVN.pdf
│       ├── cluster_corr_heatmap_PKHL_vs_XXCD.pdf
│       ├── cluster_exp_FSLD_annot.pdf
│       ├── cluster_exp_KSFB_annot.pdf
│       ├── cluster_exp_NGPL_annot.pdf
│       ├── cluster_exp_PBVN_annot.pdf
│       ├── cluster_exp_PKHL_annot.pdf
│       ├── cluster_exp_PKHL_clusters.pdf
│       ├── cluster_exp_XXCD_annot.pdf
│       ├── lda_FSLD.predictedBy.PBVN_harmonizedTo_FSLD_posterior.pdf
│       ├── lda_KSFB.predictedBy.PKHL_posterior.pdf
│       ├── lda_NGPL.predictedBy.KSFB_harmonizedTo_NGPL_posterior.pdf
│       ├── lda_PBVN.predictedBy.PKHL_posterior.pdf
│       ├── lda_XXCD.predictedBy.PKHL_posterior.pdf
│       ├── spleen_tissue_FSLD_annots.pdf
│       ├── spleen_tissue_KSFB_annots.pdf
│       ├── spleen_tissue_NGPL_annots.pdf
│       ├── spleen_tissue_PBVN_annots.pdf
│       ├── spleen_tissue_PKHL_annots.pdf
│       ├── spleen_tissue_PKHL_clusters.pdf
│       ├── spleen_tissue_XXCD_annots.pdf
│       ├── spleen_tsne_datasets_FSLD_annots_beforeHarmonyTsne.pdf
│       ├── spleen_tsne_datasets_KSFB_annots_beforeHarmonyTsne.pdf
│       ├── spleen_tsne_datasets_NGPL_annots_beforeHarmonyTsne.pdf
│       ├── spleen_tsne_datasets_PBVN_annots_beforeHarmonyTsne.pdf
│       ├── spleen_tsne_datasets_PKHL_annots_beforeHarmonyTsne.pdf
│       ├── spleen_tsne_datasets_XXCD_annots_beforeHarmonyTsne.pdf
│       ├── spleen_tsne_datasets_all_beforeHarmonyTsne.pdf
│       └── spleen_tsne_datasets_separate_beforeHarmonyTsne.pdf
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
│   ├── squidpy_co_occurance_in_R.Rmd
│   ├── squidpy_imc.ipynb
│   ├── squidpy_pkhl.ipynb
│   └── squidpy_slideseq.ipynb
├── tests
└── vignettes
    └── tutorial.Rmd
```


