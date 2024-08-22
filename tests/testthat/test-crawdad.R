test_that("comprehensive testing for shuffling, get results, get subsets, etc", {
  
  # data(sim)
  # ncores <- 2
  # 
  # ## convert to SP
  # cells <- toSP(pos = sim[,c("x", "y")],
  #               celltypes = sim$type)
  # expect_equal(class(cells$geometry)[1], "sfc_POINT")
  # 
  # ## Make shuffled background
  # shuffle.list <- makeShuffledCells(cells,
  #                                  scales = c(100, 200, 500, 800, 1000),
  #                                  perms = 1,
  #                                  ncores = ncores,
  #                                  seed = 1,
  #                                  verbose = TRUE)
  # expect_equal(shuffle.list$`100`$`1`[1:10], c(`1` = "D", `2` = "D", `3` = "D", `4` = "D", `5` = "A", `6` = "D", 
  #                                              `7` = "D", `8` = "A", `9` = "B", `10` = "D"))
  # ## Run pairwise analysis
  # results <- findTrends(cells,
  #                      dist = 100,
  #                      shuffle.list = shuffle.list,
  #                      ncores = ncores,
  #                      verbose = TRUE)
  # expect_gt(results$A["1000", "D"], -5.47)
  # 
  # ## melt the results into a data.frame
  # dat <- meltResultsList(results)
  # expect_equal(dat[1,], structure(list(scale = 100,
  #                                      neighbor = structure(1L,
  #                                                           levels = c("A","B", "C", "D"),
  #                                                           class = "factor"),
  #                                      Z = 0,
  #                                      reference = "A",
  #                                      id = NA),
  #                                 row.names = 1L,
  #                                 class = "data.frame"))
  # 
  # ## filter trends
  # results.coloc <- filterCoTrends(results = results, alpha = 0.05)
  # results.sep <- filterSepTrends(results = results, alpha = 0.05)
  # results.change <- filterChangeTrends(results = results, alpha = 0.05)
  # expect_equal(colnames(results.coloc$A), c("A", "B"))
  # expect_equal(colnames(results.sep$A), c("C", "D"))
  # expect_null(colnames(results.change$A))
  # 
  # ## Defining subsets
  # binomMat <- binomialTestMatrix(cells,
  #                               neigh.dist = 100,
  #                               ncores = ncores,
  #                               verbose = TRUE)
  # expect_gt(binomMat[1,"A"], 0.996)
  # 
  # subset.list <- selectSubsets(binomMat,
  #                             cells$celltypes,
  #                             sub.type = "near",
  #                             sub.thresh = 0.05,
  #                             ncores = ncores,
  #                             verbose = TRUE)
  # expect_equal(subset.list$`C_near_B`[1:5], c("25", "43", "87", "227", "264"))
  # 
  # ## make a temporary cell type annotation factor with labeling cells that are of a specific subset
  # annots_temp <- selectLabels(df = cells,
  #                              com = cells$celltypes,
  #                              subset_list = subset.list,
  #                              cellIDs = c("A", "B", "C", "D"),
  #                              subsetIDs = c("C_near_B"))
  # expect_equal(annots_temp[25:30], structure(c(`25` = 4L, `26` = 5L, `27` = 5L, `28` = 5L, `29` = 1L, 
  #                                              `30` = 1L), levels = c("A", "B", "C", "C_near_B", "D"), class = "factor"))
  # 
  # ## get neighbors
  # neighCells <- getNeighbors(cells = cells,
  #                           reference.ids = subset.list[["C_near_B"]],
  #                           removeRef = TRUE, ## whether to keep the reference cells in the output or to remove them and just look at the neighbors.
  #                           dist = 100,
  #                           returnSP = FALSE)
  # expect_equal(head(neighCells), structure(c(`1` = NA, `2` = NA, `3` = NA, `4` = NA, `5` = 1L, 
  #                                            `6` = NA), levels = c("A", "B", "C", "D"), class = "factor"))
  # 
  
})