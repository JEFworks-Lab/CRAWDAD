# Functions to analyse Akoya datasets from HuBMAP
# Version 0.99.1

library(MASS)
library(RANN)
library(igraph)
library(Seurat)
library(ggplot2)
library(gplots)
library(gridExtra)
library(Matrix)


#' Take Akoya output for a given CODEX dataset, extract data of interest, and convert it to a Seurat object
#' 
#' @description In the Seurat object, two Assays made:
#'     RawExpr - the raw intensity values
#'     VolnormExpr - intensity values normalized by cell volume
#'    
#'     In each assay, counts are the values, and data is the log-transformed values
#'     
#'     For the expression, it removed Blanks, Empty, and DAPI measurements.
#'     
#'     Built to be used on the compensated.csv output files from Akoya.
#' 
#' @param path path to csv file
#' @param name name of the sample, a misc tag (default: NULL)
#' @param volNorm set how to normalize using cell areas (this log transforms them first) (default: "log")
#' @param verbose verbosity (default: TRUE)
#' @param plot not used yet but perhaps show histograms of expression
#'
#' @return a Seurat object of the dataset.
#'     Meta.data includes x and y positions and cell areas. misc slot has information about the commands used.
#' 
#' @example
#' obj <- buildObject(path = /path/to/csv, name = "Tissue X", volNorm = "log")
#'
buildObject <- function(path,
                        name = NULL,
                        volNorm = "log",
                        verbose = TRUE,
                        plot = TRUE){
  
  if(verbose){
    message("Data Path:", "\n", path)
  }
  results <- read.csv2(file = path, header = TRUE, sep = ",") 
  
  # cell_ids <- results[,"cell_id.cell_id"]
  cell_ids <- rownames(results)
  
  ## get cell positions
  pos <- results[,grepl(pattern = "x.x|y.y", x = colnames(results), ignore.case = TRUE)]
  rownames(pos) <- cell_ids
  colnames(pos) <- c("x", "y")
  ## Based on the image on the website, y coordinates should be reversed.
  pos[,"y"] <- pos[,"y"] * -1
  
  ## cell areas
  cellsize <- results[,grepl(pattern = "size.size", x = colnames(results), ignore.case = TRUE)]
  names(cellsize) <- cell_ids
  
  ## raw gene expression values
  ## contains the measurements for each cycle
  ## these include each protein but also multiple dapi, blanks, and empty readings
  gexp_raw <- results[,grepl(pattern = "cyc", x = colnames(results), ignore.case = TRUE)]
  ## need to convert to numeric
  gexp_raw <- apply(gexp_raw, 2, as.numeric)
  rownames(gexp_raw) <- cell_ids
  
  ## remove the dapi, blanks, empty readings
  gexp <- gexp_raw[,!grepl(pattern = "DAPI|BLANK|EMPTY|HOECHST|HAND|DRAQ5", x = colnames(gexp_raw), ignore.case = TRUE)]
  
  ## fix names for the filtered gexp. Keep for the raw so know the channel and cycle information
  ## this removed the cyc and channel information
  colnames(gexp) <- stringr::str_split_fixed(string = colnames(gexp), pattern = "\\.", n = 2)[,2]
  ## for older codex, names have additional "." like "CD45.RO", so remove to compare to other datasets
  colnames(gexp) <- gsub("\\.", "", colnames(gexp))
  
  ## normalize expression by cell area
  gexp_volNorm <- gexp/log10(cellsize)
  
  ## metadata table for cluster community labels
  metadata <- data.frame(row.names = cell_ids,
                         "cellarea" = cellsize,
                         "x" = pos[,"x"],
                         "y" = pos[,"y"])
  
  if(verbose){
    message("Creating Seurat object with RawExpr assay...")
  }
  obj <- Seurat::CreateSeuratObject(counts = t(gexp), assay = "RawExpr", project = name, meta.data = metadata)
  Seurat::Misc(object = obj, slot = "name") <- name
  Seurat::Misc(object = obj, slot = "data_path") <- path
  Seurat::Misc(object = obj, slot = "orig_ids") <- colnames(results[,grepl(pattern = "cyc0", x = colnames(results))])
  ## obj@assays$RawExpr@counts is the original data
  ## obj@assays$RawExpr@data will be result if you do Seurat::NormalizeData
  ## but let's use this slot to hold just the log10-transformed data
  obj <- Seurat::SetAssayData(object = obj, slot = "data", assay = "RawExpr",
                              new.data = as(log10(as.matrix(obj@assays$RawExpr@counts) + 1), "dgCMatrix"))
  
  if(verbose){
    message("Adding cell volume normalized expression VolnormExpr assay...")
  }
  ## add in normalized expression
  # divide each column (cell) by it's area
  ## transformation of the area
  if(volNorm == "log"){
    message("VolnormExpr via / log10(cellarea)")
    t <- log10(obj@meta.data$cellarea)
    Seurat::Misc(object = obj, slot = "VolnormExpr") <- "VolnormExpr via / log10(cellarea)"
  } else {
    message("VolnormExpr via / cellarea")
    t <- obj@meta.data$cellarea
    Seurat::Misc(object = obj, slot = "VolnormExpr") <- "VolnormExpr via / cellarea"
  }
  m <- sweep(obj@assays$RawExpr@counts, MARGIN=2, FUN="/", STATS=t) ## margin 2 means column-wise
  norm_assay <- Seurat::CreateAssayObject(counts = m)
  # add this assay to the previously created Seurat object
  obj[["VolnormExpr"]] <- norm_assay
  ## obj@assays$VolnormExpr@counts is the original data
  ## obj@assays$VolnormExpr@data will be result if you do Seurat::NormalizeData
  ## but let's use this slot to hold just the log10-transformed data
  obj <- Seurat::SetAssayData(object = obj, slot = "data", assay = "VolnormExpr",
                              new.data = as(log10(as.matrix(obj@assays$VolnormExpr@counts) + 1), "dgCMatrix"))
  
  if(verbose){
    message("Number of cells:", "\n", nrow(obj@meta.data))
    message("Number of proteins:", "\n", nrow(obj@assays$RawExpr@counts))
    message("Protein names:")
    message(paste(rownames(obj@assays$RawExpr@counts), " "))
  }
  
  return(obj)
  
}


#' Import getClusters from MERINGUE
#'
#' Note that when the dataset is very large, the memory requirements also get very large
#' Example: 150k cells and k=50 cells required almost 48G. Why is this?
#'
getClusters <- function (pcs, k,
                         method = igraph::cluster_louvain,
                         weight = FALSE,
                         verbose = TRUE,
                         details = FALSE) {
  
  if (verbose) {
    print("finding approximate nearest neighbors ...")
    print(paste0("k = ", k))
  }
  # nearest neighbors in PC space
  nn = RANN::nn2(pcs, k = k) ## KNN
  names(nn) <- c('idx', 'dists')
  
  if(weight) {
    if(verbose) {
      print('using transcriptional distance weighting')
    }
    weight <- 1/(1+ as.vector(nn$dists))
  } else {
    if(verbose) {
      print('using equal weighting')
    }
    weight <- rep(1, nrow(pcs))
  }
  
  if (verbose) {
    print("calculating clustering ...")
  }
  nn.df = data.frame(from = rep(1:nrow(nn$idx), k),
                     to = as.vector(nn$idx),
                     weight = weight
  )
  g <- igraph::graph_from_data_frame(nn.df, directed = FALSE)
  g <- igraph::simplify(g)
  km <- method(g)
  if (verbose) {
    mod <- igraph::modularity(km)
    if (mod < 0.3) {
      print("WARNING")
    }
    print(paste0("graph modularity: ", mod))
  }
  com <- km$membership
  names(com) <- rownames(pcs)
  com <- factor(com)
  if (verbose) {
    print("identifying cluster membership ...")
    print(table(com))
  }
  if (details) {
    return(list(com = com, mod = mod, g = g))
  }
  else {
    return(com)
  }
}


#' Determine cluster memberships of cells based on protein expressions
#' 
#' @description Uses MERINGUE::getClusters to get cluster memberships. Does this for a subset of cells to speed.
#'     Uses the memberships of the subset to impute the memberships of the other cells.
#'     Can use MASS::lda to build a linear model to impute memberships.
#'     Or uses RANN::nn2 to use nearest neighbors to impute memberships.
#'     
#'     This function is really used as an internal function in getCommunities()
#'     
#'     Note that if all cells used in the subset, then getClusters() will be for all cells too,
#'     and arguably more "accurate" because all cells determined this way and not imputed.
#'     Will be a lot slower though.
#' 
#' @param mat cell by protein expression matrix
#' @param K number of nearest neighbors for MERINGUE::getClusters
#' @param nsubsample number of cells to subsample
#' @param method clustering method (default:igraph::cluster_walktrap)
#' @param seed set seed primarily for subsampling
#' @param vote if TRUE, use the RANN::nn2 instead of MASS:lda to impute memberships
#' @param verbose verbosity
#'
#' @return A list that contains
#' \itemize{
#' \item all: the imputed communities for all cells
#' \item sub: the communities for cell subset from MERIGNUE::getClusters.
#' }
#' 
#' @example
#' com <- getApproxComMembership(mat, K = 50, nsubsample = nrow(mat) * 0.5)
#'
getApproxComMembership <- function(mat, K, nsubsample = nrow(mat) * 0.5, method = igraph::cluster_louvain,
                                    seed = 0, vote = FALSE, verbose = TRUE){
  
  if (verbose) {
    print(paste0("Subsampling from ", nrow(mat), " cells to ",
                 nsubsample, " ... "))
  }
  set.seed(seed)
  subsample <- sample(rownames(mat), nsubsample)
  if (verbose) {
    print("Identifying cluster membership for subsample ... ")
  }
  pcs.sub <- mat[subsample, ]
  
  com.sub <- getClusters(pcs.sub, k = K, method = method, verbose = TRUE)
  
  if (verbose) {
    print("Imputing cluster membership for rest of cells ... ")
  }
  if (vote) {
    data <- mat[subsample, ]
    query <- mat[setdiff(rownames(mat), subsample), ]
    knn <- RANN::nn2(data, query, k = K)[[1]]
    rownames(knn) <- rownames(query)
    com.nonsub <- unlist(apply(knn, 1, function(x) {
      nn <- rownames(data)[x]
      nn.com <- com.sub[nn]
      return(names(sort(table(nn.com), decreasing = TRUE)[1]))
    }))
    com.all <- factor(c(com.sub, com.nonsub)[rownames(mat)])
  }
  else {
    df.sub <- data.frame(celltype = com.sub, pcs.sub)
    model <- MASS::lda(celltype ~ ., data = df.sub)
    df.all <- data.frame(mat)
    model.output <- stats::predict(model, df.all)
    com.all <- model.output$class
    names(com.all) <- rownames(df.all)
    if (verbose) {
      print("LDA model accuracy for subsample ...")
      print(table(com.all[names(com.sub)] == com.sub))
    }
  }
  return(list("all" = com.all,
              "sub" = com.sub))
}


#' Determine cluster memberships of cells based on protein expressions and adds it as a meta.data column
#' 
#' @description Uses getApproxComMembership() to determine cell communities.
#' 
#' @param object the Seurat object
#' @param assay select the expression Assay (default: "VolnormExpr")
#' @param mat.name select the expression assay to use (counts or data) (default: "data")
#' @param numneigh number of nearest neighbors for MERINGUE::getClusters (default 50)
#' @param sample.frac porportion of cells to subsample. If 1, uses all cells and
#'    all communities are from MERINGUE::getClusters(). If less than 1,
#'    then some communities imputed by MASS::lda
#' @param method clustering method (default:igraph::cluster_walktrap)
#' @param seed set seed primarily for subsampling
#'
#' @return a Seurat object of the dataset now with the assigned communities as a column in the meta.data
#' 
#' @example
#' obj <- getCommunities(obj, assay = "VolnormExpr", mat.name = "data", numneigh = 50)
#'
getCommunities <- function(object, assay = "VolnormExpr", mat.name = "data",
                           numneigh = 50,
                           sample.frac = 1,
                           method = igraph::cluster_louvain,
                           seed = 0){
  
  ## extract matrix of interest, and transpose so cells are rows
  m <- t(as.matrix(methods::slot(object[[assay]], mat.name)))
  message("size of selected matrix: ", dim(m)[1], " ", dim(m)[2])
  
  nsub <- nrow(m) * sample.frac
  message("cells to be subsampled: ", nsub)
  
  com <- getApproxComMembership(mat = m,
                                K = numneigh,
                                nsubsample = nsub,
                                method = method,
                                seed = seed,
                                vote = FALSE,
                                verbose = TRUE)
  
  ## if getClusters done on all the cells, then com$sub are the communities computed
  ## for all the cells at once and no predictions done
  ## if using just a subset, and vote=FALSE, use MASS::lda to predict communities
  ## and want to use com$all, which are the predicted communities and here com$sub
  ## is just the getClusters for the subset
  if(sample.frac == 1){
    message("Using communities for all cells from getClusters...")
    communities <- com$sub
  } else {
    message("Using communities predicted by LDA model...")
    communities <- com$all
  }
  communities <- communities[order(as.numeric(names(communities)))]
  
  ## add a metadata column to the object
  column <- paste0("com_nn", numneigh, "_", assay, "_", mat.name)
  object[[column]] <- communities
  
  return(object)
  
}


#' Compute the protein expression for each cluster
#' 
#' @description Sums together the expression of each protein for cells in each cluster.
#'     Option to average the expression by dividing by number of cells in each cluster.
#'     Option to scale expression of each protein across clusters
#' 
#' @param object the Seurat object
#' @param clusters a column of clusters in the meta.data
#' @param assay select the expression Assay (default: "VolnormExpr")
#' @param mat select the expression assay to use ("counts" or "data") (default: "data")
#' @param ave option to average the expression
#' @param scaled option to scale expression of each protein across clusters
#' @param plot return heatmap of scaled cluster protein expression
#' 
#' @return cluster protein expression matrix.
#' 
#' @examples 
#' clusterExpressionMat <- clusterExpressions(obj, clusters = "colInMetadata", assay = "VolnormExpr", mat = "data")
#' 
clusterExpressions <- function(object, clusters,
                               assay = "VolnormExpr", mat = "data",
                               ave = TRUE,
                               scaled = TRUE,
                               plot = TRUE){
  
  c <- object@meta.data[, clusters]
  names(c) <- rownames(object@meta.data)
  
  ## extract matrix of interest, proteins are rows and cells are columns
  g <- as.matrix(methods::slot(object[[assay]], mat))
  
  mm <- model.matrix(~ 0 + c)
  rownames(mm) <- names(c)
  cluster_gexp <- g %*% mm
  
  if(ave){
    ## divide the summed expression each each protein for a cluster by the total cells in the cluster
    ## columns here are clusters, so divide each column by total cells in cluster
    cluster_gexp <- sweep(x = cluster_gexp, MARGIN = 2, STATS = table(c), FUN = "/") ## margin 2 means column-wise
  }
  
  if(scaled){
    ## scale basically centers and scales each column; so basically
    ## values in a given column are converted to z-scores so on same scale.
    ## here, each protein is basically a difference scale. So transpose
    ## so each protein is a column and then scale so each cluster compared to others
    ## on the same scale for each protein
    cluster_gexp <- scale(t(cluster_gexp))
  }
  
  if(plot){
    
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    gplots::heatmap.2(cluster_gexp,
                      density.info = "none",
                      trace = "none",
                      scale = "row", ## also scale across the rows
                      col = cm.colors(256),#sommer::jet.colors(n = 256), #cm.colors(256),
                      key.title = "",
                      key.xlab = "scaled expression",
                      cexRow = 0.8,
                      cexCol = 0.8,
                      adjRow = c(0.7,NA),
                      adjCol = c(NA,0.5),
                      offsetCol = 0
                      # margins = c(4.1, 5), # bottom margin, left margin for names
                      # lmat = rbind(4:3, 2:1), # positions in 2x2 grid where to plot stuff
                      # lhei = c(0.6, 2.0), # sizes of the heights of each section in the plotting grid
                      # lwid = c(0.7, 1.9) # sizes of the widths of each section in the plotting grid
    )
  }
  
  return(cluster_gexp)
  
}


#' Visualize expression of a given protein on the tissue for a Seurat object
#' 
#' @description uses the x and y position information and marker expression values
#'     of a selected assay in the Seurat object
#' 
#' @param object the Seurat object
#' @param marker a row of gene expression for cells in the selected matrix
#' @param assay select the expression Assay (default: "VolnormExpr")
#' @param mat select the expression assay to use ("counts" or "data") (default: "data")
#' @param title title of plot (default: NULL)
#' @param axisAdj how much to increase axis ranges. If tissue, 100 okay, if embedding, 1 ok (default: 100)
#' @param size size of points (default: 0.01)
#' @param a alpha of points (default: 1; no transparency)
#' @param legendTitle title for legend colorbar (default: NULL) 
#' 
#' @return plot of gene expression in cells
#' 
#' @examples 
#' vizExpressionObj(object = HBM389.PKHL.936,
#'     marker = "CD21",
#'     assay = "VolnormExpr", mat = "data",
#'     legendTitle = bquote(log[10]("total expression + 1")))
#' 
vizExpressionObj <- function(object, marker,
                          assay = "VolnormExpr", mat = "data",
                          axisAdj = 100, s = 0.01, a = 1, title = NULL, legendTitle = NULL){
  
  ## if object is seurat S4 object, else assume matrix and factor already
  ## will need to make this check better in future
  ## maybe have embeddings stored in object? Tricky if embedding is for mult datasets
  if(typeof(object) == "S4"){
    pos <- object@meta.data[, c("x", "y")]
    ## extract matrix of interest, proteins are rows and cells are columns
    g <- as.matrix(methods::slot(object[[assay]], mat))
    g <- g[marker,]
  } else {
    stop("object needs to be an S4 Seurat object")
  }
  
  dat <- data.frame("x" = pos[,"x"],
                    "y" = pos[,"y"],
                    "Expression" = g)
  
  plt <- ggplot2::ggplot(data = dat) +
    # ggplot2::geom_point(ggplot2::aes(x = x, y = y,
    #                                  color = Expression), size = s, alpha = a) +
    
    scattermore::geom_scattermore(ggplot2::aes(x = x, y = y,
                                               color = Expression), pointsize = s, alpha = a,
                                  pixels=c(1000,1000)) +
    
    viridis::scale_color_viridis(option = "A", direction = 1)

  plt <- plt + ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$y)-axisAdj, max(dat$y)+axisAdj)) +
               ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$x)-axisAdj, max(dat$x)+axisAdj) ) +
    
    ggplot2::labs(title = title,
                  x = "x",
                  y = "y") +
    
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                   axis.text.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 12, colour = "black"),
                   legend.title = ggplot2::element_text(size = 12, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major.y =  ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(size = 1, colour = "black")
                   # legend.position="none"
    ) +
    
    ggplot2::guides(color = ggplot2::guide_colorbar(#title = bquote(log[10]("total expression + 1")),
                                                    title = legendTitle,
                                                    title.position = "left",
                                                    title.hjust = 0.5,
                                                    ticks.colour = "black",
                                                    ticks.linewidth = 2,
                                                    frame.colour= "black",
                                                    frame.linewidth = 2,
                                                    label.hjust = 0,
                                                    barheight = 8
    )) +
    
    ggplot2::coord_equal()
  
  plt
  
}


#' Visualize expression of a given protein on the tissue for an input matrix
#' 
#' @description Needs a position matrix with x any y coordianted and an expression matrix.
#'     Rows should be cells
#' 
#' @param mat cell x marker expression matrix
#' @param pos cell x coord matrix
#' @param marker a column of gene expression for cells in the selected matrix
#' @param title title of plot (default: NULL)
#' @param axisAdj how much to increase axis ranges. If tissue, 100 okay, if embedding, 1 ok (default: 100)
#' @param size size of points (default: 0.01)
#' @param a alpha of points (default: 1; no transparency)
#' @param legendTitle title for legend colorbar (default: NULL) 
#' @param virdisCol viridis palette to use (default: "A", aka magma )
#' @param virdisDir viridis palette direction to use (default: -1, reverse direction)
#' @param customPalette argument option to input own custom color palette that will overwrite the default viridis (default: NULL)
#' 
#' @return plot of marker expression in cells
#' 
#' @examples 
#' vizExpressionMat(mat, pos, marker = "CD21")
#' 
vizExpressionMat <- function(mat, pos, marker,
                             axisAdj = 100, s = 0.01, a = 1, title = NULL, legendTitle = NULL, virdisCol = "A", virdisDir = -1, customPalette = NULL){
  
  ## make sure same cells and in same order
  cells <- intersect(rownames(pos), rownames(mat))
  
  if(length(cells) == 0){
    stop("There are no cells shared between pos and mat wrt their names", "\n",
         "Make sure that at least some rownames are shared between them.")
    
  }
  
  message(length(cells), " cells shared between pos and mat.")
  
  pos <- pos[cells,]
  mat <- mat[cells,]
  
  g <- mat[,marker]
  
  dat <- data.frame("x" = pos[,"x"],
                    "y" = pos[,"y"],
                    "Expression" = g)
  
  plt <- ggplot2::ggplot(data = dat) +
    # ggplot2::geom_point(ggplot2::aes(x = x, y = y,
    #                                  color = Expression), size = s, alpha = a) +
    scattermore::geom_scattermore(ggplot2::aes(x = x, y = y,
                                               color = Expression), pointsize = s, alpha = a,
                                  pixels=c(1000,1000))
  
  if(is.null(customPalette)){
    plt <- plt + viridis::scale_color_viridis(option = virdisCol, direction = virdisDir, na.value = "green")
  } else {
    plt <- plt + customPalette
  }
  
  plt <- plt + ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$y)-axisAdj, max(dat$y)+axisAdj)) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$x)-axisAdj, max(dat$x)+axisAdj) ) +
    
    ggplot2::labs(title = title,
                  x = "x",
                  y = "y") +
    
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                   axis.text.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 12, colour = "black"),
                   legend.title = ggplot2::element_text(size = 12, colour = "black", angle = 90),
                   legend.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major.y =  ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(size = 1, colour = "black")
                   # legend.position="none"
    ) +
    
    ggplot2::guides(color = ggplot2::guide_colorbar(#title = bquote(log[10]("total expression + 1")),
                                                    title = legendTitle,
                                                    title.position = "left",
                                                    title.hjust = 0.5,
                                                    ticks.colour = "black",
                                                    ticks.linewidth = 1,
                                                    frame.colour= "black",
                                                    frame.linewidth = 1,
                                                    label.hjust = 0,
                                                    barheight = 8
    )) +
    
    ggplot2::coord_equal()
  
  plt
  
}


transparentCol <- function(color, percent = 50, name = NULL) {
  ## Get RGB values for named color
  rgb.val <- grDevices::col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- grDevices::rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                          maxColorValue = 255,
                          alpha = (100 - percent) * 255 / 100,
                          names = name)
  
  ## Save the color
  invisible(t.col)
}

