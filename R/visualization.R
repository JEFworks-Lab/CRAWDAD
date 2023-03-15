#' Plot trends
#'
#' @description Plot panel of Z-score trends for each reference and neighbor cell-type pairs.
#'
#' @param results list or data.frame; the information about the resolution, Z-score, reference and the neighbor cell. It can be the result directly obtained by the findTrends function or the melted version created by the `meltResultsList()` function.
#' @param idcol character; if results are a data.frame, this is the column that contains the additional feature to plot multiple trend lines with
#' @param legend boolean to produce legend, if results are a melted data.frame with "idcol" column (default: FALSE)
#' @param ... additional plotting parameters for base R plotting. Fed into "lines()" in script
#'
#' @return nothing
#'
#' @export
plotTrends <- function(results,
                       idcol = "id",
                       legend = FALSE,
                       ...){
  
  
  ## setup to check if original list output from `findTrends`, and plot one way
  ## or if a melted dataframe with ids from merging multiple trend analyses, check if dataframe and plot that way
  
  ## if in original list format from `findTrends`:
  if(inherits(results, "list")){
    message("results detected to be a list")
    
    par(mfrow=c(length(names(results)), length(names(results))),
        mar=rep(4,4))
    
    ## for each reference cell type, ie a dataframe in the list..
    sapply(names(results), function(ct1) {
      # print(ct1)
      results.norm <- results[[ct1]]
      results.norm[is.nan(results.norm)] <- NA
      results.norm[is.infinite(results.norm)] <- NA
      
      ## for each neighbor cell type...
      sapply(colnames(results.norm), function(ct2) {
        rg <- max(abs(results.norm[, ct2]), na.rm = TRUE)
        resolutions <- rownames(results.norm)
        
        if(is.infinite(rg)){
          rg <- 2.0
        }
        
        ## instantiate a plot and plot trend
        plot(resolutions, results.norm[,ct2],
             type="l", lwd = 2,
             main=paste0(ct1,' ref \n', ct2, " neighbors"),
             ylim=c(-rg, rg),
             xlab="resolution", ylab="Z", ...)
        
        ## threshold lines
        abline(h = -2, col='red')
        abline(h = 2, col='red')
      })
    })
    
    ## if a melted dataframe,
    ## will have an additional column that can serve to plot
    ## several trend lines on the same plot instance
    ## for example, ref vs neigh at different distances
  } else if(inherits(results, "data.frame")){
    message("results detected to be a data.frame")
    
    refs <- unique(results[,"reference"])
    neighs <- unique(results[,"neighbor"])
    ids <- unique(results[,idcol])
    
    cl <- rainbow(length(ids))

    par(mfrow=c(length(neighs), length(refs)),
        mar=c(4,4,4,6)) ## bot, top, left, right
    
    ## for each reference cell type...(rows)
    sapply(refs, function(ct1) {
      # print(ct1)
      results.norm <- results[results[,"reference"] == ct1,]
      results.norm[is.nan(results.norm[,"Z"]), "Z"] <- NA
      results.norm[is.infinite(results.norm[,"Z"]), "Z"] <- NA
      
      ## for each neighbor cell type...(columns)
      sapply(neighs, function(ct2) {
        results.norm.neigh <- results.norm[results.norm[,"neighbor"] == ct2,]
        
        yl <- max(abs(results.norm.neigh[, "Z"]), na.rm = TRUE)
        xl <- max(as.numeric(results.norm.neigh[,"resolution"]))
        
        if(is.infinite(yl)){
          yl <- 2.0
        }
        
        ## instantiate a plot
        plot(0, 0, type = "n",
             main=paste0(ct1,' ref \n', ct2, " neighbors"),
             cex.main=1,
             ylim=c(-yl, yl),
             xlim=c(0, xl),
             xlab="resolution", ylab="Z")
        
        ## for each id param, draw a line on plot instance
        for(i in 1:length(ids)){
          id <- ids[i]
          
          results.norm.neigh.id <- results.norm.neigh[results.norm.neigh[,idcol] == id,]
          
          lines(as.numeric(results.norm.neigh.id[,"resolution"]), results.norm.neigh.id[,"Z"],
                type="l", lwd=2, col=cl[i], ...)
        }
        
        ## threshold lines
        abline(h = -2, col='red')
        abline(h = 2, col='red')
        
        if(legend){
          legend("topright", inset=c(-0.4,0), xpd=TRUE, legend = id, col=cl, pch=20, cex=0.5, title = idcol)
        }
        
      })
    })
    
  } else {
    stop("`results` are neither a list from `findTrends()` or a melted data.frame from `meltResultsList()`")
  }
  
}


#' This one overlays each neighbor trend wrt the same reference cell type on the plot
#' 
#' @param results data.frame; the information about the resolution, Z-score, reference and the neighbor cell. It can be the result directly obtained by the melted `findTrends` output created by the `meltResultsList()` function.
#' @param ... additional plotting parameters for base R plotting. Fed into "lines()" in script
#' 
#' @return nothing
#' 
#' @export
plotTrendsOverlay <- function(results,
                              ...){
  
  
  ## setup to check if original list output from `findTrends`, and plot one way
  ## or if a melted dataframe with ids from merging multiple trend analyses, check if dataframe and plot that way
  
  ## if in original list format from `findTrends`:
  if(inherits(results, "list")){
    message("results detected to be a list")
    
    par(mfrow=c(length(names(results)), length(names(results))),
        mar=rep(4,4))
    
    ## for each reference cell type, ie a dataframe in the list..
    sapply(names(results), function(ct1) {
      # print(ct1)
      results.norm <- results[[ct1]]
      results.norm[is.nan(results.norm)] <- NA
      results.norm[is.infinite(results.norm)] <- NA
      
      ## for each neighbor cell type...
      sapply(colnames(results.norm), function(ct2) {
        rg <- max(abs(results.norm[, ct2]), na.rm = TRUE)
        resolutions <- rownames(results.norm)
        
        ## instantiate a plot and plot trend
        plot(resolutions, results.norm[,ct2],
             type="l", lwd = 2,
             main=paste0(ct1,' ref \n', ct2, " neighbors"),
             ylim=c(-rg, rg),
             xlab="resolution", ylab="Z", ...)
        
        ## threshold lines
        abline(h = -2, col='red')
        abline(h = 2, col='red')
      })
    })
    
    ## if a melted dataframe,
    ## will have an additional column that can serve to plot
    ## several trend lines on the same plot instance
    ## for example, ref vs neigh at different distances
  } else if(inherits(results, "data.frame")){
    message("results detected to be a data.frame")
    
    results <- results[,c("resolution", "neighbor", "reference", "Z")]
    
    refs <- unique(results[,"reference"])
    neighs <- unique(results[,"neighbor"])
    cl <- rainbow(length(neighs))
    
    par(mfrow=c(length(refs),1),
        mar=c(4,4,4,8)) ## bot, top, left, right
    
    ## for each reference cell type...(rows)
    sapply(refs, function(ct1) {
      # print(ct1)
      results.norm <- results[results[,"reference"] == ct1,]
      results.norm[is.nan(results.norm[,"Z"]), "Z"] <- NA
      results.norm[is.infinite(results.norm[,"Z"]), "Z"] <- NA
      
      ## set limits based on trends with other cell types, not self
      results.norm.limits <- results.norm[results.norm[,"neighbor"] != ct1,]
      yl_max <- max(results.norm.limits[, "Z"], na.rm = TRUE)
      yl_min <- min(results.norm.limits[, "Z"], na.rm = TRUE)
      xl <- max(as.numeric(results.norm.limits[,"resolution"]))
      if(is.infinite(yl_max)){
        yl_max <- 2.0
      }
      if(is.infinite(yl_min)){
        yl_min <- -2.0
      }
      
      ## instantiate a plot
      plot(0, 0, type = "n",
           main=paste0(ct1," ref"),
           cex.main=1,
           ylim=c(yl_min,yl_max),
           xlim=c(0, xl),
           xlab="resolution", ylab="Z")
      
      ## for each neighbor cell type draw a line on plot instance
      for(i in 1:length(neighs)){
        ct2 <- neighs[i]
        # ignore showing trends with self because typically very large and
        # masks relationships with other cell types
        if(ct1 != ct2){
          results.norm.neigh.id <-  results.norm[results.norm[,"neighbor"] == ct2,]
          lines(as.numeric(results.norm.neigh.id[,"resolution"]), results.norm.neigh.id[,"Z"],
                type="l", lwd=0.8, col=cl[i], ...)
        }
      }
      
      ## threshold lines
      abline(h = -1, col='black')
      abline(h = 1, col='black')

      legend("topright", inset=c(-0.4,0), xpd=TRUE, legend = neighs, col=cl, pch=20, cex=0.5, title = "neighbors")

    })
    
  } else {
    stop("`results` are neither a list from `findTrends` or a melted data.frame from `meltResultsList`")
  }
  
}


#' Visualize all clusters on the tissue
#' 
#' @description uses the x and y position information and a chosen set of communities
#' 
#' @param cells either a data.frame or sp::SpatialPointsDataFrame object with cell spatial coordinates
#' @param coms a factor of cell type labels for the cells
#' @param ofInterest a vector of specific clusters to visualize (default; NULL)
#' @param title title of plot (default: NULL)
#' @param axisAdj how much to increase axis ranges. If tissue, 100 okay, if embedding, 1 ok (default: 100)
#' @param size size of points (default: 0.01)
#' @param a alpha of points (default: 1; no transparency)
#' @param nacol color of the NA values for cells of "other" cluster (default: (transparentCol(color = "gray", percent = 50)))
#' 
#' @return plot
#' 
#' @examples
#' data(slide)
#' vizAllClusters(slide, coms = slide$celltypes)
#' 
#' @export
vizAllClusters <- function(cells, coms, ofInterest = NULL,
                           axisAdj = 1, s = 0.5, a = 1, title = NULL,
                           nacol = transparentCol(color = "gray", percent = 50)){
  
  ## if cells are a data.frame with "x" and "y" cell coordinate columns
  if( class(cells)[1] %in% c("data.frame", "matrix") ){
    pos <- cells[,c("x", "y")]
    tempCom <- factor(coms)
    names(tempCom) <- rownames(cells)
  }
  
  ## if cells are the sp::SpatialPointsDataFrame() object
  if( any(class(cells) == "sf") ){
    p <- spToDF(cells)
    pos <- p[,c("x", "y")]
    tempCom <- factor(coms)
    names(tempCom) <- rownames(cells)
  }
  
  if(!is.null(ofInterest)){
    ## goal:
    ## setup so the clusters of interest are plotted on top of everything else
    
    tempCom[which(!tempCom %in% ofInterest)] <- NA
    tempCom <- droplevels(tempCom)
    
    cluster_cell_id <- which(tempCom %in% ofInterest)
    other_cells_id <- as.vector(which(is.na(tempCom)))
    
    cluster_cols <- rainbow(n = length(ofInterest))
    names(cluster_cols) <- ofInterest
    
    dat <- data.frame("x" = pos[,"x"],
                      "y" = pos[,"y"])
    
    ## note: "Clusters" will be a variable id used to assign colors.
    ## for the "other cells" make this NA
    dat_cluster <- data.frame("x" = pos[cluster_cell_id,"x"],
                              "y" = pos[cluster_cell_id,"y"],
                              "Clusters" = as.vector(tempCom[cluster_cell_id]))
    
    dat_other <- data.frame("x" = pos[other_cells_id,"x"],
                            "y" = pos[other_cells_id,"y"],
                            "Clusters" = NA)
    
    plt <- ggplot2::ggplot() +
      
      ## plot other cells
      scattermore::geom_scattermore(data = dat_other, ggplot2::aes(x = x, y = y,
                                                                   color = Clusters), pointsize = s, alpha = a,
                                    pixels=c(1000,1000)) +
      ## cluster cells on top
      scattermore::geom_scattermore(data = dat_cluster, ggplot2::aes(x = x, y = y,
                                                                     color = Clusters), pointsize = s, alpha = a,
                                    pixels=c(1000,1000)) +
      
      ggplot2::scale_color_manual(values = cluster_cols, na.value = nacol)
    
  } else {
    
    tempCom <- droplevels(tempCom)
    dat <- data.frame("x" = pos[,"x"],
                      "y" = pos[,"y"],
                      "Clusters" = tempCom)
    
    plt <- ggplot2::ggplot(data = dat) +

      scattermore::geom_scattermore(ggplot2::aes(x = x, y = y,
                                                 color = Clusters), pointsize = s, alpha = a,
                                    pixels=c(1000,1000)) +
      
      ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom))), na.value = nacol)
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
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(size = 1, colour = "black")
                   # legend.position="none"
    ) +
    
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
    ) +
    
    ggplot2::coord_equal()
  
  plt
  
}


#' Visualize each cluster separately
#' 
#' @description Returns a gridExtra of grobs.
#'     A single plot where each panel is a different cluster highlighted on the tissue
#' 
#' @param cells either a data.frame or sp::SpatialPointsDataFrame object with cell spatial coordinates
#' @param coms a factor of cell type labels for the cells
#' @param axisAdj how much to increase axis ranges. If tissue, 100 okay, if embedding, 1 ok (default: 100)
#' @param size size of points (default: 0.01)
#' @param a alpha of points (default: 1; no transparency)
#' @param nacol color of the NA values for cells of "other" cluster (default: (transparentCol(color = "gray", percent = 50)))
#' @param clustcol color of the cells of the given cluster (default: red)
#' 
#' @return plot where each panel is a different celltype
#' 
#' @examples 
#' data(slide)
#' vizEachCluster(slide, coms = slide$celltypes, ofInterest = c("Bergmann", "Purkinje"))
#' 
#' @export
vizEachCluster <- function(cells, coms, axisAdj = 1, s = 0.5, a = 1,
                           nacol = transparentCol(color = "gray", percent = 50),
                           clustcol = "red"){
  
  ## if cells are a data.frame or matrix  with "x" and "y" cell coordinate columns
  if( class(cells)[1] %in% c("data.frame", "matrix") ){
    pos <- cells[,c("x", "y")]
    ctemp <- factor(coms)
    names(ctemp) <- rownames(cells)
  }
  
  ## if cells are the sp::SpatialPointsDataFrame() object
  if( any(class(cells) == "sf") ){
    p <- spToDF(cells)
    pos <- p[,c("x", "y")]
    ctemp <- factor(coms)
    names(ctemp) <- rownames(cells)
  }
  
  ps <- lapply(levels(ctemp), function(c) {
    
    tempCom <- ctemp
    tempCom[which(!tempCom %in% c(c))] <- NA
    tempCom <- droplevels(tempCom)
    
    cluster_cell_id <- which(tempCom == c)
    other_cells_id <- as.vector(which(is.na(tempCom)))
    
    dat <- data.frame("x" = pos[,"x"],
                      "y" = pos[,"y"])
    
    dat_cluster <- data.frame("x" = pos[cluster_cell_id,"x"],
                              "y" = pos[cluster_cell_id,"y"],
                              "Clusters" = c)
    
    dat_other <- data.frame("x" = pos[other_cells_id,"x"],
                            "y" = pos[other_cells_id,"y"],
                            "Clusters" = NA)
    
    plt <- ggplot2::ggplot() +
      
      ## plot other cells
      scattermore::geom_scattermore(data = dat_other, ggplot2::aes(x = x, y = y,
                                                                   color = Clusters), pointsize = s, alpha = a,
                                    pixels=c(1000,1000)) +
      ## cluster cells on top
      scattermore::geom_scattermore(data = dat_cluster, ggplot2::aes(x = x, y = y,
                                                                     color = Clusters), pointsize = s, alpha = a,
                                    pixels=c(1000,1000)) +
      
      ggplot2::scale_color_manual(values = c("red"), na.value = nacol) +
      
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$y)-axisAdj, max(dat$y)+axisAdj)) +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$x)-axisAdj, max(dat$x)+axisAdj) ) +
      
      ggplot2::labs(title = c, #paste0("c", c, ", volume normalized"),
                    x = "x",
                    y = "y") +
      
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, color = "black"),
                     axis.text.y = ggplot2::element_text(size=10, color = "black"),
                     axis.title.y = ggplot2::element_text(size=10),
                     axis.title.x = ggplot2::element_text(size=10),
                     axis.ticks.x = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size=10),
                     legend.text = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     legend.background = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     panel.grid.major.y =  ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(size = 1, colour = "black"),
                     plot.margin = ggplot2::unit(c(1,1,1,1), "pt"), # change default margins around each plot
                     legend.position="none"
      ) +
      
      ggplot2::coord_equal()
    
    plt
    
  })
  
  numplots <- length(levels(ctemp))
  numcols <- 4
  numrows <- ceiling(numplots/numcols)
  numpanels <- numcols * numrows
  
  panel_layout <- t(matrix(data = 1:numpanels, nrow = numrows, ncol = numcols))
  
  p <- gridExtra::grid.arrange(
    grobs = ps,
    layout_matrix = panel_layout
  )
  
  p
}


## dat = the data.frame of pvals, trend cluster assignments, etc for each pairwise combo. i.e. PKHL
## could subset this to get specific interactions
## clusters = the name of the column that has the clusters to color by
## nc = the number of unique clusters to color by.
## default should be number of unique entries in the clusters column
## but its possible that one dataset might not have a cluster because all 
## datasets clustered together so a cluster could be specific to a particular
## datasets
## colors = by default used rainbow(nc), but can specifically change to vector 
## to be used in scale_color_manual

#' Plot trends with ggplot2
#' 
#' @description The input data.frame should be the results list from `findTrends()` that has been melted into a data.frame using `meltResultsList()`.
#' 
#' @param dat data.frame; the information about the resolution, Z-score, reference and the neighbor cell.
#' @param id column name that contains an additional feature to color trend lines (default: "id")
#' @param yaxis column that has significance value across resolutions (default: "Z")
#' @param sig.thresh threshold for significance, ie Z score significance threshold (default: 2).
#' @param nc number of colors to use for labeling features in id column
#' @param colors color assignment for each of the features in id column
#' @param title plot title (default: NULL)
#' @param facet boolean to facet wrap on reference and neighbor cell types. (default: TRUE)
#' @param lines boolean to plot lines (default: TRUE)
#' @param points boolean to plot points (default: TRUE)
#' 
#' @export
vizTrends <- function(dat, id = "id", yaxis = "Z",
                      sig.thresh = 2, # -log10(0.05/nrow(dat)), ## sig thresh for num tests
                      nc = length(unique(dat[[id]])),
                      colors = rainbow(nc),
                      title = NULL,
                      facet = TRUE,
                      lines = TRUE,
                      points = TRUE){
  
  plt <- ggplot2::ggplot(data = dat)
  if(points){
    plt <- plt + ggplot2::geom_point(ggplot2::aes(x = resolution, y = .data[[yaxis]], color = .data[[id]]), size = 1.5)
  }
  if(lines){
    plt <- plt + ggplot2::geom_path(ggplot2::aes(x = resolution, y = .data[[yaxis]], color = .data[[id]]), size = 1.5)
  }
  plt <- plt + ggplot2::scale_color_manual(values = colors, na.value = "black") +
    ggplot2::geom_hline(yintercept = 0, color = "black", size = 1) +
    ggplot2::geom_hline(yintercept = sig.thresh, color = "black", size = 0.6, linetype = "dotted") +
    ggplot2::geom_hline(yintercept = -sig.thresh, color = "black", size = 0.6, linetype = "dotted") +
    # ggplot2::facet_grid(neighbor ~ reference) +
    ggplot2::labs(title = title) +
    # ggplot2::scale_x_log10() +
    # ggplot2::scale_y_continuous(trans = ggallin::pseudolog10_trans) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, color = "black", angle = -90, vjust = 0.5, hjust = 0),
                   axis.text.y = ggplot2::element_text(size=12, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   plot.background = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_line(size = 0.1, colour = "black"),
                   panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                   axis.line = ggplot2::element_line(size = 0, colour = "black"),
                   panel.spacing = ggplot2::unit(0.1, "lines"),
                   strip.text = ggplot2::element_text(size = 12),
                   legend.title = ggplot2::element_blank(),
                   # legend.position="none"
    )
  
  if(facet){
    plt <- plt + ggplot2::facet_grid(neighbor ~ reference)
  }
  
  plt
  
}


#' Update cell type labels to only include specific cell types or subtypes
#' @description make a factor of selected cell type labels for visualizing specific cells with `vizAllClusters()`.
#'      Could append new entries into the `subset_list` to select and label other custom subsets of cells 
#'
#' @param df dataframe or sp::SpatialPointsDataFrame object of cells
#' @param com original factor or vector of cell type labels for cells in `df`
#' @param cellIDs vector of cell type labels to include in output (default: NA)
#' @param subset_list list of subsets from `selectSubsets()` (default; NA)
#' @param subsetIDs vector of susbet cell type labels to include in output (names of vectors in `subset_list`) (default: NA)
#'
#' @return factor of specific cell type labels for visualizing with `vizAllClusters()` in the parameter: `clusters`
#' 
#' @export
selectLabels <- function(df,
                         com,
                         cellIDs = NA,
                         subset_list = NA,
                         subsetIDs = NA
){
  
  ## get vector to append cell annotations of interest
  annots_temp <- rep(NA, length(rownames(df)))
  names(annots_temp) <- rownames(df)
  
  ## append the neighbor cells
  if(!is.na(cellIDs[1])){
    for(neigh_id in cellIDs){
      annots_temp[com == neigh_id] <- neigh_id
    }
  } else {
    message("`cellIDs` set to NA")
  }
  
  ## get cell ids that are part of subset
  ## these added after the neigbors in case
  ## you want to plot all CD4 T cells first
  # then color the ones that are a subset
  ## note that the order will be important
  ## because labels are overwritten in this way
  if( !is.na(subsetIDs[1]) & !is.na(subset_list[1]) ){
    for(subsetID in subsetIDs){
      cells_temp <- as.numeric(subset_list[[subsetID]])
      ## append the subset label
      annots_temp[cells_temp] <- subsetID
    }
  } else {
    message("`subsetIDs` and `subset_list` set to NA")
  }
  
  annots_temp <- as.factor(annots_temp)
  return(annots_temp)
}


#' Function to make transparent colors
#' 
#' @noRd
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