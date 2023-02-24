#' Test if intra sample vs inter samples trends are significantly different
#' 
#' @description Computes summed distance between z scores at each resolution for two given trends.
#'     Gets summed distances for combinations of intra patient trends, and inter patient trends.
#'     T-test between intra and inter distances to test for significant difference.
#'     Also can return heatmap matrix of distances between each dataset and distributions.
#'     The idea is that if inter trends are similar, then their distances will be small like the intra trends.
#'     But if the inter trends are different, their distances should be larger than the intra trends.
#'     Likewise, if the intra trends are highly variable, then might be expected to also have large differences between inter as well.
#' 
#' @param samples list of sample set names for each patient ex: list(c("PKHL", "XXCD"),c("KSFB", "NGPL"),c("PBVN", "FSLD"))
#' @param refID name of reference cell type for trend of interest
#' @param neighID name of neighbor cell type for trend of interest
#' @param zscores table of zscores for each cell type combo (rows) for each sample vs each resolution tested (columns)
#' @param heatmap return heatmap of distances between sample trends (boolean, default: TRUE)
#' @param distplot return plot of distirbutions of intra and inter trend distances (boolean, default: TRUE)
#'
#' @return pvalues for each patient
#' 
#' @examples
#' \dontrun{
#' diffTrendTesting(list(c("PKHL", "XXCD"),c("KSFB", "NGPL"),c("PBVN", "FSLD")),
#'                  refID = "Ki67 proliferating",
#'                  neighID = "CD4 Memory T cells",
#'                  zscores = combined_zscores)
#' }
#'
diffTrendTesting <- function(samples, refID, neighID, zscores, heatmap = TRUE, distplot = TRUE){
  
  r <- refID
  n <- neighID
  id <- paste0(r, " vs ", n)
  m <- expand.grid(rep(list(unlist(samples)),2))
  
  ## get table of summed distances for each trend combination
  distances <- unlist(lapply(1:nrow(m), function(i){
    t1 <- paste0(id, "_", as.vector(m[i,1]))
    t2 <- paste0(id, "_", as.vector(m[i,2]))
    t1_zscores <- zscores[t1,]
    t2_zscores <- zscores[t2,]
    sum(abs(t1_zscores - t2_zscores))
  }))
  m$dist <- distances
  
  ## get list of inter trend comparisons for each set of patient datasets
  s <- 1:length(samples)
  inter.pairs <- lapply(s, function(ix){
    patient <- samples[[which(s %in% ix)]]
    others <- unlist(samples[which(!s %in% ix)])
    var1 <- c()
    for(sample in patient){
      var1 <- c(var1, rep(sample, length(others)))
    }
    var2 <- rep(others, length(patient))
    inter_pairs <- cbind(var1, var2) 
  })
  
  ## vector of distances for the intra trends for each patient
  intra.vals <- unlist(lapply(samples, function(i){
    m[m$Var1 == i[1] & m$Var2 == i[2], "dist"]
  }))
  
  ## list of distances between inter trend comparisons for each set of patient datasets
  inter.vals <- lapply(1:length(inter.pairs), function(ix){
    pairs <- inter.pairs[[ix]]
    vals <- unlist(lapply(1:nrow(pairs), function(i){
      m[m$Var1 == pairs[i,1] & m$Var2 == pairs[i,2], "dist"]
    }))
    vals
  })
  
  ## get distributions of the distances
  intra.dists <- MASS::fitdistr(intra.vals, "normal")
  inter.dists <- lapply(1:length(inter.vals), function(ix){
    MASS::fitdistr(inter.vals[[ix]], "normal")
  })
  
  ## perform T-tests
  inter.tests <- lapply(1:length(inter.vals), function(ix){
    t.test(intra.vals, inter.vals[[ix]])
  })
  
  ## plotting
  if(heatmap){
    dat <- m
    plt <- ggplot2::ggplot(data = dat) + ggplot2::geom_tile(ggplot2::aes(x = Var1, y = Var2, fill = dist)) + 
      ggplot2::scale_y_discrete(breaks = as.character(dat$Var2), 
                                labels = as.character(dat$Var2)) +
      ggplot2::labs(title = id,
                    x = "samples",
                    y = "samples")
    
    plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, 
                                                                    color = "black",
                                                                    hjust = 0.5,
                                                                    vjust = 0.5),
                                axis.text.y = ggplot2::element_text(size = 12,
                                                                    color = "black"),
                                axis.title.y = ggplot2::element_text(size = 13),
                                axis.title.x = ggplot2::element_text(size = 13),
                                plot.title = ggplot2::element_text(size = 15),
                                legend.text = ggplot2::element_text(size = 15,
                                                                    colour = "black"),
                                legend.title = ggplot2::element_text(size = 15,
                                                                     colour = "black",
                                                                     angle = 90),
                                panel.background = ggplot2::element_blank(),
                                panel.border = ggplot2::element_rect(fill = NA,
                                                                     color = "black",
                                                                     size = 2),
                                plot.background = ggplot2::element_blank()) +
      ggplot2::geom_text(ggplot2::aes(x = as.character(Var1),
                                      y = as.character(Var2),
                                      label = format(round(dist,1),nsmall = 2))) +
      ggplot2::scale_fill_gradientn(limits = c(0,max(m$dist)),
                                    breaks = c(0,max(m$dist)),
                                    colors = (grDevices::colorRampPalette(c("white","red")))(n = 209)) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Distance",
                                                     title.position = "left",
                                                     title.hjust = 0.5,
                                                     ticks.colour = "black",
                                                     ticks.linewidth = 2,
                                                     frame.colour = "black",
                                                     frame.linewidth = 2,
                                                     label.hjust = 0)) +
      ggplot2::coord_fixed()
    print(plt)
  }
  
  ## return plot of distributions of intra and inter patient trend distances
  if(distplot){
    xlim <- max(unlist(inter.vals)) * 1.1
    xvals <- seq(0,xlim,length = 500)
    
    intra <- dnorm(xvals, mean = intra.dists$estimate[[1]], sd = intra.dists$estimate[[2]])
    
    df <- do.call(cbind.data.frame, lapply(1:length(inter.dists), function(ix){
      dist <- inter.dists[[ix]]
      dnorm(xvals, mean = dist$estimate[[1]], sd = dist$estimate[[2]])
    }))
    colnames(df) <- paste0("inter_", 1:length(inter.dists))
    df[["intra"]] <- intra
    
    dat <- reshape2::melt(df)
    dat[["xvals"]] <- rep(xvals, length(inter.dists)+1)
    
    cols <- rainbow(n = length(inter.dists))
    plt <- ggplot2::ggplot(data = dat) +
      ggplot2::geom_line(ggplot2::aes(x = xvals, y = value, color = variable)) +
      ggplot2::scale_color_manual(values = c(cols, "black")) +
      ggplot2::labs(title = id,
                    x = "distance",
                    y = "density") +
      
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
                     panel.grid.major.y =  ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(size = 1, colour = "black")
                     # legend.position="none"
      )
    print(plt)
    
  }
  
  return(inter.tests)
  
}