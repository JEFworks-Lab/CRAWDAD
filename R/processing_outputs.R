#' Melt the output list of findTrends into a dataframe
#'
#' @description The output of findTrends is a list of dataframes, where each one
#' is for a reference cell type and contains the Z scores at each resolution for the neighbor cells.
#' So here we melt this list into a single dataframe for plotting with ggplot2 and tidyverse functions.
#'
#' @param resultsList list; output from findTrends
#' @param id character; id desired, can add a column that contains an additional identifier for the results.
#' For example, you melted a resultsList that was generated from a simulation using 2 circles.
#' Then you generate another one that was done with 1 circle. The id column for each dataframe
#' can be set to "2" and "1". Then both dataframes can be combined into one final dataframe
#' Now you have identifiers that include: resolution, neighbor, reference, and simulation type.
#' Can use these for plotting and comparing different things
#'
#' @return data.frame; the resolution, neighbor, Z-score and reference for each point in the trend.
#'
#' @export
meltResultsList <- function(resultsList, id = NA){

    df <- reshape2::melt(resultsList)
    colnames(df) <- c("resolution", "neighbor", "Z", "reference")
    ## an identifier for the resultsList analysis
    df[["id"]] <- id
    return(df)

}

