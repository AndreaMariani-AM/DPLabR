#' Custom geom_boxplot
#' @description A Custom geom_boxplot for a \code{\link[ggplot2]{ggplot}} object. It has been computed to remove outliers
#'            if too extreme and make a plot unreadable
#' @param size An integer controlling the size of the boxes. Default is 0.8
#' @param width An integer controlling the width of the boxes. Default is 0.8
#' @param position A character specifying the position. Default is "identity"
#'
#' @importFrom ggplot2 stat_summary
#'
#' @return A geom layer
#
#' @export
geom_boxplot_noOut <- function(size = 0.8, width = 0.8, position = "identity") {
  ggplot2::stat_summary(fun.data = calc_boxplot_stat,
                        geom="boxplot",
                        position = position,
                        size = size,
                        width = width)
}


#' Function to substitute geom_boxplot and ignore outliers
#'
#' @param x A \code{\link[ggplot2]{ggplot}} object
#'
#' @importFrom stats quantile
#' @noRd
#' @return A vector of modified quantiles
#'
###------- Function to subsitute geom_boxplot and ignore outliers -------###
# https://stackoverflow.com/questions/25124895/no-outliers-in-ggplot-boxplot-with-facet-wrap
calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- stats::quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}
