#'@useDynLib MFM
#'@importFrom Rcpp evalCpp
#' @importFrom methods as
#' @importFrom stats dbinom model.matrix sd terms optimize
#' @importFrom utils combn write.table
#' @importFrom GUESSFM abf.calc best.snps pp.nsnp run.bvs tag tags best.models snps read.snpmod expand.tags abf2snpmod snpprior
#' @importFrom parallel mclapply
#' @importFrom data.table as.data.table setnames
#' @importFrom gtools smartbind
#' @importFrom fields image.plot
#' @importFrom grDevices heat.colors
#' @importFrom graphics axis box image par
NULL
