#' Permutation test on selection probability with 'gglasso'
#'
#' @description The permutation test on selection probability with 'gglasso' is computed based on resamplings.
#' @param x     description x an input matrix of dimension nobs x nvars. Each row represents an observation, each column corresponds to a covariate.
#' @param y     description y a response variable, where 'y = 1' corresponds to the case and 'y = 0' corresponds to the control.
#' @param group description groups an integer vector indicating group sizes or as a symmetric adjacency matrix, which characterizes the grouping or graph structure of the predictors in 'x'.
#' @param nperm the number of permutation is a user-defined parameter, with the default value set to 100.
#'
#' @importFrom mnormt rmnorm
#' @importFrom stats rnorm runif
#' @importFrom gglasso gglasso
#'
#'
#'@export


per.gglasso <- function(x,
                        y,
                        group,
                        nperm=nperm) {

  set.seed(nperm)
  u <- NULL
  sp <- as.list( 1:nperm )

  for (per in 1:nperm) {
    fits.per <- sel.gglasso(x, sample(y), group=group, K=100, N.lam=15)$maxsel
    u <- fits.per[order(fits.per[,1], decreasing=F),2]
    sp[[per]] <- u
  }

  per.sp <- do.call("cbind", sp)
  colnames(per.sp) <- paste0("Perm", 1:nperm)
  per.sp.mean <- apply(per.sp, 1, mean)

  results <- list(per.sp = per.sp, per.sp.mean = per.sp.mean)
  return(results)


}
