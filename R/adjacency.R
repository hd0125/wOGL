#' adjacency matrix function
#'
#' @description Definition of the incidence matrix \eqn{A}: \code{A[i, j]} equals 1 when group i includes variable j.
#' @param x     description x an input matrix of dimension nobs x nvars. Each row represents an observation, each column corresponds to a covariate.
#' @param group description groups an integer vector indicating group sizes or as a symmetric adjacency matrix, which characterizes the grouping or graph structure of the predictors in 'x'.
#'
#' @export

adjacency <- function(x, group){

  p <- as.integer(ncol(x))
  adjm <- matrix(0, p, p)
  if (nrow(adjm) != ncol(adjm)){
    stop("adjacency matrix must be a square matrix.")
  }

  for (grp in group) {
    for (i in 2:length(grp)) {
      from.node <- grp[i - 1]
      to.node <- grp[i]
      adjm[from.node, to.node] <- 1
      adjm[to.node, from.node] <- 1
      diag(adjm) <- 0
    }
  }

  u.adjm <- sort(unique(as.numeric(adjm)))
  if (!all(u.adjm %in% c(0, 1))) {
    stop("adjacency matrix must have only 0 and 1.")
  }

  return(adjm = adjm)
}


