#' selection probability with 'gglasso'
#'
#' @description The selection probability with 'gglasso' is computed based on resamplings.
#' @param x     description x an input matrix of dimension nobs x nvars. Each row represents an observation, each column corresponds to a covariate.
#' @param y     description y a response variable, where 'y = 1' corresponds to the case and 'y = 0' corresponds to the control.
#' @param group description groups an integer vector indicating group sizes or as a symmetric adjacency matrix, which characterizes the grouping or graph structure of the predictors in 'x'.
#' @param ...   additional arguments that can be supplied to gglasso.
#' @param alpha the penalty mixing parameter ranges from 0 to 1, with the default value set to 0.1.
#' @param psub  the proportion of subsamples used for resampling is denoted as psub, the default value is 0.5.
#' @param N.lam the number of lambda values used for resamplings is specified, with the default value set to 15.
#' @param K     the number of resamplings is a user-defined parameter, with the default value set to 100.
#' @param eps   the numerical computations, the tolerance for small numerical values is \eqn{1e-04}.
#' @param maxit the maximum number of iterations is \eqn{1e+07}.
#'
#' @importFrom mnormt rmnorm
#' @importFrom stats rnorm runif
#' @importFrom gglasso gglasso
#'
#'
#'@export

sel.gglasso <- function( x,
                         y,
                         group,
                         ...,
                         alpha=0.1,
                         psub=0.5,
                         N.lam=15,
                         K=100,
                         eps=1e-04,
                         maxit=1e+07){

  y <- ifelse(y == 0, -1, 1)
  if( all( levels( as.factor( y ) ) %in% 0:1 ) ) y = ifelse( y==0 , -1 , 1 )
  if (psub < 0.5 || psub >= 1)
    stop("the proportion of subsamples should be between 0.5 and 1.")

  grs <- sapply(group, length)
  gr.new <- rep(1:length(grs), grs)
  sz.gr <- length(gr.new)/length(unique(gr.new))
  ngr <- length(gr.new) / sz.gr

  newx <- x[,unlist(group)]
  wc <- which(y ==  1)
  wt <- which(y == -1)
  nc <- floor(length(wc) * psub)
  nt <- floor(length(wt) * psub)
  N <- min(K, choose(length(wc), nc) * choose(length(wt), nt))
  u0 <- as.integer(rnorm(1, -1000, 1000))

  for (i in 1:N) {

    set.seed(u0 + i)
    ss <- c(sample(wc, nc), sample(wt, nt))

    if (i == 1) {
      nlam <- gglasso(newx[sort(ss), ], y[sort(ss)], group=gr.new, eps=eps, maxit=maxit)$lambda
      N.lam <- min(N.lam, length(nlam))
      lam <- seq(nlam[1], min(nlam), length.out = N.lam)
      out <- matrix(0, ngr, length(lam))
    }

    strt_time <- Sys.time()
    g <- gglasso(newx[sort(ss), ], y[sort(ss)], group=gr.new , lambda = lam , eps=eps, maxit=maxit)
    if(i %% 10 == 1){
      print(Sys.time() - strt_time)
      print(i)
    }

    beta <- as.matrix(g$beta)
    wh.sig <- apply( beta , 2 , function(Beta) unique( gr.new[Beta!=0] ) )
    out0 <- sapply( wh.sig , function(x) as.numeric( 1:ngr %in% x ) )
    out <- out + out0
  }

  if (is.null(names(group)))
    rownames(out) <- paste("gr", 1:ngr, sep = "")
  else rownames(out) <- names(group)
  colnames(out) <- paste("s", 1:length(lam), sep = "")
  beta <- as.matrix(out/N)
  maxs <- apply(beta, 1, max)
  u <- order(maxs, decreasing = TRUE)
  mat <- cbind(u, maxs[u])
  rownames(mat) <- NULL
  colnames(mat) <- c("group", "sel.prob")
  qhat <- sum(out)/(length(alpha)*length(lam)*N)
  return(list( beta=beta, maxsel=mat, lambda=lam, qhat=qhat, K=N))
}


