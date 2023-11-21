#' Weighted Overlapping Group Lasso for the selection of gene sets as groups
#'
#' @description Gene set analysis aims to identify differentially expressed gene sets, often ignoring genetic network structure, which is less effective for sparse signals. Weighted Overlapping Group Lasso leverages network knowledge to identify interconnected genes, combining network-based regularization with overlapping group lasso, using l2-norm of regression coefficients for individual genes can play a role of the weight of gene sets for group selection.
#'
#' @param x     description x an input matrix of dimension nobs x nvars. Each row represents an observation, each column corresponds to a covariate.
#' @param y     description y a response variable, where 'y = 1' corresponds to the case and 'y = 0' corresponds to the control.
#' @param group description groups an integer vector indicating group sizes or as a symmetric adjacency matrix, which characterizes the grouping or graph structure of the predictors in 'x'.
#' @param adjm  the incidence matrix \code{A}: \code{A[i, j]} equals 1 when group i includes variable j.
#' @param nperm the number of permutation tests.
#' @param stra  a vector of consecutive integers is used to indicate the stratum of each observation. Each stratum must contain exactly one case and at least one control. If not specified, 'pclogit' will perform a standard logistic regression.
#' @param nfold the default number of folds is 5. While 'nfold' can go up to the sample size, it is not advisable for large datasets.
#' @param alpha the penalty mixing parameter ranges from 0 to 1, with the default value set to 0.1.
#' @param nlam  the number of lambda values, with the default value set to 100.
#' @param N.lam the number of lambda values used for resamplings is specified, with the default value set to 15.
#' @param K     the number of resamplings is a user-defined parameter, with the default value set to 100.
#' @param sgnc  regression coefficients' signs can be provided if the 'group' is specified as either a list of group sizes or an adjacency matrix.
#'
#' @importFrom mnormt rmnorm
#' @importFrom stats rnorm runif
#' @importFrom Matrix Matrix
#' @importFrom glmnet glmnet
#' @importFrom pclogit pclogit
#' @importFrom gglasso gglasso
#'
#'
#' @return A matrix includes both the order of weighted selection probabilities arranged from the highest to the lowest and the corresponding weighted selection probabilities.
#' @details More detailed information, please refer to the provided reference below.
#' @references H. Sun and S. Wang (2012) Penalized Logistic Regression for High-dimensional DNA Methylation Data with Case-Control Studies, Bioinformatics 28(10), 1368-1375
#' @references H. Sun and S. Wang (2012) Network-based Regularization for Matched Case-Control Analysis of High-dimensional DNA Methylation Data, manuscript
#' @references Yang Yi and Hui Zou (2015) A fast unified algorithm for solving group-lasso penalize learning problems, Statistics and Computing 25, 1129-1141.
#'
#' @examples
#'
#' n <- 200
#' p <- 1000
#' x <- matrix(rnorm(n*p), nrow=n, ncol=p)
#' y <- c(rep(0, n/2), rep(1, n/2))
#'
#' # a total of 10 groups, each consisting of different number of and overlapping members
#' group <- list(gr1 = c(1:31), gr2 = c(1, 17:54),
#'               gr3 = c(1, 42:61), gr4 = c(1, 47:76),
#'               gr5 = c(1, 65:92), gr6 = c(1, 78:108),
#'               gr7 = c(1, 82:125), gr8 = c(1, 94:140),
#'               gr9 = c(1, 106:143), gr10 = c(1, 118:160))
#'
#' # an adjacency matrix
#' adjm <- adjacency(x=x, group=group)
#'
#' # weighted overlapping group lasso
#' wOGL <- wOGL(x=x,y=y, group=group, adjm=adjm, nperm=10, stra=NULL,
#'              nfold=5, alpha=0.1, nlam=100, N.lam=15, K=100, sgnc=NULL)
#'
#'
#' @export
#'

wOGL <- function( x,
                  y,
                  group,
                  adjm,
                  nperm=nperm,
                  stra=NULL,
                  nfold=5,
                  alpha=0.1,
                  nlam=100,
                  N.lam=15,
                  K=100,
                  sgnc=NULL) {


  x <- as.matrix(x)
  n <- as.integer(nrow(x))
  p <- as.integer(ncol(x))
  if (!is.matrix(x) ) stop("x must be a matrix")
  if (nrow(x) != length(y)) stop("x and y must have the same length")

  grs <- sapply(group, length)
  ovpgrs <- rep(1:length(grs), grs)

  lambdas <- pclogit(x, y, stra=NULL, alpha=alpha, group=adjm, nlam=nlam)$lambda

  foldid <- sample(1:nfold, length(y), replace=T)
  mat <- matrix(NA, nrow=nfold, ncol=nlam)

  for(i in 1:nfold){

    which <- foldid == i
    x.tran <- x[!which, ]
    y.tran <- y[!which]

    x.vald <- x[which, ]
    y.vald <- y[which]

    g0 <- pclogit(x.tran, y.tran, stra=NULL, alpha=alpha, group=adjm, lambda=lambdas)

     for(j in 1:nlam){
      b0 <- (g0$b0)[j]
      b1 <- (g0$beta)[, j]
      p.pred <- apply(x.vald, 1, function(x) 1/(1+exp(-sum(x*b1)-b0)))
      y.pred <- ifelse(p.pred > 0.5, 1, 0)
      mat[i, j] <- mean(y.pred != y.vald)
     }
  }

  mats <- apply(mat, 2, mean)
  opt.lam <- lambdas[which.min(mats)]


  lam.1se <- glmnet::cv.glmnet(x, y, family="binomial", nfolds=5, alpha=0)$lambda.1se
  sgnc <- glmnet(x, y, family="binomial", lambda=lam.1se)$beta

  if (!is.null(sgnc)) {
    sgnc <- sign(sgnc)
    if (length(sgnc) != ncol(x)) stop("the length of sgnc should be equal to p.")
  }

  g1 <- pclogit(x, y, stra=NULL, alpha=alpha, group=adjm, lambda=opt.lam, sgnc=sgnc)
  g2 <- per.gglasso(x, y, group=group, nperm=nperm)
  ogl.sp <- g2$per.sp.mean


  wj <- tapply(g1$beta[unlist(group)], ovpgrs, function(x) sqrt(sum(x^2)))
  wjm <- wj/max(wj)

  wogl.sp <- wjm * ogl.sp
  w.order <- order(wogl.sp, runif(length(group)), decreasing=TRUE)
  wmat <- cbind(w.order, wogl.sp[w.order])
  colnames(wmat) <- c("groups", "wsel.prob")

  return(wogl.mat = wmat)
}





