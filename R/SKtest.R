#' @name SKtest
#' @title High dimensional MANOVA
#' @description Testing the equality of multi-sample high dimensional mean vectors using the testing procedure by Srivastava and Kubokawa (2013). For the two-sample test, see also Srivastava and Du (2008).
#'
#' @param X.list list of \eqn{g (\ge 2)} matrices from the samples.
#'
#' @return value of test statistic and p-value.
#'
#' @references Srivastava, M. and Du, M. (2008). A test for the mean vector with fewer observations than the dimension.  \emph{Journal of Multivariate Analysis}, 99:386–402.
#' @references Srivastava, M. and Kubokawa, T. (2013). Tests for multivariate analysis of variance in high dimension under non-normality.  \emph{Journal of Multivariate Analysis}, 115:204–216.
#'
#' @export
#'
#'
#' @examples
#' # generating a list of 3-sample data
#' set.seed(123)
#' X.list <- rdatalist(g = 3)
#' str(X.list)
#' SKtest(X.list)
#'


SKtest <-function(X.list){
  g <- length(X.list)
  pvec <- unlist(lapply(X.list, ncol))
  if (all(pvec != mean(pvec))){stop("dimensions are not equal!")}
  p <- pvec[1]
  nvec <- unlist(lapply(X.list, nrow))
  n <- sum(nvec)

  X.mean <- lapply(X.list, colMeans)
  X.diff <- (mapply(function(X){(X-X.mean[[g]])}, X.mean))[, -g]
  if(g == 2){
    B <- X.diff^2/sum(1/nvec)
  } else{
    A <- (diag(1/nvec[-g]) + matrix(1/nvec[g], (g-1), (g-1)))
    B <- apply(X.diff, 1, function(x)t(x)%*%solve(A)%*%x)
  }

  X.cov <-  lapply(X.list, stats::cov)
  nS <- matrix(0, p, p)
  for(i in 1:g){
    nS = nS + X.cov[[i]]*(nvec[i]-1)
  }
  S.pooled <- nS/(n-g)
  Ds <- diag(S.pooled)
  SK.stat <- sum(B/Ds)-(n-g)*p*(g-1)/(n-g-2)

  R <- stats::cov2cor(S.pooled)
  sigsq <- 2*(1+sum(R^2)/(p^(3/2)))*(g-1)*(sum(R^2)-p^2/(n-g))

  Zstat <- SK.stat/sqrt(sigsq)
  pval <- stats::pnorm(Zstat, lower.tail=FALSE)

  return(list(statistic=Zstat, p.value=pval))
}




################################################################################
## Note: when g = 2, the follwing test gives the same result:
# set.seed(123)
# X.list <- rdatalist(g = 2)
# SKtest(X.list)
# highmean::apval_Sri2008(X.list[[1]], X.list[[2]])
# SHT::mean2.2008SD(X.list[[1]], X.list[[2]])
# pnorm(highD2pop::SK.test(X.list[[1]], X.list[[2]])$TSvalue, lower.tail = F) #one-tailed
################################################################################



