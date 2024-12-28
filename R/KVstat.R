#' @name KVstat
#' @title High dimensional MANOVA
#' @description Statistics of KVtest.
#'
#' @param X.list list of \eqn{g (\ge 2)} matrices from the samples.
#'
#' @return value of \eqn{F_n}, center estimates \eqn{\hat\xi_{n1}} and \eqn{\hat\xi_{n2}}, selected window width \eqn{L}, and scale estemates \eqn{\hat\zeta_n} and \eqn{\hat\tau_n}.
#'
#' @references Gregory, K., Carroll, R., Baladandayuthapani, V., and Lahiri, S. (2015). A two-sample test for equality of means in high dimension.  \emph{Journal of the American Statistical Association}, 110(510):837–849.
#' @references Zhang, H. and Wang, H. (2021). A more powerful test of equality of high-dimensional two-sample means.  \emph{Computational Statistics & Data Analysis}, 164:107318.
#' @references Kong, X., Villasante-Tezanos, A., Fardo, D.. and Harrar, S. (submitted). Generalized composite multi-sample tests for high-dimensional data. \emph{}
#' @export
#'
#'
#' @examples
#' # generating a list of 3-sample data
#' set.seed(123)
#' X.list <- rdatalist(g = 3)
#' str(X.list)
#' KVstat(X.list)
#'
#'


KVstat <- function(X.list){
  g <- length(X.list)
  pvec <- unlist(lapply(X.list, ncol))
  if (all(pvec != mean(pvec))){stop("dimensions are not equal!")}
  p <- pvec[1]
  nvec <- unlist(lapply(X.list, nrow))

  X.mean <- t(mapply(colMeans, X.list))
  X.var <-  t(mapply(function(X){apply(X, 2, stats::var)}, X.list))

  MSTk <- colSums((scale(X.mean, center=TRUE, scale = FALSE))^2)/(g-1)
  MSEk <- colSums(X.var/nvec)/g

  F.stat <- (MSTk/MSEk)
  Fn <- mean(F.stat)

  xi_n1 <- mean(sapply(1:p, CenterEst, X.list=X.list, npower = 1))
  xi_n2 <- mean(sapply(1:p, CenterEst, X.list=X.list, npower = 2))

  gamma.cov <- CovEst(X.list = X.list, lag.max = p-1)
  L <- pathviewr::find_curve_elbow(data.frame(index = 1:p, gamma = gamma.cov)) - 1
  zeta_n <- sum(c(1, 2*parz(L))*gamma.cov[1:(L+1)])

  gamma.acf <- as.vector(stats::acf(F.stat, lag.max = L, plot = FALSE, type = "covariance")$acf)
  tau_n <- gamma.acf[1] + 2*sum(gamma.acf[2:(L+1)])


  return(c(Fn = Fn, xi_n1 =xi_n1, xi_n2 = xi_n2, zeta_n = zeta_n, tau_n = tau_n, L = L))

}

