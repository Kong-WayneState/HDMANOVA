#' @name KVstat
#' @title High dimensional MANOVA
#' @description Statistics of KVtest.
#'
#' @param X.list list of \eqn{g (\ge 2)} matrices from the samples.
#' @param method process to be used for estimating scale parameter.\describe{
#' \item{"GCT"}{generalized component test (Gregory et. al., 2015);}
#' \item{"MPT"}{more powerful test (Zhang and Wang, 2021).}
#' }
#' @param weight.fun window weight function, "parz" for Parzen weight and "trapez" for trapezoid weight.
#' @param L value of the window width. If \code{NA}, the elbow method is used for searching it.
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


KVstat <- function(X.list, method = c("MPT","GCT"), weight.fun = "parz", Lg = NA, Lm = NA){
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

  out <- c(Fn = Fn, xi_n1 =xi_n1, xi_n2 = xi_n2)

  if("GCT" %in% method){
    if(is.na(Lg)){
      gamma.acf <- as.vector(stats::acf(F.stat, lag.max = p-1, plot = FALSE, type = "covariance")$acf)
      Lg <- pathviewr::find_curve_elbow(data.frame(index = 1:p, gamma = gamma.acf)) - 1
      gct.var.est <- sum(c(1, 2*get(weight.fun)(Lg))*gamma.acf[1:(Lg+1)])
      out <- c(out, zeta_n = sqrt(gct.var.est), Lg = Lg)
    }else{
      gamma.acf <- as.vector(stats::acf(F.stat, lag.max = Lg, plot = FALSE, type = "covariance")$acf)
      gct.var.est <- sum(c(1, 2*get(weight.fun)(Lg))*gamma.acf[1:(Lg+1)])
      out <- c(out, zeta_n = sqrt(gct.var.est), Lg = Lg)
    }}

  if("MPT" %in% method){
    if(is.na(Lm)){
      gamma.cov <- CovEst(X.list = X.list, lag.max = p-1)
      Lm <- pathviewr::find_curve_elbow(data.frame(index = 1:p, gamma = gamma.cov)) - 1
      mpt.var.est <- gamma.cov[1] + 2*sum(gamma.cov[2:(Lm+1)])
      out <- c(out, tau_n = sqrt(mpt.var.est), Lm = Lm)
      }else{
      gamma.cov <- CovEst(X.list = X.list, lag.max = Lm)
      mpt.var.est <- gamma.cov[1] + 2*sum(gamma.cov[2:(Lm+1)])
      out <- c(out, tau_n = sqrt(mpt.var.est), Lm = Lm)
    }}

  return(out)

}

