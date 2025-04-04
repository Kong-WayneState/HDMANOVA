#' @name KVtest
#' @title High dimensional MANOVA
#' @description Testing the equality of multi-sample high dimensional mean vectors using the testing procedure by Kong, et. al. (Submitted). For the two-sample test, see also Gregory et. al (2015) or Zhang and Wang (2021).
#'
#' @param X.list list of \eqn{g (\ge 2)} matrices from the samples.
#' @param method process to be used for estimating scale parameter.\describe{
#' \item{"GCT"}{generalized component test (Gregory et. al., 2015);}
#' \item{"MPT"}{more powerful test (Zhang and Wang, 2021).}
#' }
#' @param weight.fun window weight function, "parz" for Parzen weight and "trapez" for trapezoid weight.
#' @param L value of the window width. If \code{NA}, the elbow method is used for searching it.
#' @param npower power (negative) 0, 1, or 2 of n at which the center is estimated.
#'
#' @return value of test statistic, p-value and window width.
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
#' KVtest(X.list)
#' KVtest(X.list, method = "MPT", L = NA)
#'


KVtest <- function(X.list, method = "GCT", weight.fun = "parz", L = 10, npower = 0){
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
  KV.stat <- mean(F.stat)

  if(method == "GCT"){
    # acf is /p (general used)
    # Gregory_2015 using /(p-s) in the paper,
    # but using /p in their R function, (we use this)
    if(is.na(L)){
      gamma <- as.vector(stats::acf(F.stat, lag.max = p-1, plot = FALSE, type = "covariance")$acf)
      L <- pathviewr::find_curve_elbow(data.frame(index = 1:p, gamma = gamma)) - 1
      var.est <- sum(c(1, 2*get(weight.fun)(L))*gamma[1:(L+1)])
    }else{
      gamma <- as.vector(stats::acf(F.stat, lag.max = L, plot = FALSE, type = "covariance")$acf)
      var.est <- sum(c(1, 2*get(weight.fun)(L))*gamma)
    }
  }
  if(method == "MPT"){
    if(is.na(L)){
      gamma <- CovEst(X.list = X.list, lag.max = p-1)
      L <- pathviewr::find_curve_elbow(data.frame(index = 1:p, gamma = gamma)) - 1
      # L = 1
      # while(L < floor(2*sqrt(p))){
      #   if(gamma[L+1] > min((g-1)/sqrt(p), 1/(g-1)/sqrt(n))){
      #     L = L + 1
      #   }else break
      #}
      var.est <- gamma[1] + 2*sum(gamma[2:(L+1)])
    }else{
      gamma <- CovEst(X.list, lag.max = L)
      var.est <- gamma[1] + 2*sum(gamma[2:(L+1)])
    }
  }

  center.est <- mean(sapply(1:p, CenterEst, X.list=X.list, npower = npower))

  Zstat <- sqrt(p)*(KV.stat - center.est)/sqrt(var.est)
  pval <- stats::pnorm(Zstat, lower.tail = FALSE)  #One-tailed

  return(list(Statistic = Zstat, p.value = pval, window=L))

}

#' @keywords internal
#' @noRd
parz <- function(L){
  sapply(1:L, function(i){ifelse(i<(L/2), 1-6*(i/L)^2+6*(i/L)^3, 2*(1-(i/L))^3)})
}

#' @keywords internal
#' @noRd
trapez <- function(L){
  sapply(1:L, function(i){ifelse(i<=ceiling(L/2), 1, 1-(i-ceiling(L/2))/(L-ceiling(L/2)))})
}

#' @keywords internal
#' @noRd
CenterEst <- function(k, X.list, npower = 0){
  # norder = 0, 1, 2
  if (npower == 0) {
    return(ceter.est = 1)
  }

  x.list <- lapply(X.list, function(X)X[, k])
  g = length(x.list)
  nvec = mapply(length, x.list)
  n = sum(nvec)
  lambda <- n/nvec

  mu_unbiased <- function(xsamp){
    n <- length(xsamp)
    m1 <- mean(xsamp)
    m2 <- mean((xsamp - m1)^2)
    m3 <- mean((xsamp - m1)^3)
    m4 <- mean((xsamp - m1)^4)
    m5 <- mean((xsamp - m1)^5)

    ube_mu2 = n/(n-1)*m2
    ube_mu3 = n^2/((n-1)*(n-2))*m3
    ube_mu4 = n*(n^2-2*n+3)/((n-1)*(n-2)*(n-3))*m4 - 3*n*(2*n-3)/((n-1)*(n-2)*(n-3))*m2^2
    ube_mu5 = n^2*(n^2-5*n+10)/((n-1)*(n-2)*(n-3)*(n-4))*m5 - 10*n^2/((n-1)*(n-3)*(n-4))*m2*m3

    return(c(ube_mu2, ube_mu3, ube_mu4, ube_mu5))
  }

  mu.est <- mapply(mu_unbiased, x.list)
  tau.sq <- sum(mu.est[1, ]*lambda)

  if (npower == 1) {
    c1 <- 2*tau.sq^(-2)*sum(mu.est[1, ]^2*(lambda^3))
    c2 <- 2*tau.sq^(-3)*g*stats::var(mu.est[2, ]*(lambda^2))
    return(center.est = 1 +  n^(-1)*(c1 + c2))
  }

  if (npower == 2) {
    c1 <- 2*tau.sq^(-2)*sum(mu.est[1, ]^2*(lambda^3))
    c2 <- 2*tau.sq^(-3)*g*stats::var(mu.est[2, ]*(lambda^2))

    d1 = 2*tau.sq^(-2)*sum(mu.est[1,]^2*(lambda^4))
    d2 = 8*tau.sq^(-3)*sum((lambda^5)*(2*mu.est[1,]^3-mu.est[3,]*mu.est[1,]))
    d3 = 6*tau.sq^(-4)*(sum((mu.est[3, ]-mu.est[1, ]^2)*(lambda^3))*sum(mu.est[1, ]^2*(lambda^3)) -
                          g/(g-1)*sum((lambda^6)*(mu.est[4,]-6*mu.est[2,]*mu.est[1,])*mu.est[2,]) +
                          1/(g-1)*sum((lambda^2)*mu.est[2,])*sum((lambda^4)*(mu.est[4,]-6*mu.est[2,]*mu.est[1,])))
    d4 = tau.sq^(-5)*12*sum((mu.est[3,]-mu.est[1, ]^2)*(lambda^3))*g*stats::var(lambda^2*mu.est[2,])

    return(center.est = 1  +  n^(-1)*(c1 + c2) +  n^(-2)* (d1 + d2 + d3 + d4))
  }
}

#' @keywords internal
#' @noRd
CovEst <- function(X.list, lag.max){
  g <- length(X.list)
  p <- ncol(X.list[[1]])
  nvec <- unlist(lapply(X.list, nrow))
  n <- sum(nvec)

  lambda <- n/nvec
  X.cov <-  lapply(X.list, stats::cov)
  X.cov.lam <- Map("*", X.cov, lambda)
  nS <- Reduce("+", X.cov.lam)
  X2.cov.lam <- Map("*", X.cov.lam, X.cov.lam)
  nS2 <- Reduce("+", X2.cov.lam)

  numerator <- 2*nS2 + 2/((g-1)^2)*(nS^2-nS2)

  D <- colSums(matrix(unlist(lapply(X.cov, diag)),ncol=p, byrow=T)*lambda)
  denom <- D %*% t(D)
  est <- numerator / denom

  gamma <- sapply(0:(lag.max), function(r){sum(est * outer(1:p, 1:p, function(x, y){y - x == r}))/p})

  return(gamma)
}




################################################################################
## Note: when g = 2, the follwing test gives the same result:
# set.seed(123)
# p <- 100
# X.list <- rdatalist(g = 2, p = 100, err = "WD")
# KVtest(X.list, method = "MPT", L=ceiling(p^(3/8)))
# pnorm(highDmean::zwl_test(X.list[[1]], X.list[[2]])$statistic, lower.tail = F) #one-tailed
#
# KVtest(X.list, method = "GCT", L=10)
# pnorm(highD2pop::GCT.test(X.list[[1]], X.list[[2]], r=10, ntoorderminus=0)$TSvalue, lower.tail = F) #one-tailed
################################################################################

