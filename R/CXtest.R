#' @name CXtest
#' @title High dimensional MANOVA
#' @description Testing the equality of multi-sample high dimensional mean vectors using the testing procedure by Cai and Xia (2014). For the two-sample test, see also Cai, Liu and Xia (2014).
#'
#' @param X.list list of \eqn{g (\ge 2)} matrices from the samples.
#' @param A transformation matrix. \describe{
#' \item{"Identity"}{the \eqn{p\times p} identity matrix is used;}
#' \item{"Omega"}{the CLIME estimator (Cai, Liu and Luo, 2011) for precision matrix is used when the precision matrix is known to be sparse;}
#' \item{"Sigma"}{the adaptive thresholding estimator (Cai and Liu, 2011) for covariance  matrix is used when the precision matrix is not known to be sparse.}
#' }
#' @param pkg R package ("flare" or "fastclime") used for CLIME estimator when \code{A = "Omega"}.
#' @param biased a logical value indicating whether to use biased estimator for the pooled (co)variance.
#'
#' @return value of test statistic, p-value.
#'
#' @references Cai, T., Liu, W., and Luo, X. (2011). A constrained \eqn{l_1} minimization approach to sparse precision matrix estimation. \emph{Journal of the American Statistical Association}, 106:594–607.
#' @references Cai, T. and Liu, W. (2011). Adaptive thresholding for sparse covariance matrix estimation. \emph{Journal of the American Statistical Association}, 106:672–684.
#' @references Cai, T., and Xia, Y. (2014). High-dimensional sparse MANOVA.  \emph{Journal of Multivariate Analysis}, 131:174–196.
#' @references Cai, T., Liu, W., and Xia, Y. (2014). Two-sample test of high dimensional means under dependence.  \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology}, 76(2):349–372.
#' @export
#'
#'
#' @examples
#' # generating a list of 3-sample data
#' set.seed(123)
#' X.list <- rdatalist(g = 3)
#' str(X.list)
#' CXtest(X.list)
#' CXtest(X.list, A = "Omega")
#'


CXtest <- function(X.list, A = "Identity", pkg = "flare", biased=FALSE){
  # delta is used only if sparse = "covariance"
  g <- length(X.list)
  pvec <- unlist(lapply(X.list, ncol))
  if (all(pvec != mean(pvec))){stop("dimensions are not equal!")}
  p <- pvec[1]
  nvec <- unlist(lapply(X.list, nrow))
  n <- sum(nvec)

  X.cov <- lapply(X.list, stats::cov)
  X.ncov <- Map("*", X.cov, (nvec-1))
  nS <- Reduce("+", X.ncov)
  if(biased){S.pooled <- nS/n}else{S.pooled <- nS/(n-g)}

  if (A == "Identity"){
    Omega.est = diag(p)
  } else if (A == "Omega"){ # estimating Omega
    if (pkg == "fastclime"){
      # mean2.2014CLX in SHT use unbiased, CLX paper used biased
      # CLIME: fastclime (not available)
      # sink(tempfile())
      # Omega.list <- fastclime::fastclime(S.pooled)
      # Omega.select <- fastclime::fastclime.selector(Omega.list$lambdamtx, Omega.list$icovlist, lambda=(log(p)/n))
      # sink()
      # Omega.est <- Omega.select$icov
    }else if (pkg == "flare"){
      # CLIME: flare (very slow)
      X.ctr <- lapply(X.list, scale, center=TRUE, scale = FALSE)
      alldata <- do.call(rbind, X.ctr)
      sink(tempfile())
      Omega.list <- flare::sugm(alldata, method="clime")
      Omega.select <- flare::sugm.select(Omega.list, criterion = "cv")
      sink()
      Omega.est <- Omega.select$opt.icov
    }
  } else if (A == "Sigma") { # estimating Sigma
    X.ctr <- lapply(X.list, scale, center=TRUE, scale = FALSE)
    lambda.mat.list <- lapply(X.ctr, function(X){((t(X^2)%*%(X^2))+nrow(X)*S.pooled^2-2*(t(X)%*%X)*S.pooled)/n})
    lambda.mat <- 2*sqrt(Reduce("+",lambda.mat.list)*log(p)/n)
    Sigma.hat <- S.pooled * (abs(S.pooled) >= lambda.mat)

    Sigma.eig <- base::eigen(Sigma.hat)
    s <- Sigma.eig$values
    s.new <- sapply(s, function(x){ifelse(x<=0, log(p)/n, x)}) #adjust for non-definite case
    v = Sigma.eig$vectors
    Sigma.est = v %*% diag(s.new) %*% t(v)
    Omega.est = solve(Sigma.est)
  }

  X.omega <- lapply(X.list, function(X){X%*%Omega.est})
  X.o.mean <- t(mapply(colMeans, X.omega))
  X.o.var <- t(mapply(function(X){apply(X, 2, stats::var)}, X.omega))
  if(biased){# CLX used biased, CX used unbiased
    s.pooled <- colSums(X.o.var*(nvec-1))/n
  }else{
    s.pooled <- colSums(X.o.var*(nvec-1))/(n-g)
  }

  tsq.pw <- function(xvec, nvec){# pairwise squared t-statistics
    g <- length(xvec)
    pairwise <- upper.tri(matrix(0, nrow=g, ncol=g))

    xdiff <- outer(xvec, xvec, "-")[pairwise]
    nprod <- outer(nvec, nvec, "*")[pairwise]
    nsum <- outer(nvec, nvec, "+")[pairwise]

    tsq <- nprod/nsum*(xdiff^2)
    return(tsq)
  }

  tsq.pw.mat <- t(t(apply(X.o.mean, 2, tsq.pw, nvec=nvec))/s.pooled)# b*p

  if (g==2){
    Mstat <- max(tsq.pw.mat)  # b=1, Sigma0=1, sig.sq=1, d=1, H=1
    pval <- 1 - exp(-exp(-(Mstat-2*log(p)+log(log(p)))/2)/sqrt(pi))
  }else{
    Mstat <- max(colSums(tsq.pw.mat))
    b <- g*(g-1)/2
    Sigma0 <- matrix(0, b, b)
    row <- 0
    for (i in 1:(g-1)){
      for (j in (i+1):g){
        row = row+1
        col = row
        for (k in i:(g-1)){
          if(k == i){
            for(l in j:g){
              if (l==j) {
                Sigma0[row, col] = 1/2
              }else{
                Sigma0[row, col] = sqrt((nvec[j]*nvec[l])/((nvec[i]+nvec[j])*(nvec[k]+nvec[l])))
              }
              col = col + 1
            }
          }else{
            for (l in (k+1):g){
              if (j == k){
                Sigma0[row, col] = -sqrt((nvec[i]*nvec[l])/((nvec[i]+nvec[j])*(nvec[k]+nvec[l])))
              } else if (j == l){
                Sigma0[row, col] = sqrt((nvec[i]*nvec[k])/((nvec[i]+nvec[j])*(nvec[k]+nvec[l])))
              } else{
                Sigma0[row, col] = 0
              }
              col = col + 1
            }
          }
        }
      }
    }
    Sigma0 <- Sigma0 + t(Sigma0)
    evalue <- eigen(Sigma0)$values
    sig.sq <- evalue[1]
    d <- sum(evalue == sig.sq)
    H <- 1/prod(sqrt(1-evalue[-(1:d)]/sig.sq))
    pval <- 1 - exp(-H/gamma(d/2)*exp(-(Mstat/sig.sq-2*log(p)-(d-2)*log(log(p)))/2))
  }

  return(list(statistic=Mstat, p.value=pval))

}




################################################################################
## Note: when g = 2, the follwing test gives the same result:
# set.seed(123)
# X.list <- rdatalist(g = 2, err = "WD")
# CXtest(X.list, A = "Identity")
# highmean::apval_Cai2014(X.list[[1]], X.list[[2]], eq.cov = TRUE)
#
# # CLIME estimator for precision matrix
# set.seed(12)
# CXtest(X.list, A = "Omega", pkg = "flare", biased = TRUE)
# set.seed(12)
# SHT::mean2.2014CLX(X.list[[1]], X.list[[2]],  precision="sparse", cov.equal = TRUE)
#
# # adaptive thresholding estimator for covariance matrix
# CXtest(X.list, A = "Sigma", biased = TRUE)
# HDtest::CLX(t(X.list[[1]]), t(X.list[[2]]), alpha=0.05, DNAME = "A")
# SHT::mean2.2014CLX(X.list[[1]], X.list[[2]], precision="unknown", cov.equal = TRUE)  #the adjusted eigenvalue is too big.
################################################################################

