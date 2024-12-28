#' @name HBtest
#' @title High dimensional MANOVA
#' @description Testing the equality of multi-sample high dimensional mean vectors using the testing procedure by Hu et. al. (2017) and Kong and Harrar (2021). For the two-sample test, see also Chen and Qin (2010).
#'
#' @param X.list list of \eqn{g (\ge 2)} matrices from the samples.
#' @param est method to be used for estimating tr(\eqn{\Sigma_i^2}). \describe{
#' \item{"BS"}{Bai and Saranadata, (1996);}
#' \item{"LC"}{Li and Chen, (2012);}
#' \item{"CQ"}{Chen and Qin, (2010).}
#' }
#'
#' @return value of test statistic and p-value.
#
#' @references Bai, Z. and Saranadasa, H. (1996). Effect of high dimension: By an example of a two sample problem. \emph{Statistica Sinica}, 6:311–329.
#' @references Li, J., Chen, S. (2012). Two sample tests for high-dimensional covariance matrices.  \emph{The Annals of Statistics}, 40(2):908–940.
#' @references Chen, S. and Qin, Y. (2010). A two-sample test for high-dimensional data with applications to gene-set testing.  \emph{The Annals of Statistics}, 38(2):808–835.
#' @references Hu, J., Bai, Z., Wang, C., and Wang, W. (2017). On testing the equality of high dimensional mean vectors with unequal covariance matrices.  \emph{Annals of the Institute of Statistical Mathematics}, 69(2):365–387.
#' @references Kong, X. and Harrar, S. (2021), High-dimensional MANOVA under weak conditions. \emph{Statistics: A Journal of Theoretical and Applied Statistics} 55:321–349.
#' @export
#'
#'
#' @examples
#' # generating a list of 3-sample data
#' set.seed(123)
#' X.list <- rdatalist(g = 3)
#' str(X.list)
#' HBtest(X.list)
#'


HBtest <- function(X.list, est = "LC"){
  g <- length(X.list)
  pvec <- unlist(lapply(X.list, ncol))
  if (all(pvec != mean(pvec))) {stop("dimensions are not equal!")}
  p <- pvec[1]
  nvec <- unlist(lapply(X.list, nrow))

  X.mean <- t(mapply(colMeans, X.list))
  X.var <- t(mapply(function(X){apply(X, 2, stats::var)}, X.list))

  HB.stat <- sum(colSums((scale(X.mean, center=TRUE, scale = FALSE))^2))*g -
    sum(colSums(X.var/nvec))*(g-1)

  X.cov <-  lapply(X.list, stats::cov)
  Sn.list <- Map("/", X.cov, nvec)
  Sn <- Reduce("+", Sn.list)
  Sn2.list <- Map("*", Sn.list, Sn.list)
  S2n <- Reduce("+", Sn2.list)

  if (est == "BS"){
    trSig1 <- function(S, n){
      (n-1)/(n*(n+1)*(n-2))*(sum(S^2) - sum(diag(S))^2/(n-1))
    }
    sigsq <- 2*(g-1)^2*sum(mapply(trSig1, X.cov, nvec)) + 2*sum(Sn^2-S2n)
  }

  if(est == "LC"){
    trSig2 <- function(X, n){
      A <- X%*%t(X)
      diag(A) <- 0
      c <- (sum(A^2)*(n-1)*(n-2) - 2*sum((rowSums(A))^2)*(n-1) +
             (sum(A))^2)/(n^2*(n-1)^2*(n-2)*(n-3))
      return(c)
    }
    sigsq <- 2*(g-1)^2*sum(mapply(trSig2, X.list, nvec)) + 2*sum(Sn^2-S2n)
  }
  if(est == "CQ") {
  trSig3 <- function(X, n){
    A <- X%*%t(X)
    diag(A) <- 0
    c <- (sum(A^2)*(n-1)^2-sum((rowSums(A))^2)*(2*n-1)+
           (sum(A))^2)/(n^2*(n-1)^2*(n-2)^2)
    return(c)
  }
  sigsq <- 2*(g-1)^2*sum(mapply(trSig3, X.list, nvec)) + 2*sum(Sn^2-S2n)
  }

  Zstat <- HB.stat/sqrt(sigsq)
  pval <- stats::pnorm(Zstat, lower.tail = FALSE)

  return(list(Statistic = Zstat, p.value = pval))
}




################################################################################
## Note: when g = 2, the follwing test gives the same result:
# set.seed(123)
# X.list <- rdatalist(g = 2, p = 500, nvec = c(200,250))
# HBtest(X.list, est = "CQ")
# highmean::apval_Chen2010(X.list[[1]], X.list[[2]], eq.cov=F)
# pnorm(highD2pop::ChenQin.test(X.list[[1]], X.list[[2]])$ChQ, lower.tail = F) # one-tailed
# pnorm(HDtest::CQ2(t(X.list[[1]]), t(X.list[[2]]), DNAME = "A")$statistics, lower.tail = F) #one-tailed
## As our code uses the simplified form, it is much faster when n or p is large.
################################################################################
