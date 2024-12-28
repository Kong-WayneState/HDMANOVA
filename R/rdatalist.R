#' @name rdatalist
#' @title Simulate multi-sample of high-dimensional data for MANOVA.
#' @description The data are generated from various time series dependence structures and innovation distributions:
#' \deqn{X_{ij}=\mu_{i}+\varepsilon_{ij},\quad 1\le i \le g,\ 1\le j\le n_i.}
#'
#'
#' @param g number of samples, \eqn{g \ge 2}.
#' @param p dimension.
#' @param nvec a vector of sample sizes, length \eqn{g}. If length 1, the same sample size will be used.
#' @param beta proportion of non-equal elements in mean \eqn{\mu_g}.
#' @param delta departure magnitudes in mean \eqn{\mu_g}.
#' @param dist innovation distribution. If length 1, the same distribution will be used. \describe{
#' \item{"normal"}{Normal(mean = 0, sd = 1);}
#' \item{"t"}{t(df = 3);}
#' \item{"cauchy"}{Cauchy(location = 0, scale = 0.1);}
#' \item{"gamma"}{centered Gamma(shape = 4, scale = 2).}
#' }
#' @param err dependence structures for error processes \eqn{\{\varepsilon_{ijk}: k = 1,\ldots, p\}}. If length 1, the same dependence structure will be used. \describe{
#'  \item{"IND"}{(Independent) independent generation from the given innovation distribution.}
#'  \item{"WD"}{(Weakly Dependent) random generation of ARMA(2, 2) model from the given innovation distribution, where autoregressive coefficients \eqn{\phi = (0.4, −0.1)} and moving-average coefficients \eqn{\theta = (0.2, 0.3)};}
#'  \item{"SD"}{(Strongly Dependent) random generation of AR(1) model from the given innovation distribution, where autoregressive coefficients \eqn{\phi = 0.9};}
#'  \item{"LD"}{(Long-range Dependent) \eqn{\varepsilon_{ij}=U^t\eta_{ij}}, where \eqn{\eta_{ij}} is independently generated \eqn{p}-dim vector from the given innovation distribution and \eqn{U^tU = A} is the Cholesky factorization of the \eqn{p \times p} matrix \eqn{A=(a_{st})}, \eqn{a_{st}=0.5[(|s-t|+1)^{1.4}+(|s-t|-1)^{1.4}-2|s-t|^{1.4}]}.}
#'}
#'
#'
#' @return a list of \eqn{n_i \times p} sample matrices for \eqn{i = 1,\ldots, g}.
#'
#' @export
#'
#'
#' @examples
#' set.seed(123)
#' X.list <- rdatalist()
#' str(X.list)
#' X.list <- rdatalist(g = 3,
#'                     p = 300,
#'                     nvec = c(80, 90, 100),
#'                     beta = 0.1,
#'                     delta = 0.5,
#'                     dist = c("t", "normal","normal"),
#'                     err = c("IND","WD","SD")
#'                     )
#' str(X.list)
#'

rdatalist <- function(g = 3, p = 100, nvec = 50, beta = 0, delta = 0,
                     dist = "normal", err = "IND"){
  mulist <- lapply(1:(g-1), function(i){rep(0, p)})
  mulist[[g]] <- c(rep(delta, round(p*beta)), rep(0, p-round(p*beta)))
  if(length(nvec) == 1) nvec = rep(nvec, g)
  if(length(dist) == 1) dist = rep(dist, g)
  if(length(err) == 1) err = rep(err, g)

  OneSamp <-function(n, dist, err, mu){
    if (err == "LD"){
      rgenfun <- eval(parse(text=paste0("rld", "_", dist)))
      X <- rgenfun(n=n, mu=mu, H=0.7)
    }else{
      if (err == "IND"){
        arma <- list(order = c(0, 0, 0))
      } else if (err == "WD"){
        arma <- list(ar = c(.4, -.1), ma = c(0.2, 0.3))
      } else if (err == "SD"){
        arma <- list(ar = 0.9)
      }
      rgenfun <- eval(parse(text=paste0("rarma", "_", dist)))
      X <- rgenfun(n=n, mu=mu, arma=arma)
    }
    return(X)
  }


  X.list = lapply(1:g, function(i){OneSamp(n=nvec[i], dist = dist[i], err = err[i], mu=mulist[[i]])})

  return(X.list)

}


#############################################################

#' @keywords internal
#' @noRd
# ARMA(2,2) with N(0, 1) Innovations
rarma_normal <- function(n, mu, arma, center=0, sig=1){
  arma_normal<- function(nsim, p, center, sig, arma){
    e = stats::arima.sim(n = p, model = arma,
                  rand.gen = function(n, ...) stats::rnorm(n, mean=center, sd = sig))
    return(e)
  }
  e = sapply(rep(1,n), arma_normal, p=length(mu), center=center, sig=sig, arma=arma)
  X = t(mu + e)
  return(X)
}

#' @keywords internal
#' @noRd
# ARMA(2,2) with t(3) Innovations
rarma_t <- function(n, mu, arma, q=3){
  arma_t<- function(nsim, p, q, arma){
    e = stats::arima.sim(n = p, model = arma,
                  rand.gen = function(n, ...) stats::rt(n, df = q))
    return(e)
  }
  e = sapply(rep(1,n), arma_t, p=length(mu), q=q, arma=arma)
  X = t(mu + e)
  return(X)
}

#' @keywords internal
#' @noRd
# ARMA(2,2) with Cauchy(0, 0.1) Innovations
rarma_cauchy <- function(n, mu, arma, center=0, sig=0.1){
  arma_cauchy <- function(nsim, p, center, sig, arma){
    e = stats::arima.sim(n = p, model = arma,
                  rand.gen = function(n, ...) stats::rcauchy(n, location=center, scale=sig))
    return(e)
  }

  e = sapply(rep(1,n), arma_cauchy, p=length(mu), center=center, sig=sig, arma=arma)
  X = t(mu + e)
  return(X)
}

#' @keywords internal
#' @noRd
# ARMA(2,2) with centered Gamma(4, 2) Innovations
rarma_gamma <- function(n, mu, arma, a=4, s=2){
  arma_gamma <- function(nsim, p, a, s, arma){
    e = stats::arima.sim(n = p, model = arma,
                  rand.gen = function(n, ...) (stats::rgamma(n, shape=a, scale = s) - a*s))
    return(e)
  }

  e = sapply(rep(1,n), arma_gamma, p=length(mu), a=a, s=s, arma=arma)
  X = t(mu + e)
  return(X)
}

#' @keywords internal
#' @noRd
# Long-range dependent N(0, 1)
rld_normal <- function(n, mu, H = .7, center=0, sig=1){
  p = length(mu)
  R <- matrix(0, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      k <- abs(i - j)
      R[i, j] = R[j, i] = .5 * ((k + 1) ^ (2 * H) + (k - 1) ^ (2 * H) - 2 * k ^ (2 * H))
    }
  }
  diag(R) <- 1
  U <- chol(R)
  rand.gen = stats::rnorm(n * p, mean = center, sd = sig)
  e <- matrix(rand.gen, nrow = n) %*% U
  X = t(mu + t(e))
  return(X)
}

#' @keywords internal
#' @noRd
# Long-range dependent t(3)
rld_t <- function(n, mu, H = .7, q = 3){
  p = length(mu)
  R <- matrix(0, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      k <- abs(i - j)
      R[i, j] = R[j, i] = .5 * ((k + 1) ^ (2 * H) + (k - 1) ^ (2 * H) - 2 * k ^ (2 * H))
    }
  }
  diag(R) <- 1
  U <- chol(R)
  rand.gen = stats::rt(n * p, df = q)
  e <- matrix(rand.gen, nrow = n) %*% U
  X = t(mu + t(e))
  return(X)
}

#' @keywords internal
#' @noRd
# Long-range dependent Cauchy(0, 0.1)
rld_cauchy <- function(n, mu, H = .7, center=0, sig=0.1){
  p = length(mu)
  R <- matrix(0, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      k <- abs(i - j)
      R[i, j] = R[j, i] = .5 * ((k + 1) ^ (2 * H) + (k - 1) ^ (2 * H) - 2 * k ^ (2 * H))
    }
  }
  diag(R) <- 1
  U <- chol(R)
  rand.gen = stats::rcauchy(n * p, location=center, scale=sig)
  e <- matrix(rand.gen, nrow = n) %*% U
  X = t(mu + t(e))
  return(X)
}

#' @keywords internal
#' @noRd
# Long-range dependent centered Gamma(4, 2)
rld_gamma <- function(n, mu, H = .7, a=4, s=2){
  p = length(mu)
  R <- matrix(0, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      k <- abs(i - j)
      R[i, j] = R[j, i] = .5 * ((k + 1) ^ (2 * H) + (k - 1) ^ (2 * H) - 2 * k ^ (2 * H))
    }
  }
  diag(R) <- 1
  U <- chol(R)
  rand.gen = stats::rgamma(n * p, shape=a, scale = s) - a*s
  e <- matrix(rand.gen, nrow = n) %*% U
  X = t(mu + t(e))
  return(X)
}



