#' Statistical power of phi-divergence test.
#' @param alpha - type-I error rate.
#' @param n - dimension parameter, i.e. the number of input statitics to construct phi-divergence statistic.
#' @param s - phi-divergence parameter. s = 2 is the higher criticism statitic.s = 1 is the Berk and Jones statistic.
#' @param beta - search range parameter. Search range = (1, beta*n). Beta must be between 1/n and 1.
#' @param method - different alternative hypothesis, including mixtures such as, "gaussian-gaussian", "gaussian-t", "t-t", "chisq-chisq", and "exp-chisq". By default, we use Gaussian mixture.
#' @param eps - mixing parameter of the mixture.
#' @param mu - mean of non standard Gaussian model.
#' @param df - degree of freedom of t/Chisq distribution and exp distribution.
#' @param delta - non-cenrality of t/Chisq distribution.
#' @return Power of phi-divergence test.
#' @details We consider the following hypothesis test,
#' \deqn{H_0: X_i\sim F, H_a: X_i\sim G}
#' Specifically, \eqn{F = F_0} and \eqn{G = (1-\epsilon)F_0+\epsilon F_1}, where \eqn{\epsilon} is the mixing parameter, \eqn{F_0} and \eqn{F_1} is
#' speified by the "method" argument:
#'
#' "gaussian-gaussian": \eqn{F_0} is the standard normal CDF and \eqn{F = F_1} is the CDF of normal distribution with \eqn{\mu} defined by mu and \eqn{\sigma = 1}.
#'
#' "gaussian-t": \eqn{F_0} is the standard normal CDF and \eqn{F = F_1} is the CDF of t distribution with degree of freedom defined by df.
#'
#' "t-t": \eqn{F_0} is the CDF of t distribution with degree of freedom defined by df and \eqn{F = F_1} is the CDF of non-central t distribution with degree of freedom defined by df and non-centrality defined by delta.
#'
#' "chisq-chisq": \eqn{F_0} is the CDF of Chisquare distribution with degree of freedom defined by df and \eqn{F = F_1} is the CDF of non-central Chisquare distribution with degree of freedom defined by df and non-centrality defined by delta.
#'
#' "exp-chisq": \eqn{F_0} is the CDF of exponential distribution with parameter defined by df and \eqn{F = F_1} is the CDF of non-central Chisqaure distribution with degree of freedom defined by df and non-centrality defined by delta.
#' @seealso \code{\link{stat.phi}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#' @examples
#' #If the alternative hypothesis Gaussian mixture with eps = 0.1 and mu = 1.2:#
#' power.phi(0.05, n=10, s=2, beta=0.5, eps = 0.1, mu = 1.2)
#' @export
#' @importFrom stats pchisq pnorm pt qchisq qexp qnorm qt
power.phi <- function(alpha, n, s, beta, method = "gaussian-gaussian",
                      eps=0, mu=0, df=1, delta=0){

  t = 28
  q = suppressWarnings(qphi(p=1-alpha, diag(10), s, k0=1, k1=beta*n))
  m = floor(beta*n)
  if(m<1){
    stop("Search Range is too small: beta*n < 1")
  }
  if(m>n){
    stop("Search Range is too large: beta*n > n")
  }
  t = min(t,m)
  M = seq(m-1,m-t,-1)
  ep = 10^(-15) # root-finding left endpoint
  boundary = rep(0,t)


  for (i in 1:t){
    boundary[i] = phi.f.inv.i(q, i, n, s)
  }
  boundary_end = phi.f.inv.i(q, m, n, s)

  if(method=="gaussian-gaussian"){
    g <- function(x)(1-eps)*x + eps*(1-pnorm(qnorm(1-x)-mu)) #Gaussian Mixture Alternative
  }
  if(method=="gaussian-t"){
    g <- function(x)(1-eps)*x + eps*(1-pt(qnorm(1-x),df)) #t alternative
  }
  if(method=="t-t"){
    g <- function(x)(1-eps)*x + eps*(1-pt(qt(1-x,df,0),df,delta)) #t mixture
  }
  if(method=="chisq-chisq"){
    g <- function(x)(1-eps)*x + eps*(1-pchisq(qchisq(1-x,df,0),df,delta)) #Chisq Mixture
  }
  if(method=="exp-chisq"){
    g <- function(x)(1-eps)*x + eps*(1-pchisq(qexp(1-x,df),df,delta)) #Chisq-EXP Mixture
  }
  # if(method=="gGaussian-gGaussian"){
  #   g <- function(x)(1-eps)*x + eps*(1-pnormp(qnormp(1-x,0,sigmap,p),mu,sigmap,p)) #Generalized Normal Mixture
  # }
  #####
  um = g(boundary_end)
  u = g(boundary)
  pp = exp(lfactorial(n) - lfactorial(n-(0:(t-1))))*pbeta(um,M+1,n-m+1,lower.tail=F)
  a = rep(1,t)
  if(t>1){
    for (i in 1:(t-1)){
      d=rep(0,i)
      for (k in 1:i){
        d[k]=dpois(k,u[i])*exp(u[i])
      }
      a[i+1]=-a[i:1]%*%d
    }
    p = 1-pp[1:t]%*%a[1:t]
    if(p>=pbeta(u[2], 2, n-1)-u[1]*n*(1-pbeta(u[2],1,n-1))&&p<=1){
      return(drop(p))
    }else{
      #warning("q is too small, p is only an upper bound", call. = FALSE)
      #return(1-pbeta(u[2], 2, n-1)+u[1]*n*(1-pbeta(u[2],1,n-1)))
      warning("q is too small. The algorithm fails due to loss of significant digits", call. = FALSE)
      return(1)
    }
  }else{
    p = pbeta(u[1], 1, n)
    return(drop(p))
  }
}
