#' CDF of phi-divergence statitic under the null hypothesis.
#' @param q - quantile, must be a scalar.
#' @param n - dimension parameter, i.e. the number of input statitics to construct phi-divergence statistic.
#' @param s - phi-divergence parameter. s = 2 is the higher criticism statitic.s = 1 is the Berk and Jones statistic.
#' @param beta - search range parameter . Beta must be between 1/n and 1.
#' @return The left-tail probability of the null distribution of  phi-divergence statistic at the given quantile.
#' @seealso \code{\link{stat.phi}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Jager, Leah; Wellner, Jon A. "Goodness-of-fit tests via phi-divergences". Annals of Statistics 35 (2007).
#' @examples
#' pval <- runif(10)
#' hcstat <- stat.phi(pval, 2, 0.5)$value
#' pphi(q=hcstat, n=10, s=2, beta=0.5)
#' @export
#' @importFrom stats dpois pbeta
pphi <- function(q, n, s, beta){
  t = 28
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

  #####
  um = boundary_end
  u = boundary
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
      return(1-drop(p))
    }else{
      #warning("q is too small, p is only an upper bound", call. = FALSE)
      #return(1-pbeta(u[2], 2, n-1)+u[1]*n*(1-pbeta(u[2],1,n-1)))
      warning("q is too small. The algorithm fails due to loss of significant digits", call. = FALSE)
      return(0)
    }
  }else{
    p = pbeta(u[1], 1, n)
    return(1-drop(p))
  }
}
