#' Quantile of phi-divergence statitic under the null hypothesis.
#' @param p -  a scalar left probability that defines the quantile.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param s - the phi-divergence test parameter.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @return Quantile of phi-divergence statistics.
#' @seealso \code{\link{stat.phi}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#' @examples
#' ## The 0.05 critical value of HC statistic when n = 10:
#' qphi(p=.95, M=diag(10), k0=1, k1=5, s=2, onesided=FALSE)
#' @export
qphi <- function(p, M, k0, k1, s=2, onesided=FALSE){
  lower = -1
  upper = 10
  maxnumIter = 500
  numIter = 1
  p_cal_lower = suppressWarnings(pphi(lower, M, k0, k1, s, onesided))
  p_cal_upper = suppressWarnings(pphi(upper, M, k0, k1, s, onesided))
  while((p_cal_lower-p)*(p_cal_upper-p)>0 && numIter < maxnumIter){
    if(p_cal_lower>p){
      lower = lower*1.1
      p_cal_lower = suppressWarnings(pphi(lower, M, k0, k1, s, onesided))
    }
    if(p_cal_upper < p){
      upper = upper*1.1
      p_cal_upper = suppressWarnings(pphi(upper, M, k0, k1, s, onesided))
    }
    numIter = numIter + 1
  }
  q = mean(c(lower,upper))
  p_cal = suppressWarnings(pphi(q, M, k0, k1, s, onesided))
  error = (p_cal-p)/p
  numIter = 1
  while(abs(error) > 1e-5 && numIter < maxnumIter){
    if(error > 0){
      upper = q
      q = mean(c(lower,upper))
      p_cal = suppressWarnings(pphi(q, M, k0, k1, s, onesided))
      error = (p_cal-p)/p
    }else{
      lower = q
      q = mean(c(lower,upper))
      p_cal = suppressWarnings(pphi(q, M, k0, k1, s, onesided))
      error = (p_cal-p)/p
    }
    numIter = numIter + 1
  }
  if(numIter < maxnumIter){
    return(q)
  }else{
    return("Bisection fails, choose another lower or upper bound to try again")
  }
}
