#' Quantile of phi-divergence statitic under the null hypothesis.
#' @param p -  a scalar left probability that defines the quantile.
#' @param n - dimension parameter, i.e. the number of input statitics to construct phi-divergence statistic.
#' @param s - phi-divergence parameter. s = 2 is the higher criticism statitic.s = 1 is the Berk and Jones statistic.
#' @param beta - search range parameter . Beta must be between 1/n and 1.
#' @return Quantile of phi-divergence statistics.
#' @seealso \code{\link{stat.phi}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#' @examples
#' ## The 0.05 critical value of HC statistic when n = 10:
#' qphi(p=.95, n=10, s=2, beta=0.5)
#' @export
qphi <- function(p, n, s, beta){
  lower = -1
  upper = 10
  maxnumIter = 500
  numIter = 1
  p_cal_lower = pphi(lower, n, s, beta)
  p_cal_upper = pphi(upper, n, s, beta)
  while((p_cal_lower-p)*(p_cal_upper-p)>0 && numIter < maxnumIter){
    if(p_cal_lower>p){
      lower = lower*1.1
      p_cal_lower = pphi(lower, n, s, beta)
    }
    if(p_cal_upper < p){
      upper = upper*1.1
      p_cal_upper = pphi(upper, n, s, beta)
    }
    numIter = numIter + 1
  }
  q = mean(c(lower,upper))
  p_cal = pphi(q, n, s, beta)
  error = (p_cal-p)/p
  numIter = 1
  while(abs(error) > 1e-5 && numIter < maxnumIter){
    if(error > 0){
      upper = q
      q = mean(c(lower,upper))
      p_cal = pphi(q, n, s, beta)
      error = (p_cal-p)/p
    }else{
      lower = q
      q = mean(c(lower,upper))
      p_cal = pphi(q, n, s, beta)
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
