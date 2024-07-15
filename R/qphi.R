#' Quantile of phi-divergence statistic under the null hypothesis.
#' @param p -  a scalar left probability that defines the quantile.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param s - the phi-divergence test parameter.
#' @param t - numerical truncation parameter.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @param method - default = "ecc": the effective correlation coefficient method in reference 2. "ave": the average method in reference 3, which is an earlier version of reference 2. The "ecc" method is more accurate and numerically stable than "ave" method.
#' @param ei - the eigenvalues of M if available. 
#' @param err_thr - the error threshold. The default value is 1e-4. 
#' @return Quantile of the phi-divergence statistics.
#' @seealso \code{\link{stat.phi}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and power of optimal signal-detection statistics in finite case", IEEE Transactions on Signal Processing (2020) 68, 1021-1033
#' 2. Hong Zhang and Zheyang Wu. "The general goodness-of-fit tests for correlated data", Computational Statistics & Data Analysis (2022) 167, 107379
#' 3. Hong Zhang and Zheyang Wu. "Generalized Goodness-Of-Fit Tests for Correlated Data", arXiv:1806.03668.
#' @examples
#' qphi(p=.95, M=diag(10), k0=1, k1=5, s=2, onesided=FALSE)
#' qphi(p=1-1e-3, M=diag(10), k0=1, k1=5, s=2, onesided=FALSE)
#' qphi(p=1-1e-3, M=diag(10), k0=1, k1=5, s=2, onesided=FALSE, err_thr = 1e-6)
#' qphi(p=1-1e-5, M=diag(10), k0=1, k1=5, s=2, onesided=FALSE)
#' qphi(p=1-1e-5, M=diag(10), k0=1, k1=5, s=2, onesided=FALSE, err_thr = 1e-6)
#' qphi(p=1-1e-5, M=diag(10), k0=1, k1=5, s=2, onesided=FALSE, err_thr = 1e-8)
#' @export
qphi <- function(p, M, k0, k1, s=2, t=30, onesided=FALSE, method="ecc", ei=NULL, err_thr=1e-4){
  lower = -1
  upper = 10
  maxnumIter = 500
  numIter = 1
  p_cal_lower = suppressWarnings(pphi(lower, M, k0, k1, s, t, onesided, method, ei))
  p_cal_upper = suppressWarnings(pphi(upper, M, k0, k1, s, t, onesided, method, ei))
  while((p_cal_lower-p)*(p_cal_upper-p)>0 && numIter < maxnumIter){
    if(p_cal_lower>p){
      lower = lower*1.1
      p_cal_lower = suppressWarnings(pphi(lower, M, k0, k1, s, t, onesided, method, ei))
    }
    if(p_cal_upper < p){
      upper = upper*1.1
      p_cal_upper = suppressWarnings(pphi(upper, M, k0, k1, s, t, onesided, method, ei))
    }
    numIter = numIter + 1
  }
  q = mean(c(lower,upper))
  p_cal = suppressWarnings(pphi(q, M, k0, k1, s, t, onesided, method, ei))
  error = (p_cal-p)/p
  numIter = 1
  while(abs(error) > err_thr && upper-lower>1e-4 && numIter < maxnumIter){
    if(error > 0){
      upper = q
      q = mean(c(lower,upper))
      p_cal = suppressWarnings(pphi(q, M, k0, k1, s, t, onesided, method, ei))
      error = (p_cal-p)/p
    }else{
      lower = q
      q = mean(c(lower,upper))
      p_cal = suppressWarnings(pphi(q, M, k0, k1, s, t, onesided, method, ei))
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
