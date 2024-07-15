#' Quantile of Berk-Jones statistic under the null hypothesis.
#' @param p - a scalar left probability that defines the quantile.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @param method - default = "ecc": the effective correlation coefficient method in reference 2. "ave": the average method in reference 3, which is an earlier version of reference 2. The "ecc" method is more accurate and numerically stable than "ave" method.
#' @param ei - the eigenvalues of M if available. 
#' @param err_thr - the error threshold. The default value is 1e-4. 
#' @return Quantile of BJ statistics.
#' @seealso \code{\link{stat.bj}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and power of optimal signal-detection statistics in finite case", IEEE Transactions on Signal Processing (2020) 68, 1021-1033
#' 2. Hong Zhang and Zheyang Wu. "The general goodness-of-fit tests for correlated data", Computational Statistics & Data Analysis (2022) 167, 107379
#' 3. Hong Zhang and Zheyang Wu. "Generalized Goodness-Of-Fit Tests for Correlated Data", arXiv:1806.03668.
#' @examples
#' ## The 0.05 critical value of BJ statistic when n = 10:
#' qbj(p=.95, M=diag(10), k0=1, k1=5, onesided=FALSE)
#' qbj(p=1-1e-5, M=diag(10), k0=1, k1=5, onesided=FALSE, err_thr=1e-8)
#' @export
qbj <- function(p, M, k0, k1, onesided=FALSE, method="ecc", ei=NULL, err_thr=1e-4){
  return(qphi(p, M, k0, k1, s=1, onesided=onesided, method=method, ei=ei, err_thr=err_thr))
}
