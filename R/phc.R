#' CDF of Higher Criticism statistic under the null hypothesis.
#' @param q - quantile, must be a scalar.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @param method - default = "ecc": the effective correlation coefficient method in reference 2. "ave": the average method in reference 3, which is an earlier version of reference 2. The "ecc" method is more accurate and numerically stable than "ave" method.
#' @param ei - the eigenvalues of M if available. 
#' @seealso \code{\link{stat.hc}} for the definition of the statistic.
#' @return The left-tail probability of the null distribution of HC statistic at the given quantile.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and power of optimal signal-detection statistics in finite case", IEEE Transactions on Signal Processing (2020) 68, 1021-1033
#' 2. Hong Zhang and Zheyang Wu. "The general goodness-of-fit tests for correlated data", Computational Statistics & Data Analysis (2022) 167, 107379
#' 3. Hong Zhang and Zheyang Wu. "Generalized Goodness-Of-Fit Tests for Correlated Data", arXiv:1806.03668.
#' @examples
#' pval <- runif(10)
#' hcstat <- stat.phi(pval, s=2, k0=1, k1=5)$value
#' phc(q=hcstat, M=diag(10), k0=1, k1=10)
#' @export
phc <- function(q, M, k0, k1, onesided=FALSE, method="ecc", ei=NULL){
  return(pphi(q=q, M=M, k0=k0, k1=k1, s=2, onesided=onesided, method=method, ei=ei))
}
