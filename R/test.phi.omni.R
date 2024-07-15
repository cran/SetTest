########################################################################
#' calculate the right-tail probability of omnibus phi-divergence statistics under general correlation matrix.
#' @param prob - vector of input p-values.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param K0 - vector of search range starts (from the k0th smallest p-value).
#' @param K1 - vector of search range ends (at the k1th smallest p-value).
#' @param S - vector of the phi-divergence test parameters.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @param method - default = "ecc": the effective correlation coefficient method in reference 2. "ave": the average method in reference 3, which is an earlier version of reference 2. The "ecc" method is more accurate and numerically stable than "ave" method.
#' @param ei - the eigenvalues of M if available. 
#' @return p-value of the omnibus test. 
#' @return p-values of the individual phi-divergence test.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and power of optimal signal-detection statistics in finite case", IEEE Transactions on Signal Processing (2020) 68, 1021-1033
#' 2. Hong Zhang and Zheyang Wu. "The general goodness-of-fit tests for correlated data", Computational Statistics & Data Analysis (2022) 167, 107379
#' 3. Hong Zhang and Zheyang Wu. "Generalized Goodness-Of-Fit Tests for Correlated Data", arXiv:1806.03668.
#' @examples
#' M = matrix(0.3,10,10) + diag(1-0.3, 10)
#' test.phi.omni(runif(10), M=M, K0=rep(1,2), K1=rep(5,2), S=c(1,2))
#' @export

test.phi.omni <- function(prob, M, K0, K1, S, onesided=FALSE, method="ecc", ei = NULL){
  result = stat.phi.omni(p=prob, M=M, K0=K0, K1=K1, S=S, onesided=onesided, method=method, ei = ei)
  return(list(pvalue=pphi.omni(result$minp, M=M, S=S, K0=K0, K1=K1, onesided=onesided, method=method, ei = ei), 
              indiv_pval=result$pval))
}
