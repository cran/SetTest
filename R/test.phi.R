########################################################################
#' Multiple comparison test using phi-divergence statistics.
#' @param prob - vector of input p-values.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param s - phi-divergence parameter. s = 2 is the higher criticism statitic.s = 1 is the Berk and Jones statistic.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @param method - default = "ecc": the effective correlation coefficient method in reference 2. "ave": the average method in reference 3, which is an earlier version of reference 2. The "ecc" method is more accurate and numerically stable than "ave" method.
#' @param ei - the eigenvalues of M if available. 
#' @return pvalue - The p-value of the phi-divergence test.
#' @return phistat - phi-diergence statistic.
#' @return location - the order of the input p-values to obtain phi-divergence statistic.
#' @seealso \code{\link{stat.phi}} for the definition of the statistic.v
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and power of optimal signal-detection statistics in finite case", IEEE Transactions on Signal Processing (2020) 68, 1021-1033
#' 2. Hong Zhang and Zheyang Wu. "The general goodness-of-fit tests for correlated data", Computational Statistics & Data Analysis (2022) 167, 107379
#' 3. Hong Zhang and Zheyang Wu. "Generalized Goodness-Of-Fit Tests for Correlated Data", arXiv:1806.03668.
#' 4. Leah Jager and Jon Wellner. "Goodness-of-fit tests via phi-divergences". Annals of Statistics 35 (2007).
#' @examples
#' stat.test = rnorm(20) # Z-scores
#' p.test = 2*(1 - pnorm(abs(stat.test)))
#' test.phi(p.test, M=diag(20), s = 0.5, k0=1, k1=10)
#' test.phi(p.test, M=diag(20), s = 1, k0=1, k1=10)
#' test.phi(p.test, M=diag(20), s = 2, k0=1, k1=10)
#' @export
test.phi <- function(prob, M, k0, k1, s=2, onesided=FALSE, method="ecc", ei = NULL){
  n = length(prob)
  result = stat.phi(prob, s, k0, k1)
  return(list(pvalue=1-pphi(result$value, M=M, s=s, k0=k0, k1=k1, onesided=onesided, method=method, ei = ei), 
              phistat=result$value, 
              location=which(result$value==max(result$value))))
}
