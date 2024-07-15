#' Multiple comparison test using Berk and Jones (BJ) statitics.
#' @param prob - vector of input p-values.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @param method - default = "ecc": the effective correlation coefficient method in reference 2. "ave": the average method in reference 3, which is an earlier version of reference 2. The "ecc" method is more accurate and numerically stable than "ave" method.
#' @param ei - the eigenvalues of M if available. 
#' @return pvalue - the p-value of the Berk-Jones test.
#' @return bjstat - the Berk-Jones statistic.
#' @return location - the order of the input p-values to obtain BJ statistic.
#' @seealso \code{\link{stat.bj}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and power of optimal signal-detection statistics in finite case", IEEE Transactions on Signal Processing (2020) 68, 1021-1033
#' 2. Hong Zhang and Zheyang Wu. "The general goodness-of-fit tests for correlated data", Computational Statistics & Data Analysis (2022) 167, 107379
#' 3. Hong Zhang and Zheyang Wu. "Generalized Goodness-Of-Fit Tests for Correlated Data", arXiv:1806.03668.
#' 4. Leah Jager and Jon Wellner. "Goodness-of-fit tests via phi-divergences". Annals of Statistics 35 (2007).
#' @examples
#' test.bj(runif(10), M=diag(10), k0=1, k1=10)
#' #When the input are statistics#
#' stat.test = rnorm(20)
#' p.test = 2*(1 - pnorm(abs(stat.test)))
#' test.bj(p.test, M=diag(20), k0=1, k1=10)
#' @export
test.bj <- function(prob, M, k0, k1, onesided=FALSE, method="ecc", ei = NULL){
  n = length(prob)
  result = stat.phi(prob, s=1, k0=k0, k1=k1)
  return(list(pvalue=1-pphi(result$value, M=M, s=1, k0=k0, k1=k1, onesided=onesided, method=method, ei = ei), 
              bjstat=result$value, 
              location=which(result$value==max(result$value))))
}
