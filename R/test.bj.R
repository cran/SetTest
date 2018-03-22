#' Multiple comparison test using Berk and Jones (BJ) statitics.
#' @param prob - vector of input p-values.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @return pvalue - the p-value of the Berk-Jones test.
#' @return bjstat - the Berk-Jones statistic.
#' @return location - the order of the input p-values to obtain BJ statistic.
#' @seealso \code{\link{stat.bj}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Jager, Leah; Wellner, Jon A. "Goodness-of-fit tests via phi-divergences". Annals of Statistics 35 (2007).
#'
#' 3. Berk, R.H. & Jones, D.H. Z. "Goodness-of-fit test statistics that dominate the Kolmogorov statistics". Wahrscheinlichkeitstheorie verw Gebiete (1979) 47: 47.
#' @examples
#' test.bj(runif(10), M=diag(10), k0=1, k1=10)
#' #When the input are statistics#
#' stat.test = rnorm(20)
#' p.test = 2*(1 - pnorm(abs(stat.test)))
#' test.bj(p.test, M=diag(20), k0=1, k1=10)
#' @export
test.bj <- function(prob, M, k0, k1, onesided=FALSE){
  n = length(prob)
  result = stat.phi(prob, s=1, k0=k0, k1=k1)
  return(list(pvalue=1-pphi(result$value, M=M, s=1, k0=k0, k1=k1, onesided=onesided), bjstat=result$value, location=which(result$value==max(result$value))))
}
