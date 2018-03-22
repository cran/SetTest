#' Quantile of Berk-Jones statitic under the null hypothesis.
#' @param p - a scalar left probability that defines the quantile.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @return Quantile of BJ statistics.
#' @seealso \code{\link{stat.bj}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#'
#' 3. Berk, R.H. & Jones, D.H. Z. "Goodness-of-fit test statistics that dominate the Kolmogorov statistics". Wahrscheinlichkeitstheorie verw Gebiete (1979) 47: 47.
#' @examples
#' ## The 0.05 critical value of BJ statistic when n = 10:
#' qbj(p=.95, M=diag(10), k0=1, k1=5, onesided=FALSE)
#' @export
qbj <- function(p, M, k0, k1, onesided=FALSE){
  return(qphi(p, M, k0, k1, s=1, onesided=FALSE))
}
