#' CDF of Berk-Jones statitic under the null hypothesis.
#' @param q - quantile, must be a scalar.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @return The left-tail probability of the null distribution of B-J statistic at the given quantile.
#' @seealso \code{\link{stat.bj}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#'
#' 3. Berk, R.H. & Jones, D.H. Z. "Goodness-of-fit test statistics that dominate the Kolmogorov statistics". Wahrscheinlichkeitstheorie verw Gebiete (1979) 47: 47.
#' @examples
#' pval <- runif(10)
#' bjstat <- stat.phi(pval, s=1, k0=1, k1=10)$value
#' pbj(q=bjstat, M=diag(10), k0=1, k1=10)
#' @export
pbj <- function(q, M, k0, k1, onesided=FALSE){
  return(pphi(q=q, M=M, k0=k0, k1=k1, s=1, onesided=onesided))
}
