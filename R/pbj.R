#' CDF of Berk-Jones statitic under the null hypothesis.
#' @param q - quantile, must be a scalar.
#' @param n - dimension parameter, i.e. the number of input statitics to construct BJ statistic.
#' @param beta - search range parameter. Beta must be between 1/n and 1.
#' @return The left-tail probability of the null distribution of B-J statistic at the given quantile.
#' @seealso \code{\link{stat.bj}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal.
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#'
#' 3. Berk, R.H. & Jones, D.H. Z. "Goodness-of-fit test statistics that dominate the Kolmogorov statistics". Wahrscheinlichkeitstheorie verw Gebiete (1979) 47: 47.
#' @examples
#' pval <- runif(10)
#' bjstat <- stat.phi(pval, 1, 0.5)$value
#' pbj(q=bjstat, n=10, beta=0.5)
#' @export
pbj <- function(q, n, beta){
  return(pphi(q, n, 1, beta))
}
