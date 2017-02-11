#' Quantile of Berk-Jones statitic under the null hypothesis.
#' @param p - a scalar left probability that defines the quantile.
#' @param n - dimension parameter, i.e. the number of input statitics to construct HC statistic.
#' @param beta - search range parameter . Beta must be between 1/n and 1.
#' @return Quantile of BJ statistics.
#' @seealso \code{\link{stat.bj}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#'
#' 3. Berk, R.H. & Jones, D.H. Z. "Goodness-of-fit test statistics that dominate the Kolmogorov statistics". Wahrscheinlichkeitstheorie verw Gebiete (1979) 47: 47.
#' @examples
#' ## The 0.05 critical value of BJ statistic when n = 10:
#' qbj(p=.95, n=10, beta=0.5)
#' @export
qbj <- function(p, n, beta){
  return(qphi(p, n, s=1, beta))
}
