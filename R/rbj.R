#' Random number generation of Berk and Jones (BJ) statitic.
#' @param k - number of observations.
#' @param n - dimension parameter, i.e. the number of input statitics to construct BJ statistic.
#' @param beta - search range parameter . Beta must be between 1/n and 1.
#' @return Random number from BJ statistics. The length of the result is determined by k for rbj.
#' @seealso \code{\link{stat.bj}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Jager, Leah; Wellner, Jon A. "Goodness-of-fit tests via phi-divergences". Annals of Statistics 35 (2007).
#' 3. Berk, R.H. & Jones, D.H. Z. "Goodness-of-fit test statistics that dominate the Kolmogorov statistics". Wahrscheinlichkeitstheorie verw Gebiete (1979) 47: 47.
#' @examples
#' ## Generate a vector of 10 random numbers from BJ statistic when n = 100:
#' rbj(k=10, n=100, beta=1)
#' @export
#' @importFrom stats runif
rbj <- function(k, n, beta){
  return(rphi(k, n, 1, beta))
}
