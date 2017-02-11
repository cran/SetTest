#' Random number generation of Higher Criticism (HC) statitic.
#' @param k - number of observations.
#' @param n - dimension parameter, i.e. the number of input statitics to construct HC statistic.
#' @param beta - search range parameter . Beta must be between 1/n and 1.
#' @seealso \code{\link{stat.hc}} for the definition of the statistic.
#' @return Random number from HC statistics. The length of the result is determined by k for rhc.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#' @examples
#' ## Generate a vector of 10 random numbers from HC statistic when n = 100:
#' rhc(k=10, n=100, beta=1)
#' @export
#' @importFrom stats runif
rhc <- function(k, n, beta){
  return(rphi(k, n, 2, beta))
}
