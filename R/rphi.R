#' Random number generation of phi-divergence statitics.
#' @param k - number of observations.
#' @param n - dimension parameter, i.e. the number of input statitics to construct phi-divergence statistic.
#' @param s - phi-divergence parameter. s = 2 is the higher criticism statitic.s = 1 is the Berk and Jones statistic.
#' @param beta - search range parameter . Beta must be between 1/n and 1.
#' @return Random number from phi-divergence statistics. The length of the result is determined by k for rphi.
#' @seealso \code{\link{stat.phi}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Jager, Leah; Wellner, Jon A. "Goodness-of-fit tests via phi-divergences". Annals of Statistics 35 (2007).
#' @examples
#' ## Generate a vector of 10 random numbers from Berk and Jones statistic when n = 100:
#' rphi(k=10, n=100, s=1, beta=1)
#' @export
#' @importFrom stats runif
rphi <- function(k, n, s, beta){
  return(unlist(lapply(1:k,function(x)stat.phi(runif(n), s, beta)$value)))
}
