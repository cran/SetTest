#' Statistical power of Higher Criticism test.
#' @param alpha - type-I error rate.
#' @param n - dimension parameter, i.e. the number of input statitics to construct Higher Criticism statistic.
#' @param beta - search range parameter. Search range = (1, beta*n). Beta must be between 1/n and 1.
#' @param method - different alternative hypothesis, including mixtures such as, "gaussian-gaussian", "gaussian-t", "t-t", "chisq-chisq", and "exp-chisq". By default, we use Gaussian mixture.
#' @param eps - mixing parameter of the mixture.
#' @param mu - mean of non standard Gaussian model.
#' @param df - degree of freedom of t/Chisq distribution and exp distribution.
#' @param delta - non-cenrality of t/Chisq distribution.
#' @return Power of HC test.
#' @details We consider the following hypothesis test,
#' \deqn{H_0: X_i\sim F, H_a: X_i\sim G}
#' Specifically, \eqn{F = F_0} and \eqn{G = (1-\epsilon)F_0+\epsilon F_1}, where \eqn{\epsilon} is the mixing parameter, \eqn{F_0} and \eqn{F_1} is
#' speified by the "method" argument:
#'
#' "gaussian-gaussian": \eqn{F_0} is the standard normal CDF and \eqn{F = F_1} is the CDF of normal distribution with \eqn{\mu} defined by mu and \eqn{\sigma = 1}.
#'
#' "gaussian-t": \eqn{F_0} is the standard normal CDF and \eqn{F = F_1} is the CDF of t distribution with degree of freedom defined by df.
#'
#' "t-t": \eqn{F_0} is the CDF of t distribution with degree of freedom defined by df and \eqn{F = F_1} is the CDF of non-central t distribution with degree of freedom defined by df and non-centrality defined by delta.
#
#' "chisq-chisq": \eqn{F_0} is the CDF of Chisquare distribution with degree of freedom defined by df and \eqn{F = F_1} is the CDF of non-central Chisquare distribution with degree of freedom defined by df and non-centrality defined by delta.
#'
#' "exp-chisq": \eqn{F_0} is the CDF of exponential distribution with parameter defined by df and \eqn{F = F_1} is the CDF of non-central Chisqaure distribution with degree of freedom defined by df and non-centrality defined by delta.
#' @seealso \code{\link{stat.hc}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#' @examples
#' power.hc(0.05, n=10, beta=0.5, eps = 0.1, mu = 1.2)
#' @export
#' @importFrom stats dpois pbeta
power.hc <- function(alpha, n, beta, method = "gaussian-gaussian",
                      eps=0, mu=0, df=1, delta=0){
  return(power.phi(alpha, n, 2, beta, method, eps, mu, df, delta))
}
