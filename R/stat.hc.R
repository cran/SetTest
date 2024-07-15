#' Construct Higher Criticism (HC) statistics.
#' @param p - vector of input p-values.
#' @param k0 - search range left end parameter. Default k0 = 1.
#' @param k1 - search range right end parameter. Default k1 = 0.5*number of input p-values.
#' @return value - HC statistic constructed from a vector of p-values.
#' @return location - the order of the p-values to obtain HC statistic.
#' @return stat - vector of marginal HC statistics.
#' @details Let \eqn{p_{(i)}}, \eqn{i = 1,...,n} be a sequence of ordered p-values, the higher criticism statistic
#' \deqn{HC = \sqrt{n} \max_{1 \leq i\leq \lfloor \beta n \rfloor} [i/n - p_{(i)}] /\sqrt{p_{(i)}(1 - p_{(i)})}}
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#' @examples
#' stat.hc(runif(10))
#' #When the input are statistics#
#' stat.test = rnorm(20)
#' p.test = 1 - pnorm(stat.test)
#' stat.hc(p.test, k0 = 1, k1 = 10)
#' @export
stat.hc <- function(p, k0=1, k1=NA){
  result = stat.phi(p, 2, k0, k1)
  return(list(value=result$value, location=result$location, stat=result$stat))
}
