#' Construct Berk and Jones (BJ) statitics.
#' @param p - vector of input p-values.
#' @param beta - search range parameter . Beta must be between 1/n and 1.
#' @return value - BJ statistic constructed from a vector of p-values.
#' @return location - the order of the p-values to obtain BJ statistic.
#' @return stat - vector of marginal BJ statistics.
#' @details Let \eqn{p_{(i)}}, \eqn{i = 1,...,n} be a sequence of ordered p-values, the Berk and Jones statistic
#' \deqn{BJ = \sqrt{2n} \max_{1 \leq i\leq \lfloor \beta n \rfloor} (-1)^j \sqrt{i/n * \log(i/n/p_{(i)}) + (1-i/n) * \log((1-i/n)/(1-p_{(i)}))}}
#' and when \eqn{p_{(i)} > i/n}, \eqn{j = 1}, otherwise \eqn{j = 0}.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Jager, Leah; Wellner, Jon A. "Goodness-of-fit tests via phi-divergences". Annals of Statistics 35 (2007).
#'
#' 3. Berk, R.H. & Jones, D.H. Z. "Goodness-of-fit test statistics that dominate the Kolmogorov statistics". Wahrscheinlichkeitstheorie verw Gebiete (1979) 47: 47.
#' @examples
#' stat.bj(runif(10), beta = 0.3)
#' #When the input are statistics#
#' stat.test = rnorm(20)
#' p.test = 1 - pnorm(stat.test)
#' stat.bj(p.test, beta = 0.5)
#' @export
stat.bj <- function(p, beta = 0.5){
  result = stat.phi(p, 1, beta)
  return(list(value=result$value, location=result$location, stat=result$stat))
}
