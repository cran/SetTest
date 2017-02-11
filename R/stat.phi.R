#' Construct phi-divergence statitics.
#' @param p - vector of input p-values.
#' @param s - phi-divergence parameter. s = 2 is the higher criticism statitic.s = 1 is the Berk and Jones statistic.
#' @param beta - search range parameter . Beta must be between 1/n and 1.
#' @return value - phi-divergence statistic constructed from a vector of p-values.
#' @return location - the order of the p-values to obtain phi-divergence statistic.
#' @return stat - vector of marginal phi-divergence statistics.
#' @details Let \eqn{p_{(i)}}, \eqn{i = 1,...,n} be a sequence of ordered p-values, the phi-divergence statistic
#' \deqn{PHI = \sqrt{2n}/(s - s^2) \max_{1 \leq i\leq \lfloor \beta n \rfloor} (-1)^j \sqrt{1 - (i/n)^s  (p_{(i)})^s - (1-i/n)^{(1-s)} * (1-p_{(i)})^{(1-s)}}}
#' and when \eqn{p_{(i)} > i/n}, \eqn{j = 1}, otherwise \eqn{j = 0}.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Jager, Leah; Wellner, Jon A. "Goodness-of-fit tests via phi-divergences". Annals of Statistics 35 (2007).
#' @examples
#' stat.phi(runif(10), 2, 0.5)
#' #When the input are statistics#
#' stat.test = rnorm(20)
#' p.test = 1 - pnorm(stat.test)
#' stat.phi(p.test, s = 0.5, beta = 0.5)
#' @export
stat.phi <- function(p, s, beta){
  n = length(p)
  m = floor(beta*n)
  p = sort(p)[1:m] # if m<1, it automatically optains minimum
  thr = ((1:m)/n)
  stat = (2*(p < thr)-1)*sqrt(2*n*phi.f(thr, p, s))
  #return(max(stat))
  return(list(value=max(stat), location=which(stat==max(stat)), stat=stat))
}
