#' Construct phi-divergence statistics.
#' @param p - vector of input p-values.
#' @param s - phi-divergence parameter. s = 2 is the higher criticism statitic.s = 1 is the Berk and Jones statistic.
#' @param k0 - search range left end parameter. Default k0 = 1.
#' @param k1 - search range right end parameter. Default k1 = 0.5*number of input p-values.
#' @return value - phi-divergence statistic constructed from a vector of p-values.
#' @return location - the order of the p-values to obtain phi-divergence statistic.
#' @return stat - vector of marginal phi-divergence statistics.
#' @details Let \eqn{p_{(i)}}, \eqn{i = 1,...,n} be a sequence of ordered p-values, the phi-divergence statistic
#' \deqn{PHI = \sqrt{2n}/(s - s^2) \max_{1 \leq i\leq \lfloor \beta n \rfloor} (-1)^j \sqrt{1 - (i/n)^s  (p_{(i)})^s - (1-i/n)^{(1-s)} * (1-p_{(i)})^{(1-s)}}}
#' and when \eqn{p_{(i)} > i/n}, \eqn{j = 1}, otherwise \eqn{j = 0}.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Jager, Leah; Wellner, Jon A. "Goodness-of-fit tests via phi-divergences". Annals of Statistics 35 (2007).
#' @examples
#' stat.phi(runif(10), s = 2)
#' #When the input are statistics#
#' stat.test = rnorm(20)
#' p.test = 1 - pnorm(stat.test)
#' stat.phi(p.test, s = 0.5, k0 = 2, k1 = 5)
#' @export
stat.phi <- function (p, s, k0=1, k1=NA) {
  n = length(p)
  if(is.na(k1)){
    m = floor(0.5*n)
  }else{
    m = floor(k1)
  }

  p = sort(p)[k0:m]
  thr = ((k0:m)/n)
  stat = (2 * (p < thr) - 1) * sqrt(2 * n * phi.f(thr, p, s))
  return(list(value = max(stat), location = which.max(stat),
              stat = stat))
}
