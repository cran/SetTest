#' Multiple comparison test using phi-divergence statitics.
#' @param prob - vector of input p-values.
#' @param s - phi-divergence parameter. s = 2 is the higher criticism statitic.s = 1 is the Berk and Jones statistic.
#' @param beta - search range parameter . Beta must be between 1/n and 1.
#' @return pvalue - The p-value of the phi-divergence test.
#' @return phistat - phi-divergence statistic.
#' @return location - the order of the input p-values to obtain phi-divergence statistic.
#' @seealso \code{\link{stat.phi}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Jager, Leah; Wellner, Jon A. "Goodness-of-fit tests via phi-divergences". Annals of Statistics 35 (2007).
#' @examples
#' test.phi(runif(10), s = 0.5, beta = 0.5)
#' #When the input are statistics#
#' stat.test = rnorm(20)
#' p.test = 1 - pnorm(stat.test)
#' test.phi(p.test, s = 0.5,beta = 0.5)
#' @export
test.phi <- function(prob, s, beta){
  n = length(prob)
  result = stat.phi(prob, s, beta)
  return(list(pvalue=1-pphi(result$value, n, s, beta), phistat=result$value, location=which(result$value==max(result$value))))
}
