########################################################################
#' Multiple comparison test using phi-divergence statitics.
#' @param prob - vector of input p-values.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param s - phi-divergence parameter. s = 2 is the higher criticism statitic.s = 1 is the Berk and Jones statistic.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @return pvalue - The p-value of the phi-divergence test.
#' @return phistat - phi-diergence statistic.
#' @return location - the order of the input p-values to obtain phi-divergence statistic.
#' @seealso \code{\link{stat.phi}} for the definition of the statistic.v
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Jager, Leah; Wellner, Jon A. "Goodness-of-fit tests via phi-divergences". Annals of Statistics 35 (2007).
#' @examples
#' test.phi(runif(10), M=diag(10), s = 0.5, k0=1, k1=10)
#' #When the input are statistics#
#' stat.test = rnorm(20)
#' p.test = 2*(1 - pnorm(abs(stat.test)))
#' test.phi(p.test, M=diag(20), s = 0.5, k0=1, k1=10)
#' @export
test.phi <- function(prob, M, k0, k1, s=2, onesided=FALSE){
  n = length(prob)
  result = stat.phi(prob, s, k0, k1)
  return(list(pvalue=1-pphi(result$value, M=M, s=s, k0=k0, k1=k1, onesided=onesided), phistat=result$value, location=which(result$value==max(result$value))))
}
