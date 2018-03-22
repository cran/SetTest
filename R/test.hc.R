#' Multiple comparison test using Higher Criticism (HC) statitics.
#' @param prob - vector of input p-values.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param LS - if LS = T, then method of Li and Siegmund (2015) will be implemented.When n and q is very large, approximation method is prefered.
#' @param ZW - if ZW = T, then approximation method of Zhang and Wu will be implemented.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @return pvalue - The p-value of the HC test.
#' @return hcstat - HC statistic.
#' @return location - the order of the input p-values to obtain HC statistic.
#' @seealso \code{\link{stat.hc}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#'
#' 3. Li, Jian; Siegmund, David. "Higher criticism: p-values and criticism". Annals of Statistics 43 (2015).
#' @examples
#' pval.test = runif(10)
#' test.hc(pval.test, M=diag(10), k0=1, k1=10)
#' test.hc(pval.test, M=diag(10), k0=1, k1=10, LS = TRUE)
#' test.hc(pval.test, M=diag(10), k0=1, k1=10, ZW = TRUE)
#' #When the input are statistics#
#' stat.test = rnorm(20)
#' p.test = 2*(1 - pnorm(abs(stat.test)))
#' test.hc(p.test, M=diag(20), k0=1, k1=10)
#' @export
test.hc <- function(prob, M, k0, k1, LS = F, ZW = F, onesided=FALSE){
  n = length(prob)
  result = stat.phi(prob, 2, k0, k1)
  if(LS == F && ZW == F){
    return(list(pvalue=1-pphi(result$value, M=M, s=2, k0=k0, k1=k1, onesided=onesided), hcstat=result$value, location=which(result$value==max(result$value))))
  }else{
    beta = k1/n
    if(LS == T){
      return(list(pvalue=hcls(result$value, n, beta), hcstat=result$value, location=which(result$value==max(result$value))))
    }else{
      return(list(pvalue=hczw(result$value, n, beta), hcstat=result$value, location=which(result$value==max(result$value))))
    }
  }
}
