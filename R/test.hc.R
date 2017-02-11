#' Multiple comparison test using Higher Criticism (HC) statitics.
#' @param prob - vector of input p-values.
#' @param beta - search range parameter . Beta must be between 1/n and 1.
#' @param LS - if LS = T, then method of Li and Siegmund (2015) will be implemented.When n and q is very large, approximation method is prefered.
#' @param ZW - if ZW = T, then approximation method of Zhang and Wu will be implemented.
#' @return pvalue - The p-value of the HC test.
#' @return hcstat - HC statistic.
#' @return location - the order of the input p-values to obtain HC statistic.
#' @seealso \code{\link{stat.hc}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#'
#' 3. Li, Jian; Siegmund, David. "Higher criticism: p-values and criticism". Annals of Statistics 43 (2015).
#' @examples
#' pval.test = runif(10)
#' test.hc(pval.test, 0.5)
#' test.hc(pval.test, 0.5, LS = TRUE)
#' test.hc(pval.test, 0.5, ZW = TRUE)
#' #When the input are statistics#
#' stat.test = rnorm(20)
#' p.test = 1 - pnorm(stat.test)
#' test.hc(p.test, beta = 0.5)
#' @export
test.hc <- function(prob, beta, LS = F, ZW = F){
  n = length(prob)
  result = stat.phi(prob, 2, beta)
  if(LS == F && ZW == F){
    return(list(pvalue=1-pphi(result$value, n, 2, beta), hcstat=result$value, location=which(result$value==max(result$value))))
  }else{
    if(LS == T){
      return(list(pvalue=hcls(result$value, n, beta), hcstat=result$value, location=which(result$value==max(result$value))))
    }else{
      return(list(pvalue=hczw(result$value, n, beta), hcstat=result$value, location=which(result$value==max(result$value))))
    }
  }
}
