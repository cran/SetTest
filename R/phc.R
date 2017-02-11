#' CDF of Higher Criticism statitic under the null hypothesis.
#' @param q - quantile, must be a scalar.
#' @param n - dimension parameter, i.e. the number of input statitics to construct HC statistic.
#' @param beta - search range parameter . Beta must be between 1/n and 1.
#' @param LS - if LS = T, then method of Li and Siegmund (2015) will be implemented.When n and q is very large, approximation method is prefered.
#' @param ZW - if ZW = T, then approximation method of Zhang and Wu will be implemented.
#' @seealso \code{\link{stat.hc}} for the definition of the statistic.
#' @return The left-tail probability of the null distribution of HC statistic at the given quantile.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal
#' Signal Detection Methods in Finite Samples", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#'
#' 3. Li, Jian; Siegmund, David. "Higher criticism: p-values and criticism". Annals of Statistics 43 (2015).
#' @examples
#' pval <- runif(10)
#' hcstat <- stat.phi(pval, 2, 0.5)$value
#' phc(q=hcstat, n=10, beta=0.5)
#' @export
phc <- function(q, n, beta, LS = F, ZW = F){
  if(LS == F && ZW == F){
    return(pphi(q, n, 2, beta))
  }else{
    if(LS == T){
      return(1-hcls(q, n, beta))
    }else{
      return(1-hczw(q, n, beta))
    }
  }
}
