#' CDF of Higher Criticism statitic under the null hypothesis.
#' @param q - quantile, must be a scalar.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param LS - if LS = T, then method of Li and Siegmund (2015) will be implemented (for independence case only).
#' @param ZW - if ZW = T, then approximation method of Zhang and Wu will be implemented.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @seealso \code{\link{stat.hc}} for the definition of the statistic.
#' @return The left-tail probability of the null distribution of HC statistic at the given quantile.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#'
#' 3. Li, Jian; Siegmund, David. "Higher criticism: p-values and criticism". Annals of Statistics 43 (2015).
#' @examples
#' pval <- runif(10)
#' hcstat <- stat.phi(pval, s=2, k0=1, k1=5)$value
#' phc(q=hcstat, M=diag(10), k0=1, k1=10)
#' @export
phc <- function(q, M, k0, k1, LS = F, ZW = F, onesided=FALSE){
  if(LS == F && ZW == F){
    return(pphi(q=q, M=M, k0=k0, k1=k1, s=2, onesided=onesided))
  }else{
    n = length(M[1,])
    beta = k1/n
    if(LS == T){
      return(1-hcls(q, n, beta))
    }else{
      return(1-hczw(q, n, beta))
    }
  }
}
