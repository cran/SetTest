########################################################################
#' calculate the right-tail probability of omnibus phi-divergence statistics under general correlation matrix.
#' @param prob - vector of input p-values.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param K0 - vector of search range starts (from the k0th smallest p-value).
#' @param K1 - vector of search range ends (at the k1th smallest p-value).
#' @param S - vector of the phi-divergence test parameters.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#' @examples
#' M = matrix(0.3,10,10) + diag(1-0.3, 10)
#' test.phi.omni(runif(10), M=M, K0=rep(1,4), K1=rep(5,4), S=c(-1,0,1,2))
#' @export

test.phi.omni <- function(prob, M, K0, K1, S, onesided=FALSE){
  result = stat.phi.omni(p=prob, M=M, K0=K0, K1=K1, S=S, onesided=onesided)
  return(list(pvalue=pphi.omni(result$minp, M=M, S=S, K0=K0, K1=K1, onesided=onesided), indiv_pval=result$pval))
}
