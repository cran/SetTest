########################################################################
#' calculate the omnibus phi-divergence statistics under general correlation matrix.
#' @param p - input pvalues.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param K0 - vector of search range starts (from the k0th smallest p-value).
#' @param K1 - vector of search range ends (at the k1th smallest p-value).
#' @param S - vector of the phi-divergence test parameters.
#' @param t - numerical truncation parameter.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#' @examples
#' M = toeplitz(1/(1:10)*(-1)^(0:9)) #alternating polynomial decaying correlation matrix
#' stat.phi.omni(runif(10), M=M, K0=rep(1,4), K1=rep(5,4), S=c(-1,0,1,2))
#' @export
stat.phi.omni <- function(p, M, K0=rep(1,4), K1=rep(length(M[1,]),4), S=c(-1,0,1,2), t=30, onesided=FALSE){
  #n = length(p)
  phistat = unlist(mapply(function(x,y,z) stat.phi(p=p,x,y,z), S, K0, K1)[1,])
  pval = mapply(function(x,z1,z2,y) 1-pphi(x,M=M,z1,z2,y,t=t,onesided=onesided), x=phistat, z1=K0, z2=K1, y=S)
  if(all(is.na(pval))){
    return(list(minp=1,pval=pval))
  }else{
    return(list(minp=min(pval,na.rm=TRUE), pval=pval))
  }
}
