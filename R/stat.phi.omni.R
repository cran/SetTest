########################################################################
#' calculate the omnibus phi-divergence statistics under general correlation matrix.
#' @param p - input pvalues.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param K0 - vector of search range starts (from the k0th smallest p-value).
#' @param K1 - vector of search range ends (at the k1th smallest p-value).
#' @param S - vector of the phi-divergence test parameters.
#' @param t - numerical truncation parameter.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @param method - default = "ecc": the effective correlation coefficient method in reference 2. "ave": the average method in reference 3, which is an earlier version of reference 2. The "ecc" method is more accurate and numerically stable than "ave" method.
#' @param ei - the eigenvalues of M if available. 
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and power of optimal signal-detection statistics in finite case", IEEE Transactions on Signal Processing (2020) 68, 1021-1033
#' 2. Hong Zhang and Zheyang Wu. "The general goodness-of-fit tests for correlated data", Computational Statistics & Data Analysis (2022) 167, 107379
#' 3. Hong Zhang and Zheyang Wu. "Generalized Goodness-Of-Fit Tests for Correlated Data", arXiv:1806.03668.
#' @examples
#' p.test = runif(10)
#' M = toeplitz(1/(1:10)*(-1)^(0:9)) #alternating polynomial decaying correlation matrix
#' stat.phi.omni(p.test, M=M, K0=rep(1,2), K1=rep(5,2), S=c(1,2))
#' @export
stat.phi.omni <- function(p, M, K0=rep(1,2), K1=rep(length(M[1,]),2), S=c(1,2), t=30, onesided=FALSE, method="ecc", ei = NULL){
  phistat = unlist(mapply(function(x,y,z) stat.phi(p=p,x,y,z), S, K0, K1)[1,])
  pval = mapply(function(x,z1,z2,y) 1-pphi(x,M=M,z1,z2,y,t=t,onesided=onesided, method=method, ei=ei), x=phistat, z1=K0, z2=K1, y=S)
  if(all(is.na(pval))){
    return(list(minp=1,pval=pval))
  }else{
    return(list(minp=min(pval,na.rm=TRUE), pval=pval))
  }
}
