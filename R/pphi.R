########################################################################
#' calculate the left-tail probability of phi-divergence under general correlation matrix.
#' @param q - quantile, must be a scalar.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param s - the phi-divergence test parameter.
#' @param t - numerical truncation parameter.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @param method - default = "ecc": the effective correlation coefficient method in reference 2. "ave": the average method in reference 3, which is an earlier version of reference 2. The "ecc" method is more accurate and numerically stable than "ave" method.
#' @param ei - the eigenvalues of M if available. 
#' @return Left-tail probability of the phi-divergence statistics.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and power of optimal signal-detection statistics in finite case", IEEE Transactions on Signal Processing (2020) 68, 1021-1033
#' 2. Hong Zhang and Zheyang Wu. "The general goodness-of-fit tests for correlated data", Computational Statistics & Data Analysis (2022) 167, 107379
#' 3. Hong Zhang and Zheyang Wu. "Generalized Goodness-Of-Fit Tests for Correlated Data", arXiv:1806.03668.
#' @examples
#' M = toeplitz(1/(1:10)*(-1)^(0:9)) #alternating polynomial decaying correlation matrix
#' pphi(q=2, M=M, k0=1, k1=5, s=2)
#' pphi(q=2, M=M, k0=1, k1=5, s=2, method = "ecc")
#' pphi(q=2, M=M, k0=1, k1=5, s=2, method = "ave")
#' pphi(q=2, M=diag(10), k0=1, k1=5, s=2)
#' @export

pphi <- function(q, M, k0, k1, s=2, t=30, onesided=FALSE, method="ecc", ei = NULL)
{
  qtemp = q
  q = max(0.01, q)
  n = length(M[1,])
  nrep = n
  if(method=="ave"){
    rho = sample(M[row(M)!=col(M)], 2*nrep)
    rho = abs(rho)
    digi = rep(0.02, nrep)
    rho = ceiling(rho/digi)*digi
    rho[rho==1] = 0.99
    tbl = table(rho)
    rho_unique = as.numeric(names(tbl)) # unique rho
    rho_freq = as.numeric(tbl) # frequency of unique rho
    nrep_uniq = length(rho_unique)
  }else if(method=="ecc"){
    if(is.null(ei)){
      ei = abs(eigen(M, only.values = TRUE, symmetric = TRUE)$values)
    }
    rho_unique = (n/(n-1)*mean((abs(ei-1))^3)/(1+(n-1)^2))^(1/3)
    rho_freq = 1
  }
  
  pvalue_cal = sapply(rho_unique, function(x)pphi.rho(threshold=q, n=n, rho=x, k0=k0, k1=k1, s=s, t=t, onesided=onesided))
  res = sum(pvalue_cal*rho_freq)/sum(rho_freq)
  if(qtemp>=0.01){
    return(res)
  }else{
    warning(paste0("Left-tail prob. < ",round(res,3), " or p-value > ", 1-round(res,3)))
    return(res)
  }
}
