########################################################################
#' calculate the left-tail probability of omnibus phi-divergence statistics under general correlation matrix.
#' @param q - quantile, must be a scalar.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param K0 - vector of search range starts (from the k0th smallest p-value).
#' @param K1 - vector of search range ends (at the k1th smallest p-value).
#' @param S - vector of the phi-divergence test parameters.
#' @param t - numerical truncation parameter.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @param method - default = "ecc": the effective correlation coefficient method in reference 2. "ave": the average method in reference 3, which is an earlier version of reference 2. The "ecc" method is more accurate and numerically stable than "ave" method.
#' @param ei - the eigenvalues of M if available. 
#' @return Left-tail probability of omnibus phi-divergence statistics.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and power of optimal signal-detection statistics in finite case", IEEE Transactions on Signal Processing (2020) 68, 1021-1033
#' 2. Hong Zhang and Zheyang Wu. "The general goodness-of-fit tests for correlated data", Computational Statistics & Data Analysis (2022) 167, 107379
#' 3. Hong Zhang and Zheyang Wu. "Generalized Goodness-Of-Fit Tests for Correlated Data", arXiv:1806.03668.
#' @examples
#' M = matrix(0.3,10,10) + diag(1-0.3, 10)
#' pphi.omni(0.05, M=M, K0=rep(1,2), K1=rep(5,2), S=c(1,2))
#' @export
#' @importFrom stats loess predict
pphi.omni <- function(q, M, K0, K1, S, t=30, onesided=FALSE, method="ecc", ei = NULL)
{
  n = length(M[1,])
  if(method=="ave"){
    nrep = n
    nthr=5
    lower = max(1e-12,q*0.9); upper = min(1,q*1.1)
    thr = seq(lower,upper,length.out=nthr)
    rho = sample(M[row(M)!=col(M)], nrep)
    rho = abs(rho)
    
    digi = rep(0.02, nrep)
    rho = round(rho/digi)*digi # rho round to the nearest 0.02
    tbl = table(rho)
    rho_unique = as.numeric(names(tbl)) # unique rho
    rho_freq = as.numeric(tbl) # frequency of unique rho
    
    nrep_uniq = length(rho_unique)
    pvalue_cal = c()
    for(i in 1:nrep_uniq){
      freq = rho_freq[i]
      y = sapply(1:nthr, function(x)pphi.rho.omni(threshold=thr[x], n=n, rho=rho_unique[i], K0=K0, K1=K1, S=S, t=t, onesided=onesided))
      pvalue_cal = append(pvalue_cal,rep(y, freq))
    }
    
    data = data.frame(y=pvalue_cal, x=rep(thr, nrep))
    fit = loess(y~x, data = data)
    return(predict(fit, q))
  }else if(method=="ecc"){
    if(is.null(ei)){
      ei = abs(eigen(M, only.values = TRUE, symmetric = TRUE)$values)
    }
    rho_unique = (n/(n-1)*mean((abs(ei-1))^3)/(1+(n-1)^2))^(1/3)
    return(pphi.rho.omni(threshold=q, n=n, rho=rho_unique, K0=K0, K1=K1, S=S, t=t, onesided=onesided))
  }
}

