########################################################################
#' calculate the left-tail probability of phi-divergence under general correlation matrix.
#' @param q - quantile, must be a scalar.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param s - the phi-divergence test parameter.
#' @param t - numerical truncation parameter.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#' @examples
#' M = toeplitz(1/(1:10)*(-1)^(0:9)) #alternating polynomial decaying correlation matrix
#' pphi(q=2, M=M, k0=1, k1=5, s=2)
#' pphi(q=2, M=diag(10), k0=1, k1=5, s=2)
#' @export
#' @importFrom stats loess predict

pphi <- function(q, M, k0, k1, s=2, t=30, onesided=FALSE)
{
  qtemp = q
  q = max(0.01, q)
  NODES=seq(-4,4,length.out=40)
  n = length(M[1,])
  nrep = n
  nthr = 5
  lower = q*0.9#max(0.5,q*0.9)
  upper = q*1.1
  thr = seq(lower,upper,length.out=nthr)
  rho = sample(M[row(M)!=col(M)], nrep)
  rho = abs(rho)
  digi = rep(0.02, nrep)
  rho = round(rho/digi)*digi
  rho[rho==1] = 0.99
  tbl = table(rho)
  rho_unique = as.numeric(names(tbl)) # unique rho
  rho_freq = as.numeric(tbl) # frequency of unique rho

  nrep_uniq = length(rho_unique)
  pvalue_cal = c()
  for(i in 1:nrep_uniq){
    freq = rho_freq[i]
    y = sapply(1:nthr, function(x)pphi.rho(threshold=thr[x], n=n, rho=rho_unique[i], k0=k0, k1=k1, NODES=NODES, s=s, t=t, onesided=onesided))
    pvalue_cal = append(pvalue_cal,rep(y, freq))
  }

  data = data.frame(y=pvalue_cal, x=rep(thr, nrep))
  fit = loess(y~x, data = data)
  res = predict(fit, q)
  if(qtemp>=0.01){
    return(res)
  }else{
    warning(paste0("Left-tail prob. < ",round(res,3), " or p-value > ", 1-round(res,3)))
    return(res)
  }
}
