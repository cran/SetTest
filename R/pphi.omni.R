########################################################################
#' calculate the left-tail probability of omnibus phi-divergence statistics under general correlation matrix.
#' @param q - quantile, must be a scalar.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param K0 - vector of search range starts (from the k0th smallest p-value).
#' @param K1 - vector of search range ends (at the k1th smallest p-value).
#' @param S - vector of the phi-divergence test parameters.
#' @param t - numerical truncation parameter.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#' @examples
#' M = matrix(0.3,10,10) + diag(1-0.3, 10)
#' pphi.omni(0.05, M=M, K0=rep(1,4), K1=rep(5,4), S=c(-1,0,1,2))
#' @export
#' @importFrom stats loess predict
pphi.omni <- function(q,M, K0, K1, S, t=30, onesided=FALSE)
{
  n = length(M[1,])
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
}

