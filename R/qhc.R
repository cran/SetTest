#' Quantile of Higher Criticism statitics under the null hypothesis.
#' @param p -  a scalar left probability that defines the quantile.
#' @param M - correlation matrix of input statistics (of the input p-values).
#' @param k0 - search range starts from the k0th smallest p-value.
#' @param k1 - search range ends at the k1th smallest p-value.
#' @param onesided - TRUE if the input p-values are one-sided.
#' @param LS - if LS = T, then method of Li and Siegmund (2015) will be implemented.When n and q is very large, approximation method is prefered.
#' @param ZW - if ZW = T, then approximation method of Zhang and Wu will be implemented.
#' @return Quantile of HC statistics.
#' @seealso \code{\link{stat.hc}} for the definition of the statistic.
#' @references 1. Hong Zhang, Jiashun Jin and Zheyang Wu. "Distributions and Statistical Power of Optimal Signal-Detection Methods In Finite Cases", submitted.
#'
#' 2. Donoho, David; Jin, Jiashun. "Higher criticism for detecting sparse heterogeneous mixtures". Annals of Statistics 32 (2004).
#'
#' 3. Li, Jian; Siegmund, David. "Higher criticism: p-values and criticism". Annals of Statistics 43 (2015).
#' @examples
#' ## The 0.05 critical value of HC statistic when n = 10:
#' qhc(p=.95, M=diag(10), k0=1, k1=5, onesided=FALSE)
#' @export
qhc <- function(p, M, k0, k1, onesided=FALSE, LS = F, ZW = F){
  if(LS == F && ZW == F){
    return(qphi(p, M, k0, k1, s=2, onesided))
  }else{
    n = length(M[1,])
    if(p < 0.5){
      stop("probability is too small, approximation is not accurate")
    }else{
      if(LS == T){
        lower = 2
        upper = 10
        maxnumIter = 500
        numIter = 1
        p_cal_lower = 1- hcls(lower, n, beta)
        p_cal_upper = 1- hcls(upper, n, beta)
        while((p_cal_lower-p)*(p_cal_upper-p)>0 && numIter < maxnumIter){
          if(p_cal_lower>p){
            lower = lower*0.9
            p_cal_lower = 1- hcls(lower, n, beta)
          }
          if(p_cal_upper < p){
            upper = upper*1.1
            p_cal_upper = 1- hcls(upper, n, beta)
          }
          numIter = numIter + 1
        }
        q = mean(c(lower,upper))
        p_cal = 1- hcls(q, n, beta)
        error = (p_cal-p)/p
        numIter = 1
        while(abs(error) > 1e-5 && numIter < maxnumIter){
          if(error > 0){
            upper = q
            q = mean(c(lower,upper))
            p_cal = 1- hcls(q, n, beta)
            error = (p_cal-p)/p
          }else{
            lower = q
            q = mean(c(lower,upper))
            p_cal = 1- hcls(q, n, beta)
            error = (p_cal-p)/p
          }
          numIter = numIter + 1
        }
        return(q)
      }else{
        lower = 2
        upper = 10
        maxnumIter = 500
        numIter = 1
        p_cal_lower = 1- hczw(lower, n, beta)
        p_cal_upper = 1- hczw(upper, n, beta)
        while((p_cal_lower-p)*(p_cal_upper-p)>0 && numIter < maxnumIter){
          if(p_cal_lower>p){
            lower = lower*0.9
            p_cal_lower = 1- hczw(lower, n, beta)
          }
          if(p_cal_upper < p){
            upper = upper*1.1
            p_cal_upper = 1- hczw(upper, n, beta)
          }
          numIter = numIter + 1
        }
        q = mean(c(lower,upper))
        p_cal = 1- hczw(q, n, beta)
        error = (p_cal-p)/p
        numIter = 1
        while(abs(error) > 1e-5 && numIter < maxnumIter){
          if(error > 0){
            upper = q
            q = mean(c(lower,upper))
            p_cal = 1- hczw(q, n, beta)
            error = (p_cal-p)/p
          }else{
            lower = q
            q = mean(c(lower,upper))
            p_cal = 1- hczw(q, n, beta)
            error = (p_cal-p)/p
          }
          numIter = numIter + 1
        }
        return(q)
      }
    }
  }
}
