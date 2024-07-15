#' @importFrom stats uniroot dbinom pgamma dnorm pnorm pbeta dpois
# Inverse phi function
phi.f.inv.i <- function(q, i, n, s){
  ep = 10^(-15)
  if(i!=n){
    if(q>0){
      if (s!=1&&s!=0){
        f <- function(x)(1-(i/n)^s*x^(1-s)-(1-i/n)^(s)*(1-x)^(1-s))*2*n/(s-s^2) - q^2
        if (f(ep)<0){
          CC=0
        }else{
          CC=uniroot(f,c(ep,i/n),tol=1e-100)$root
        }
      }else if (s==1){
        f <- function(x)2*n*((i/n)*log(i/n/x)+(1-i/n)*log((1-i/n)/(1-x)))-q^2
        if (f(ep)<0){
          CC=0
        }else{
          CC=uniroot(f,c(ep,i/n),tol=1e-100)$root
        }
      }else{
        f <- function(x)2*n*((x)*log(x/(i/n))+(1-x)*log((1-x)/(1-i/n)))-q^2
        if (f(ep)<0){
          CC=0
        }
        else{
          CC=uniroot(f,c(ep,i/n),tol=1e-100)$root
        }
      }
    }else if(q<0){
      if (s!=1&&s!=0){
        f <- function(x)(1-(i/n)^s*x^(1-s)-(1-i/n)^(s)*(1-x)^(1-s))*2*n/(s-s^2) - q^2
        if (f(1-ep)<0){
          CC=1
        }else{
          CC=uniroot(f,c(i/n, 1-ep),tol=1e-100)$root
        }
      }else if (s==1){
        f <- function(x)2*n*((i/n)*log(i/n/x)+(1-i/n)*log((1-i/n)/(1-x)))-q^2
        if (f(1-ep)<0){
          CC=1
        }else{
          CC=uniroot(f,c(i/n, 1-ep),tol=1e-100)$root
        }
      }else{
        f <- function(x)2*n*((x)*log(x/(i/n))+(1-x)*log((1-x)/(1-i/n)))-q^2
        if (f(1-ep)<0){
          CC=1
        }else{
          CC=uniroot(f,c(i/n, 1-ep),tol=1e-100)$root
        }
      }
    }else{
      CC = i/n
    }
  }else{
    if(q<=0||s<=0){
      CC = 1
    }else{
      if(s==1){
        f <- function(x)-2*log(x)*n - q^2
        CC=uniroot(f,c(ep,1),tol=1e-100)$root
      }else{
        f <- function(x)(1-x^(1-s))*2*n/(s-s^2) - q^2
        CC=uniroot(f,c(ep,1),tol=1e-100)$root
      }
    }
  }
  return(CC)
}

# phi function
phi.f <- function(u, v, s){
  # if(any(sapply(u, function(x) x > 0 && x < 1)==F)||any(sapply(v, function(x) x > 0 && x < 1)==F)){
  #   warning("The 0 or 1 of u, v have been replaced by 1e-15 and 1-1e-15")
  # }
  u[u==0] = 1e-15
  u[u==1] = 1 - 1e-15
  v[v==0] = 1e-15
  v[v==0] = 1 - 1e-15
  if(s!=0&&s!=1){
    return((1-u^s*v^(1-s)-(1-u)^s*(1-v)^(1-s))/(s-s^2))
  }else if(s==1){
    return(u*log(u/v)+(1-u)*log((1-u)/(1-v)))
  }else{
    return(v*log(v/u)+(1-v)*log((1-v)/(1-u)))
  }
}

####################################################
#Li & Siegmund's approximation to HC               #
#Arguments: threshold: threshold                   #
#                beta: search range is at most beta#
#                   n: number of data              #
####################################################
hcls <- function(q,n,beta){
  PVALUE=0
  for (i in 1:(beta*n-1)){
    C=(n/(n+q^2))*(i/n+0.5*(q^2/n-q/sqrt(n)*sqrt(q^2/n+4*i*(n-i)/n^2)))
    Cprime=(n/(n+q^2))*(1-q/sqrt(n)*(1-2*i/n)/sqrt(q^2/n+4*i*(n-i)/n^2))
    A=(1-((n-i+1)*Cprime)/(n-n*C))*dbinom(i,n,C)
    PVALUE = PVALUE + A
  }
  return(PVALUE)
}


####################################################
#Zhang & Wu's approximation to HC                  #
#Arguments: threshold: threshold                   #
#                beta: search range is at most beta#
#                   n: number of data              #
####################################################
h <- function(x,k) x*pgamma(k*x,k-1,1)-pgamma(k*x,k,1)
hczw <- function(threshold,n,beta){
  PVALUE=0
  for (i in 1:(beta*n-1)){
    C=(n/(n+threshold^2))*(i/n+0.5*(threshold^2/n-threshold/sqrt(n)*sqrt(threshold^2/n+4*i*(n-i)/n^2)))
    Cprime=(n/(n+threshold^2))*(1-threshold/sqrt(n)*(1-2*i/n)/sqrt(threshold^2/n+4*i*(n-i)/n^2))
    A=(1-(n+1)/n*Cprime + h((n+1)/n*Cprime,min(beta*n-i,sqrt(n))))*dpois(i,(n+1)*C)
    PVALUE = PVALUE + A
  }
  return(PVALUE + pgamma((n+1)*C,beta*n,1))
}

########################################################################
# calculate the left-tail probability of phi-divergence under positive equal correlation.
# threshold - threshold.
# n - number of p-values.
# rho - common correlation coefficient.
# k0 - search range starts from the k0th smallest p-value.
# k1 - search range ends at the k1th smallest p-value.
# NODES - the numerical integration domain (under the standard normal pdf).
# s - the phi-divergence test parameter.
# t - numerical truncation parameter.
# onesided - TRUE if the input p-values are one-sided.
# normdn - dnorm(NODES)
pphi.rho <- function(threshold, n, rho, k0, k1, s=2, t=30, onesided=F, NODES=seq(-8,8,length.out=50),normdn=NULL){
  nNODES = length(NODES)
  t = min(t,k1-k0+1)
  u = sapply(1:k1,function(x)phi.f.inv.i(threshold, x, n, s))
  um = u[length(u)]
  if(is.null(normdn)){
    if(rho!=0){
      SUM = sum(dnorm(NODES))
      result = 0
      if(onesided==T){
        qu = qnorm(1-u)
        qum = qnorm(1-um)
        sqrho = sqrt(rho)
        sq_rho = sqrt(1-rho)
        if(length(qu)==1){
          D = t(as.matrix(sapply(1:nNODES, function(x)pnorm((qu-sqrho*NODES[x])/sq_rho,lower.tail = F))))
        }else{
          D = sapply(1:nNODES, function(x)pnorm((qu-sqrho*NODES[x])/sq_rho,lower.tail = F))
        }
        
        DM = sapply(1:nNODES, function(x)pnorm((qum-sqrho*NODES[x])/sq_rho,lower.tail = F))
        result = sum(dnorm(NODES)*sapply(1:nNODES, function(x)UnifCross(D[,x],DM[x],t=t,n=n,k0=k0,k1=k1)))/SUM

      }else{
        qu = qnorm(1-u/2)
        qum = qnorm(1-um/2)
        sqrho = sqrt(rho)
        sq_rho = sqrt(1-rho)
        if(length(qu)==1){
          D = t(as.matrix(sapply(1:nNODES, function(x)pnorm((qu-sqrho*NODES[x])/sq_rho,lower.tail = F) +
                                   pnorm((-qu-sqrho*NODES[x])/sq_rho))))
        }else{
          D = sapply(1:nNODES, function(x)pnorm((qu-sqrho*NODES[x])/sq_rho,lower.tail = F) +
                       pnorm((-qu-sqrho*NODES[x])/sq_rho))
        }
        
        DM = sapply(1:nNODES, function(x)pnorm((qum-sqrho*NODES[x])/sq_rho,lower.tail = F) +
                      pnorm((-qum-sqrho*NODES[x])/sq_rho))
        
        result = sum(dnorm(NODES)*sapply(1:nNODES, function(x)UnifCross(D[,x],DM[x],t=t,n=n,k0=k0,k1=k1)))/SUM
        
      }
      return(result)
    }else{
      return(UnifCross(u,um,t=t,n=n,k0=k0,k1=k1))
    }
  }else{
    # edit new version
    SUM = sum(normdn)
    sqrho = sqrt(rho)
    sq_rho = sqrt(1-rho)
    sqrhoNODES = sqrho*NODES
    result = 0
    if(onesided==T){
      qu = qnorm(1-u)
      qum = qnorm(1-um)
      if(length(qu)==1){
        D = t(as.matrix(sapply(1:nNODES, function(x)pnorm((qu-sqrhoNODES[x])/sq_rho,lower.tail = F))))
      }else{
        D = sapply(1:nNODES, function(x)pnorm((qu-sqrhoNODES[x])/sq_rho,lower.tail = F))
      }
      
      DM = sapply(1:nNODES, function(x)pnorm((qum-sqrhoNODES[x])/sq_rho,lower.tail = F))
      result = sum(normdn*sapply(1:nNODES, function(x)UnifCross(D[,x],DM[x],t=t,n=n,k0=k0,k1=k1)))/SUM
    }else{
      qu = qnorm(1-u/2)
      qum = qnorm(1-um/2)
      if(length(qu)==1){
        D = t(as.matrix(sapply(1:nNODES, function(x)pnorm((qu-sqrhoNODES[x])/sq_rho,lower.tail = F) +
                                 pnorm((-qu-sqrhoNODES[x])/sq_rho))))
      }else{
        D = sapply(1:nNODES, function(x)pnorm((qu-sqrhoNODES[x])/sq_rho,lower.tail = F) +
                     pnorm((-qu-sqrhoNODES[x])/sq_rho))
      }
      DM = sapply(1:nNODES, function(x)pnorm((qum-sqrhoNODES[x])/sq_rho,lower.tail = F) +
                    pnorm((-qum-sqrhoNODES[x])/sq_rho))
      result = sum(normdn*sapply(1:nNODES, function(x)UnifCross(D[,x],DM[x],t=t,n=n,k0=k0,k1=k1)))/SUM
    }
    return(result)
  }
}

#Crossing probablity of uniform order statistics
UnifCross <- function(u, um, t, n, k0, k1){
  #um = boundary_end
  #u = boundary
  m = floor(k1)
  #pp = exp(lfactorial(n) - lfactorial(n-(0:(t-1))))*pbeta(um,M+1,n-m+1,lower.tail=F)
  pp = rep(NA, t)
  pp[1] = pbeta(um,m,n-m+1,lower.tail=F)
  pp[2:t] = exp(lfactorial(n) - lfactorial(n-(k0:(k0+t-2))))*pbeta(um,m-(k0:(k0+t-2)),n-m+1,lower.tail=F)
  a = rep(1,t)
  a[2] = -dpois(k0,u[k0])*exp(u[k0])
  if(t>2){
    for (i in 2:(t-1)){
      d = dpois(c(1:(i-1), i+k0-1),u[i+k0-1])*exp(u[i+k0-1])
      # a[i+1]=-a[i:1]%*%d
      # edit new version
      a[i+1]=-sum(a[i:1]*d)
    }
    #p = 1-pp[1:t]%*%a[1:t]
    #if(p>=pbeta(u[2], 2, n-1)-u[1]*n*(1-pbeta(u[2],1,n-1))&&p<=1){
    #return(pp[1:(t-k0+1)]%*%a[1:(t-k0+1)])
    #return(drop(pp[1:t]%*%a[1:t]))
    # edit new version
    finalp = sum(pp*a)
    #}else{
    #warning("q is too small. The algorithm fails due to loss of significant digits", call. = FALSE)
    #return(0)
    #}
  }else if(t==2){
    finalp = sum(pp*a)
  }else{
    p = pbeta(u[k0], k0, n-k0+1)
    finalp = 1-drop(p)
  }
  if(finalp>1|finalp<0){
    p = pbeta(u[k0], k0, n-k0+1)
    finalp = 1-drop(p)
  }
  return(finalp)
}

qphi.rho <- function (p, n, rho,NODES, k0, k1,s,t=28,onesided=F, err_thr=1e-4) {
  lower = -1
  upper = 10
  maxnumIter = 100
  numIter = 1
  ddn = dnorm(NODES)
  p_cal_lower = pphi.rho(lower, n=n, rho=rho,NODES=NODES, k0=k0, k1=k1,s=s,t=t,onesided=onesided, normdn=ddn)
  p_cal_upper = pphi.rho(upper, n=n, rho=rho,NODES=NODES, k0=k0, k1=k1,s=s,t=t,onesided=onesided, normdn=ddn)
  while ((p_cal_lower - p) * (p_cal_upper - p) > 0 && numIter <
         maxnumIter) {
    if (p_cal_lower > p) {
      lower = lower * 1.1
      p_cal_lower = pphi.rho(lower, n=n, rho=rho,NODES=NODES, k0=k0, k1=k1,s=s,t=t,onesided=onesided, normdn=ddn)
    }
    if (p_cal_upper < p) {
      upper = upper * 1.1
      p_cal_upper = pphi.rho(upper, n=n, rho=rho,NODES=NODES, k0=k0, k1=k1,s=s,t=t,onesided=onesided, normdn=ddn)
    }
    numIter = numIter + 1
  }
  q = mean(c(lower, upper))
  p_cal = pphi.rho(q, n=n, rho=rho,NODES=NODES, k0=k0, k1=k1,s=s,t=t,onesided=onesided, normdn=ddn)
  error = (p_cal - p)/(1-p)
  numIter = 1
  while (abs(error) > err_thr && upper-lower>1e-4 && numIter < maxnumIter) {
    if (error > 0) {
      upper = q
      q = mean(c(lower, upper))
      p_cal = pphi.rho(q, n=n, rho=rho,NODES=NODES, k0=k0, k1=k1,s=s,t=t,onesided=onesided, normdn=ddn)
    }
    else {
      lower = q
      q = mean(c(lower, upper))
      p_cal = pphi.rho(q, n=n, rho=rho,NODES=NODES, k0=k0, k1=k1,s=s,t=t,onesided=onesided, normdn=ddn)
    }
    error = (p_cal - p)/(1-p)
    numIter = numIter + 1
  }
  return(q)
}



########################################################################
# calculate the right-tail probability of omnibus phi-divergence statistics under equal correlation matrix.
# threshold - threshold.
# n - number of p-values.
# rho - common correlation coefficient.
# K0 - vector of search range starts (from the k0th smallest p-value).
# K1 - vector of search range ends (at the k1th smallest p-value).
# S - vector of the phi-divergence test parameters.
# NODES - the numerical integration domain (under the standard normal pdf).
# t - numerical truncation parameter.
# onesided - TRUE if the input p-values are one-sided.
# pphi.rho.omni(0.05, n=10, rho=0.3, K0=rep(1,2), K1=rep(5,2), S=c(1,2))

pphi.rho.omni <- function(threshold, n, rho, K0, K1, S, NODES=seq(-8,8,length.out=50), t=30, onesided=F){
  t = min(t,max(K1)-min(K0)+1)
  M = floor(K1)
  Q = mapply(function(x,y,z)qphi.rho(p=1-threshold, n=n, rho=rho, NODES=NODES, k0=x, k1=y, s=z, t=t, onesided=onesided),
             K0,K1,S)
  
  id = which(Q=="Bisection fails, choose another lower or upper bound to try again")
  if(length(id)>0&length(id)<length(S)){
    Q = as.numeric(Q[-id])
    K0 = K0[-id]
    K1 = K1[-id]
    S = S[-id]
  }
  boundary =matrix(0, max(K1)-min(K0)+1, length(S))
  for (i in 1:length(S)) {
    boundary[K0[i]:K1[i],i] = mapply(function(x)phi.f.inv.i(Q[i],x,n,S[i]), K0[i]:K1[i])
  }
  u = apply(boundary,1,max)
  um = max(boundary[max(K1),])
  dnormNODES = dnorm(NODES)
  dnormz = dnormNODES/sum(dnormNODES)
  sqrhoz = sqrt(rho)*NODES
  sqrho_1 = sqrt(1-rho)
  result = 0
  if(onesided==T){
    qnormu = qnorm(1-u)
    qnormum = qnorm(1-um)
    
    for(i in 1:length(NODES)){
      d = pnorm((qnormu-sqrhoz[i])/sqrho_1,lower.tail = F)
      dm = pnorm((qnormum-sqrhoz[i])/sqrho_1,lower.tail = F)
      result = result + dnormz[i]*UnifCross(d,dm,t=t,n=n,k0=min(K0),k1=max(K1))
    }
  }else{
    qnormu = qnorm(1-u/2)
    qnormum = qnorm(1-um/2)
    for(i in 1:length(NODES)){
      d = pnorm((qnormu-sqrhoz[i])/sqrho_1,lower.tail = F) +
        pnorm((-qnormu-sqrhoz[i])/sqrho_1)
      dm = pnorm((qnormum-sqrhoz[i])/sqrho_1,lower.tail = F) +
        pnorm((-qnormum-sqrhoz[i])/sqrho_1)
      result = result + dnormz[i]*UnifCross(d,dm,t=t,n=n,k0=min(K0),k1=max(K1))
    }
  }
  return(1-result)
}
