#' @importFrom stats uniroot dbinom pgamma
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
