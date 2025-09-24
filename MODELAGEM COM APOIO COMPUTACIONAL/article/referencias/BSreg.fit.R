########################################################
# "BS quantile regression no-spatial" ##################
# Authors: L. Sanchez, V. Leiva ########################
# Date: 28-Julio-2018                               ####
########################################################
# File: Fit for the BS regression model                #
########################################################

# Required packages

library(MASS)
library(maxLik)
library(VGAM)

bsreg.fit = function(x, 
                     y,
                     link = "log") {
  # x is the design matrix with first column containing only one's
  # y is the vector of values of the response variable  
  
  n           = NROW(x)
  p           = NCOL(x)
  
  linkstr = link
  linkobj = make.link(linkstr)
  linkfun = linkobj$linkfun
  linkinv = linkobj$linkinv
  Q.eta  = linkobj$mu.eta
  
  ystar = linkfun(y)
  
  beta  = ginv(t(x)%*%x)%*%t(x)%*%ystar  # initial values for estimating beta's
  
  # #Usar este código comentado si hay errores en la estimación
  # if(linkstr=="identity"){
  #   beta = abs(ginv(t(x)%*%x)%*%t(x)%*%ystar)
  # }
  
  xbar  = mean(y)
  vart  = (n / (n - 1)) * var(y)
  r     = vart/xbar^2
  
  alphai = sqrt((2*r-2+2*sqrt(1+3*r))/(5-r)) # This expression is obtained by to express alpha
  # in terms of the variance and mean of Y ~ BS
  if( alphai=='NaN')
  {
    s1=mean(y)
    r1=1/mean(1/y)
    alphai=sqrt(2*sqrt(s1/r1)-2) # Una opción para valor de parámetro inicial
  }
  
  start = list(beta = beta, alpha = alphai)
  if (is.list(start)) 
    start = do.call("c", start)
  
  fr    = function(vp)
  {                                         
    betab = vp[1:p]
    eta   = as.vector(x%*%betab)
    Q     = linkinv(eta)
    alphab = vp[-(1:p)] 
    q    = 0.5 
    zq    = qnorm(q, mean = 0, sd = 1)
    gma_alphab = alphab * zq + sqrt(alphab^2 * zq^2+4)
    vt    = y
    sum(-0.5*log(8*pi*vt)-log(alphab)-log(gma_alphab)-0.5*log(Q)+log(gma_alphab^2/2+2*Q/vt)-
          (2*Q/(alphab^2*gma_alphab^2*vt))*(vt*gma_alphab^2/(4*Q)-1)^2)
  } 
  grr = function(vp)
  {
    betab = vp[1:p]
    eta   = as.vector(x%*%betab)
    Q    = linkinv(eta)  
    alphab = vp[-(1:p)]
    q     = 0.5 
    zq    = qnorm(q, mean = 0, sd = 1)
    gma_alphab = alphab * zq + sqrt(alphab^2 * zq^2+4)
    gma_alphabp= zq+zq^2*alphab*(1/sqrt(alphab^2*zq^2+4))
    gma_alphabpp = 4*zq^2*1/sqrt(alphab^2*zq^2+4)^3
    vt    = y
    z     =  -0.5*(1/Q)-2*(1/(alphab^2*gma_alphab^2*vt))+gma_alphab^2*vt*(1/(8*alphab^2*Q^2))+
      4*(1/(vt*gma_alphab^2+4*Q))
    b     =   -(gma_alphab+alphab*gma_alphabp)*(1/(alphab*gma_alphab))+
      2*vt*gma_alphab*gma_alphabp*(1/(vt*gma_alphab^2+4*Q))-
      (gma_alphab*gma_alphabp*alphab-gma_alphab^2)*vt*(1/(4*Q*alphab^3))-2*(1/(alphab^3))+
      4*Q*(gma_alphab+alphab*gma_alphabp)*(1/(alphab^3*gma_alphab^3*vt))
    
    rval  = cbind(t(t(x)%*%diag(Q.eta(eta))%*%z), sum(b))
    
    return(rval)
  }
  
  A = matrix(c(rep(0,p),1),1,p+1)
  B = 0
  
  # # Usar este código si hay errores en la estimación.
  # if(linkstr=="identity"){
  #   A = diag(x=1, nrow=p+1, ncol=p+1)
  #   B = as.vector(rep(0,p+1))
  # }
  
  opt = maxBFGS(fn = fr, grad = grr, start = start, constraints=list(ineqA=A, ineqB=B))
  
  if (opt$code > 0)  warning("optimization failed to converge")

  log.lik.est = opt$maximum
  estimates   = opt$estimate

  #Criterios para selección de modelos
  AIC   = - 2 * log.lik.est + 2 * (p+1)
  AICc  = AIC + (2 * (p+1) * ((p+1) + 1)) / (n - (p+1) - 1)
  BIC   = - 2 * log.lik.est + log(n) * (p+1)

  beta  = as.vector(estimates[1:p])
  eta   = as.vector(x%*%beta)
  Q     = linkinv(eta)
  alpha = estimates[-(1:p)]
  q     = 0.5 
  zq    = qnorm(q, mean = 0, sd = 1)

  aux     = matrix(1, ncol = 1L, nrow = n)

  gma_alpha = alpha * zq + sqrt(alpha^2 * zq^2+4) 
  gma_alphap= zq+alpha*zq^2*(1/sqrt(alpha^2*zq^2+4))
  gma_alphapp = 4*zq^2*(1/sqrt(alpha^2*zq^2+4)^3)

  par_alpha     = alpha
  par_beta      = 4*Q/gma_alpha^2

  Acal = ((2*gma_alphap+alpha*gma_alphapp)*(alpha*gma_alpha)-(gma_alpha+alpha*gma_alphap)^2)/
          (alpha^2*gma_alpha^2)
  Bcal = 8*Q*(gma_alphap^2+gma_alpha*gma_alphapp)
  Ccal = (1/(4*Q*alpha^4))*(alpha^2*gma_alphap^2+alpha^2*gma_alpha*gma_alphapp-alpha*gma_alpha*gma_alphap-
          3*gma_alpha*gma_alphap*alpha+3*gma_alpha^3)
  Dcal = 4*Q*((2*gma_alphap+alpha*gma_alphapp)*alpha*gma_alpha-3*(gma_alpha+alpha*gma_alphap)^2)/
          (alpha^4*gma_alpha^4)

  if(link == "log")
    {
      a  = Q
      h1 = Q^2
      h2 = -1/(Q^2*(log(Q)^3))
    }
  
  if(link == "identity")
    {
      a  = rep(1,n)
      h1 = 1
      h2 = 0
    }
  
  if(link == "sqrt")
    {
      a  = 2*sqrt(Q)
      h1 = 4*Q
      h2 = -1/(4*Q^3)
    }
    
    integ.f2 = c()
    integ.f3 = c()
    integ.f4 = c()
    integ.f5 = c()
    
  
  f2 <- function(u,j)
    {
      return((u/(u*gma_alpha^2+4*Q[j]))^2*dbisa(u, par_beta[j], par_alpha))
    }
    
  f3 <- function(u,j)
    {
      return((u/(u*gma_alpha^2+4*Q[j])^2)*dbisa(u, par_beta[j], par_alpha))
    }
    
  f4 <- function(u,j)
    {
      return((1/(u*gma_alpha^2+4*Q[j])^2)*dbisa(u, par_beta[j], par_alpha))
    }
    
  f5 <- function(u,j)
    {
      return((1/(u*gma_alpha^2+4*Q[j]))*dbisa(u, par_beta[j], par_alpha))
    }
    
  for(i in 1:n)
    {
      integ.f2[i] = integrate(f2, lower=0, upper=1e100, j=i)$value
      integ.f3[i] = integrate(f3, lower=0, upper=1e100, j=i)$value
      integ.f4[i] = integrate(f4, lower=0, upper=1e100, j=i)$value
      integ.f5[i] = integrate(f5, lower=0, upper=1e100, j=i)$value
    }
    
    
  v = (-1/(2*Q^2)+16*integ.f4 + (1/(alpha^2*Q^2))*(1+alpha^2/2))*h1-
      (1/(2*Q)+(1/(2*alpha^2*Q))*(1+alpha^2/2)-4*integ.f5)*h2
    
  s = 8*gma_alpha*gma_alphap*integ.f3 - 
      ((gma_alpha * gma_alphap * alpha- gma_alpha^2)/(alpha^3*gma_alpha^2*Q^2))*(1+alpha^2/2)-
      ((gma_alpha +alpha* gma_alphap)/(alpha^3*gma_alpha *Q))*(1+alpha^2/2)
    
  u = Acal-Bcal*integ.f3 -2*(gma_alpha^3*gma_alphapp-gma_alpha^2*gma_alphap^2)^2*
      integ.f2 + Ccal*(4*Q/gma_alpha^2)*(1+alpha^2/2)-6/alpha^4-
      Dcal*(gma_alpha^2/(4*Q))*(1+alpha^2/2)
  
  kbb  = t(x)%*%diag(as.vector(v))%*%x
  kaa  = sum(diag(as.vector(u)))
  kba  = t(x)%*%diag(as.vector(a))%*%s
  
  Diag = function(A){
    diag.A = vector()
    for(t in 1:ncol(A)){
      diag.A[t]=A[t,t]
    }
    return(as.vector(diag.A))
  }
  
  fisher = cbind(rbind(kbb, t(kba)), rbind(kba, kaa))
  se    = sqrt(Diag(solve(fisher)))
  # se     = numeric()
  # for(t in 1:ncol(fisher)){
  #   se[t] = se1[t,t] 
  # }
  
  hess = as.matrix(opt$hessian)

  if(p == 1) {
    var.explic  = x
  }  else {
            var.explic = x[,-1]
          }

  zstatbeta  = beta / se[1:p]
  zstatalpha   = alpha / se[p:p+1] 
  pvalorbeta = 2 * pnorm(abs(zstatbeta), lower.tail = F)
  pvaloralpha  = 2 * pnorm(abs(zstatalpha), lower.tail = F)

  names(beta)  = colnames(x)
  
  rval = list(coefficients = list(beta = beta, alpha = alpha), 
              fitted.values.Q = structure(Q, .Names = names(y)), 
              #optim = opt, start = start, 
              n = n, p = p, q = q, eta = eta, 
              X = x, y=y, #var.explic = var.explic,
              Hessian = hess,
              matrix.expected.Fisher = fisher, loglik = log.lik.est, 
              link = list(quantile = linkobj), 
              converged = opt$message, 
              information.criterions = list(aic = AIC,bic = BIC,aicc=AICc), 
              se = se, zstat = list(beta = zstatbeta, alpha = zstatalpha),
              pvalor = list(beta = pvalorbeta, alpha = pvaloralpha))


  return(rval)


}



#kbb = crossprod(v*x, x)
#kaa = crossprod(u*aux, aux)
#kba = crossprod(s*x, aux)


# for(i in 1:n)
# {
#   f2 <- function(u)
#   {
#     return((u/(u*gma_alpha^2+4*Q[i]))^2*dbisa(u, par_beta[i], par_alpha))
#   }
#   
#   f3 <- function(u)
#   {
#     return((u/(u*gma_alpha^2+4*Q[i])^2)*dbisa(u, par_beta[i], par_alpha))
#   }
#   
#   f4 <- function(u)
#   {
#     return((1/(u*gma_alpha^2+4*Q[i])^2)*dbisa(u, par_beta[i], par_alpha))
#   }
#   
#   f5 <- function(u)
#   {
#     return((1/(u*gma_alpha^2+4*Q[i]))*dbisa(u, par_beta[i], par_alpha))
#   }
#   
#   inte.f2[i] = integrate(f2,0,Inf)$value
#   inte.f3[i] = integrate(f3,0,Inf)$value
#   inte.f4[i] = integrate(f4,0,Inf)$value
#   inte.f5[i] = integrate(f5,0,Inf)$value
# }

# 
# v = (-1/(2*Q^2)+16*integrate(f4,0,Inf)$value + (1/(alpha^2*Q^2))*(1+alpha^2/2))*h1-
#     (1/(2*Q)+(1/(2*alpha^2*Q))*(1+alpha^2/2)-4*integrate(f5,0,Inf)$value)*h2
# 
# s = 8*gma_alpha*gma_alphap*integrate(f3,0,Inf)$value - 
#     ((gma_alpha * gma_alphap * alpha- gma_alpha^2)/(alpha^3*gma_alpha^2*Q^2))*(1+alpha^2/2)-
#     ((gma_alpha +alpha* gma_alphap)/(alpha^3*gma_alpha *Q))*(1+alpha^2/2)
# 
# u = Acal-Bcal*integrate(f3, 0, Inf)$value -2*(gma_alpha^3*gma_alphapp-gma_alpha^2*gma_alphap^2)^2*
#     integrate(f2, 0, Inf)$value + Ccal*(4*Q/gma_alpha^2)*(1+alpha^2/2)-6/alpha^4-
#     Dcal*(gma_alpha^2/(4*Q))*(1+alpha^2/2)



# else 
# {
#   hess = as.matrix(-optimHess(opt$estimate, fr, grr))
#   se   = sqrt(diag(solve(hess)))
# }
# 
# hess = as.matrix(optimHess(opt$estimate, fr, grr))



# rval = list(coefficients = list(quantile = beta, alpha = alpha), 
#           fitted.values = structure(Q, .Names = names(y)), 
#           optim = opt, start = start, n = n, p = p, q = q, eta = eta, 
#           Hessian = hess, X = x, var.explic = var.explic, 
#           logical.hessian = Fisher, loglik = log.lik, 
#           link = list(quantile = linkobj), 
#           converged = opt$convergence, 
#           information.criterions = list(aic = AIC,bic = BIC,aicc=AICc), 
#           se = se, zstat = list(beta = zstatbeta, alpha = zstatalpha),
#           pvalor = list(beta = pvalorbeta, alpha = pvaloralpha))