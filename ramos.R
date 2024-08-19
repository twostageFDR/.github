library(locfdr)
library(mixtools)
library(EnvStats)
library(truncnorm)
library(splines)


em_truncated_mixture <- function(data, c, max.iter=1000, tol=1e-1){
  # Initial parameters
  mu <- c(mean(data) * 0.85, mean(data) * 1.15)
  sigma <- c(sd(data)*1/2, sd(data)*2)
  pi <- 0.5
  
  loglik_old <- -Inf
  
  for(iter in 1:max.iter){
    
    # E-Step
    tau1 <- dtruncnorm(data, a=-c, b=c, mean=mu[1], sd=sigma[1]) * pi
    tau2 <- dtruncnorm(data, a=-c, b=c, mean=mu[2], sd=sigma[2]) * (1-pi)
    
    tau1 <- tau1 / (tau1 + tau2)
    tau2 <- 1 - tau1
    
    # M-Step
    mu[1] <- sum(tau1 * data) / sum(tau1)
    mu[2] <- sum(tau2 * data) / sum(tau2)
    sigma[1] <- sqrt(sum(tau1 * (data - mu[1])^2) / sum(tau1))
    sigma[2] <- sqrt(sum(tau2 * (data - mu[2])^2) / sum(tau2))
    pi <- mean(tau1)
    
    # Log-Likelihood
    loglik_new <- sum(log(pi * dtruncnorm(data, a=-c, b=c, mean=mu[1], sd=sigma[1]) +
                            (1-pi) * dtruncnorm(data, a=-c, b=c, mean=mu[2], sd=sigma[2])))
    
    # Check for convergence
    if(abs(loglik_new - loglik_old) < tol){
      break
    }
    
    loglik_old <- loglik_new
  }
  
  return(list(mu=mu, sigma=sigma, pi=pi))
}

find_locmax <- function(x) (which(diff(sign(diff(x))) == -2) + 1)

Lc=function(dat, trn = 0.99, df = 10 ){
  cutoff = quantile(abs(dat), trn)
  trndat = dat[abs(dat)<cutoff]
  

  find_f = locfdr(trndat, plot = 0, type = 0, bre = 1000)
  
  
  fmat = find_f$mat[,c('x','f')]
  
  gap = fmat[2,1]-fmat[1,1]
  tot = sum(gap*fmat[,2])
  
  fmat[,2] = fmat[,2]/tot
  
  f = sapply(dat, function(x){
    loc = which.min(fmat[,1]<x)
    if(length(loc)==0) loc <- nrow(fmat)
    return(fmat[loc,2])
  })
  
  
  fixeda = function(a, c){
    tdat_a = dat[abs(dat)<=c]
    em_a = em_truncated_mixture(tdat_a, c)
    mu_a = em_a$mu
    s_a = em_a$sigma
    p_a = em_a$pi
    
    f0 = dnormMix(dat, mu_a[1], s_a[1], mu_a[2], s_a[2], p_a)
    lik = log(f0/f)
    lca = sum(lik[abs(dat)>a & abs(dat)<c])
    return(lca)
  }
  
  
  
  
  emres <- em_truncated_mixture(trndat, cutoff)
  
  mu = emres$mu
  sigma = emres$sigma
  prop = emres$pi
  #find local maxima and if num(local max)==1, return local max 
  polyfit = function(cset, lca){
    lmfit = lm(lca~poly(cset,4))
    locmax = find_locmax(lmfit$fitted.values)
    nlocmax = length(locmax)
    gap = lmfit$fitted.values[2]-lmfit$fitted.values[1]
    res = ifelse(nlocmax==1, cset[locmax], F)
    return(res)
  }
  
  cfix = NA
  
  a = seq(sigma[1], 2*sigma[1], 0.25*sigma[1])
  ncset = 2.25*sigma[1]<abs(dat) & abs(dat)<1.5
  rng = seq(1, sum(ncset), 50)
  cset = sort(abs(dat)[ncset][rng])
  
  for(a in seq(sigma[1], 2*sigma[1], 0.25*sigma[1])){
    #print(a)
    lca <- sapply(cset, function(c) fixeda(a,c))
    
    check = polyfit(cset, lca)
    
    if(check){
      cfix = check
      break
    }
  }
  
  cfix <- ifelse(is.na(cfix), quantile(abs(dat), 0.8), cfix)
  
  tc = dat[abs(dat)<cfix]
  emres = em_truncated_mixture(tc, cfix)
  mu = emres$mu
  sigma = emres$sigma
  prop = emres$pi


  return(list(param = c(mu, sigma, prop, cfix), f = fmat))
}

gpar = Lc(lfc)

#m11 = (gpar$param[1])
#m22 = (gpar$param[2])
#s11 = (gpar$param[3])
#s22 = (gpar$param[4])
#eta = (gpar$param[5])
#pval2 = eta*pnorm(lfc, m11, s11)+(1-eta)*pnorm(lfc, m22, s22)
#cutoff = gpar$param[6]
