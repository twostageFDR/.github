####
library(EnvStats)
library(BiocManager)
library(swfdr)
#install.packages("devtools")
library(devtools)
install_github("fxia22/RadaFDR", force = TRUE)
library(RadaFDR)

adafdr_test <- function (p_input, x_input, K, alpha, n_itr, qt_norm, verbose, 
                         output_folder, random_state, single_core, fast_mode, n_full, 
                         covariate_type, h){
  if (missing(K)) {
    K = 5L
  }
  if (missing(alpha)) {
    alpha = 0.1
  }
  if (missing(n_itr)) {
    n_itr = 1500L
  }
  if (missing(qt_norm)) {
    qt_norm = TRUE
  }
  if (missing(verbose)) {
    verbose = FALSE
  }
  if (missing(random_state)) {
    random_state = 0L
  }
  if (missing(single_core)) {
    single_core = TRUE
  }
  if (missing(fast_mode)) {
    fast_mode = TRUE
  }
  if (missing(n_full)) {
    n_full = NULL
  }
  if (missing(covariate_type)) {
    covariate_type = NULL
  }
  if (missing(h)) {
    h = NULL
  }
  if (missing(output_folder)) {
    output_folder = NULL
  }
  md <- import("adafdr.method")
  md$adafdr_test(p_input, x_input, fast_mode = fast_mode, K = K, 
                 alpha = alpha, n_itr = n_itr, qt_norm = qt_norm, verbose = verbose, 
                 random_state = random_state, single_core = single_core, 
                n_full = n_full, covariate_type = covariate_type, 
                 h = h, output_folder = output_folder)
}

pi_0 = function(lambda, p) sum(p > lambda)/(1-lambda)/length(p)

fdr = function(mat){
  s = apply(mat, 1, function(x) x[2]/max(sum(x[1:2]),1))
  m = mean(s)
  std = sd(s)
  
  return(c(m, std))
}
power = function(mat){
  s = apply(mat, 1, function(x) x[1]/max(sum(x[c(1,3)]),1))
  m = mean(s)
  std = sd(s)
  return(c(m, std))
}

fdr_and_power = function(ar){
  l = ncol(ar)/3
  n = dim(ar)[3]
  fdrmat  = matrix(0.0, l, 4)
  powermat = matrix(0.0, l, 4)
  for(i in 1:l){
    lseq = seq(3*(i-1)+1, by = 1, length=3)
    mat1 <- matrix(ar[1,lseq,], nr = n, nc = ncol(ar), byrow = T)
    mat2 <- matrix(ar[2,lseq,], nr = n, nc = ncol(ar), byrow = T)
    fdrmat[i, 1:2] = fdr(mat1)
    fdrmat[i, 3:4] = fdr(mat2)
    powermat[i, 1:2] = power(mat1)
    powermat[i, 3:4] = power(mat2)
  }
  return(list(fdr=fdrmat, power = powermat))
}

simgen = function(seed, n, cop, p0 = 0.95, mm=3){
  
  set.seed(seed)
  
  n0 = round(n * p0)
  n1 = n - n0
  
  lfc0 = rnorm(n0)
  temp = pnorm(lfc0)
  pv01 = 2*pmin(1-temp, temp)
  pv02 = BiCopCondSim(n0, pv01, 1, obj = cop)
  lfc_sd0 = qgamma(1-pv02, 3,4)
  
  
  lfc1 = rnorm(n1,mm)*sample(c(1,-1),n1,T)
  temp = pnormMix(lfc1, -mm, 1, mm, 1, 0.5)
  pv11 = 2*pmin(1-temp, temp)
  pv12 = BiCopCondSim(n1, pv11, 1, obj = cop)
  lfc_sd1 = qgamma(1-pv12, 3,4)
  
  
  lfc = c(lfc0, lfc1)
  lfc_sd = c(lfc_sd0, lfc_sd1)
  
  temp = pnorm(lfc)
  p1 = 2*pmin(temp, 1-temp)
  p2 = pgamma(lfc_sd,3,4)
  
  idx = c(rep(F,n0),rep(T,n1))
  
  return(cbind(lfc_sd, lfc, idx))
}

test_locfdr = function(q, idx){
  rej = q<0.05
  rej2 = q<0.10
  
  fdrset1 = c(sum(rej & idx),sum(rej & !idx), sum(!rej & idx))
  fdrset2 = c(sum(rej2 & idx),sum(rej2 & !idx), sum(!rej2 & idx))
  
  return(rbind(fdrset1, fdrset2))
}

test_Storey = function(p2, idx){
  pval_storey = p2
  
  p0hat = pi_0(0.5, pval_storey)
  rej_storey = (p0hat*n*pval_storey/rank(pval_storey)<0.05)
  rej_storey2 = (p0hat*n*pval_storey/rank(pval_storey)<0.1)
  
  fdrset1 = c(sum(rej_storey & idx),sum(rej_storey & !idx), sum(!rej_storey & idx))
  fdrset2 = c(sum(rej_storey2 & idx),sum(rej_storey2 & !idx), sum(!rej_storey2 & idx))
  
  return(rbind(fdrset1, fdrset2))
}

test_H = function(p1,p2,idx,copest){
  gammaseq = seq(0.9,1,0.001)
  nrej = sapply(gammaseq, function(x){
    pv = BiCopCDF(rep(x, n), p2, obj = copest)
    pval_H = ifelse(p1>x, p1, pv)
    p0hat = pi_0(0.5, pval_H)
    false_H = rank(pval_H)
    rej_H = (p0hat*n*pval_H/false_H<0.05)
    return(sum(rej_H))
  })
  
  gamma = gammaseq[which.max(nrej)]
  pv = BiCopCDF(rep(gamma, n), p2, obj = copest)
  pval_H = ifelse(p1>gamma, p1, pv)
  p0hat = min(pi_01(0.5, pval_H),1)
  false_H = rank(pval_H)
  rej_H = (p0hat*n*pval_H/false_H<0.05)
  rej_H2 = (p0hat*n*pval_H/false_H<0.1)
  
  fdrset1 = c(sum(rej_H & idx),sum(rej_H & !idx), sum(!rej_H & idx))
  fdrset2 = c(sum(rej_H2 & idx),sum(rej_H2 & !idx), sum(!rej_H2 & idx))
  
  return(list(fdrset = rbind(fdrset1, fdrset2), gamma = gamma))
}

test_S = function(p1,p2,idx, copest){
  
  pvalS = BiCopHfunc1(p1, p2, obj = copest)
  p0hat = min(pi_01(0.01, pvalS),1)
  rej_S = (p0hat*n*pvalS/rank(pvalS)<0.05)
  rej_S2 = (p0hat*n*pvalS/rank(pvalS)<0.1)
  
  fdrset1 = c(sum(rej_S & idx),sum(rej_S & !idx), sum(!rej_S & idx))
  fdrset2 = c(sum(rej_S2 & idx),sum(rej_S2 & !idx), sum(!rej_S2 & idx))
  
  return(rbind(fdrset1, fdrset2))
}

test_BL = function(p1, p2, idx){
  qvalues <- lm_qvalue(p2, X=p1)
  newp = qvalues$qvalues
  test = p.adjust(newp,"BH")
  
  rej1 = test < 0.05
  rej2 = test < 0.10
  
  fdrset1 = c(sum(rej1 & idx), sum(rej1&!idx), sum(!rej1&idx))
  fdrset2 = c(sum(rej2 & idx), sum(rej2&!idx), sum(!rej2&idx))
  
  return(rbind(fdrset1, fdrset2))
}

test_IHW = function(p1, p2, idx){
  dat = data.frame(p2 = p2, p1 = p1)
  res_ihw <- ihw(p2 ~ p1, alpha = 0.05)
  pval_ihw = adj_pvalues(res_ihw)
  
  rej1 = pval_ihw < 0.05
  rej2 = pval_ihw < 0.10
  fdrset1 = c(sum(rej1 & idx), sum(rej1&!idx), sum(!rej1&idx))
  fdrset2 = c(sum(rej2 & idx), sum(rej2&!idx), sum(!rej2&idx))
  return(rbind(fdrset1, fdrset2))
}

test_AdaFDR = function(p1, p2, idx){
  p1 = as.matrix(p1)
  rej1 <- adafdr_test(p2, p1, alpha = 0.05)$decision
  rej2 <- adafdr_test(p2, p1, alpha = 0.10)$decision

  fdrset1 = c(sum(rej1 & idx), sum(rej1&!idx), sum(!rej1&idx))
  fdrset2 = c(sum(rej2 & idx), sum(rej2&!idx), sum(!rej2&idx))
  return(rbind(fdrset1, fdrset2))
}

test_all = function(simdat, copnum = 23){
  
  x = simdat[,1]
  y = simdat[,2]
  idx = (simdat[,3]==1)
  
  n = length(x)
  gg = locfdr(y, plot = 0)
  g = gg$fp0[3,]
  mu = g[1]
  s = g[2]
  temp = pnorm(y, mu, s)
  p2 = 2*pmin(temp, 1-temp)
  p1 = rank(x)/n
  
  corr = wdm(p1, p2, 'kendall')
  
  if(copnum == 23){
    if(corr>0){
      copest = BiCop(3, tau = corr)
    }else{
      copest = BiCop(23, tau = corr)
    }
  }else{
    copest = BiCopEst(p1, p2, copnum)
  }
  
  res_locfdr = test_locfdr(gg$fdr, idx)
  res_storey = test_Storey(p2, idx)
  res_Hlst = test_H(p1, p2, idx, copest)
  res_H = res_Hlst$fdrset
  gamma = res_Hlst$gamma
  res_S = test_S(p1, p2, idx, copest)
  res_BL = test_BL(p1, p2, idx)
  res_ihw = test_IHW(p1, p2, idx)
  res_AdaFDR = test_AdaFDR(p1,p2,idx)
  
  res = cbind(res_locfdr, res_storey, res_H, res_S, res_BL, res_ihw, res_AdaFDR)
  
  #res1 = c(res_locfdr[1,], res_storey[1,], res_H[1,], res_S[1,], res_BL[1,], res_ihw[1,])
  #res2 = c(res_locfdr[2,], res_storey[2,], res_H[2,], res_S[2,], res_BL[2,], res_ihw[2,])
  
  return(res)
}

