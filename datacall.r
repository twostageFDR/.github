####data call####
library(locfdr)
source("ramos.R")

myfile<-read.csv("set4delta.csv")
n = nrow(myfile)
myfile$snrKO<-1/myfile$cvKO
myfile$snrWT<-1/myfile$cvWT

x11<-snrKO<-myfile$snrKO
x12<-snrWT<-myfile$snrWT
x1 = pmax(x11,x12)
maxsnr = x1

lfc<- myfile$realdata

ko = myfile[,paste0('KO',1:3)]
wt = myfile[,paste0('WT',1:3)]
snrko = apply(ko, 1, function(x) mean(x)/sd(x))
snrwt = apply(wt, 1, function(x) mean(x)/sd(x))
set4 = data.frame(wt,ko)


m11 = 0
s11 = sqrt(0.004)
m22 = -0.002
s22 = sqrt(0.042)
eta = 0.615
cutoff = 0.4

pnormmix = function(x, m, s, eta) eta*pnorm(x, m[1], s[1]) + (1-eta)*pnorm(x, m[2], s[2]) 
dnormmix = function(x, m, s, eta) eta*dnorm(x, m[1], s[1]) + (1-eta)*dnorm(x, m[2], s[2]) 

y2 = eta*pnorm(lfc, m11, s11)+(1-eta)*pnorm(lfc, m22,s22)
p2 = 2*pmin(y2, 1-y2)

Fc = mean(abs(lfc)<cutoff)
F0c = pnormmix(cutoff, c(m11,m22), c(s11,s22), c(eta, 1-eta)) - pnormmix(-cutoff, c(m11,m22), c(s11,s22), c(eta, 1-eta))
p0hat = min(1, Fc/F0c)

fmat = gpar$f

f = sapply(lfc, function(x){
  loc = which.min(fmat[,1]<x)
  if(length(loc)==0) loc <- nrow(fmat)
  return(fmat[loc,2])
})

f0 = dnormmix(lfc, c(m11,m22), c(s11,s22), c(eta, 1-eta))

rej_locfdr = (p0hat*f0/f < 0.05)
rej_locfdr2 = (p0hat*f0/f < 0.1)


n = nrow(myfile)
pval_null = eta*pnorm(lfc, m11, s11)+(1-eta)*pnorm(lfc, m22, s22)
pval2 = 2*pmin(pval_null, 1-pval_null)

small = 0.0000000001

com1 = rep(1:3, each = 9)
com2 = rep(rep(1:3, each = 3), 3)
com3 = rep(1:3, 9)

comset = cbind(com1, com2, com3)

comc1 = rep(1:27,each = 27)
comc2 = rep(1:27, 27)

combin = cbind(comset[comc1,], comset[comc2,])


kom<-as.matrix(ko)
wtm<-as.matrix(wt)
mm = cbind(kom, wtm)

group_boot = function(d){
  x = d[1:3]
  y = d[4:6]
  
  bootx = apply(combin[,1:3], 1, function(a) mean(x[a]))
  booty = apply(combin[,4:6], 1, function(a) mean(y[a]))
  
  dseq = log(bootx) - log(booty)
  
  return(sd(dseq))
}

lfc_sd <- apply(mm, 1, group_boot)
pval1 = rank(lfc_sd)/n


hxt = c("YHR096C", "YJL214W", "YEL069C", "YDL245C", "YJR158W", "YDR536W")
hxtloc = sapply(hxt, function(x) grep(x, myfile$label4))
ishxt = rep(FALSE, n)
ishxt[hxtloc] = TRUE


### DESeq2
library(DESeq2)
set4r = round(set4)

colData <- data.frame(
  condition = factor(rep(c("Control", "Treatment"), each = 3))
)
rownames(colData) <- colnames(set4r)

dds <- DESeqDataSetFromMatrix(countData = set4r,
                              colData = colData,
                              design = ~ condition)

ddsr = DESeq(dds)
# Conducting a differential expression analysis
deRes <- as.data.frame(results(ddsr))
pval_deseq = deRes$padj
lfc_sd_deseq = deRes$lfcSE
pval1_deseq = rank(lfc_sd_deseq)/n
