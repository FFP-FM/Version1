



library(mvtnorm)
source("FP_sup.R")


#load data 
dat <- data.load(pheno="pheno.csv",marker="genotype.csv",time=1:31)

#Null hypothesis
H0 <- mle_curve(pheno=dat$phenotype,times=dat$time)

#Alternative hypothesis
ret <- mle_H1(dat,times=dat$time)

head(ret)

#P-value of all markers

P <- ret[,2]

P
