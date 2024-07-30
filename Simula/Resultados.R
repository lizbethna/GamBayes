### RESULTADOS

##################################################

setwd <- "~/Documents/Sabatico/Examples/Simula_f1increasing 231028/Resultados N20 re1 casoC/"
dir <- "~/Documents/Sabatico/Examples/Simula_f1increasing 231028/Resultados N20 re1 casoC/"

##################################################


PARAM1 <- read.csv(paste0(dir,"PARAM.lme.non.re1.csv"))
PARAM2 <- read.csv(paste0(dir,"PARAM.lme.incr.re1.csv"))
PARAM3 <- read.csv(paste0(dir,"PARAM.lme.add.non.re1.csv"))
PARAM4 <- read.csv(paste0(dir,"PARAM.lme.add.incr.re1.csv"))

nombres1 = c("i","b0","b1","b2", "invtau2","tau2","tau", 
             "invsig2[1]","invsig2[2]","sig2[1]","sig2[2]", "sigma[1]","sigma[2]")
nombres2 = c("i","b0","b1","b2", "invtau2","tau2","tau", 
             "invsig2[1]","invsig2[2]","sig2[1]","sig2[2]", "sigma[1]","sigma[2]")
nombres3 = c("i","b0","b1[1]","b1[2]","b1[3]","b1[4]","b1[5]","b1[6]", 
             "b2[1]","b2[2]","b2[3]","b2[4]","b2[5]","b2[6]", 
             "invtau2","tau2","tau", 
             "lambda[1]","lambda[2]","rho[1]","rho[2]",
             "invsig2[1]","invsig2[2]","sig2[1]","sig2[2]", "sigma[1]","sigma[2]")
nombres4 = c("i","b0","b1[1]","b1[2]","b1[3]","b1[4]","b1[5]","b1[6]", 
             "b2[1]","b2[2]","b2[3]","b2[4]","b2[5]","b2[6]", 
             "invtau2","tau2","tau", 
             "lambda[1]","lambda[2]","rho[1]","rho[2]",
             "invsig2[1]","invsig2[2]","sig2[1]","sig2[2]", "sigma[1]","sigma[2]")


ERROR1 <- read.csv(paste0(dir,"ERROR.lme.non.re1.csv"))
ERROR2 <- read.csv(paste0(dir,"ERROR.lme.incr.re1.csv"))
ERROR3 <- read.csv(paste0(dir,"ERROR.lme.add.non.re1.csv"))
ERROR4 <- read.csv(paste0(dir,"ERROR.lme.add.incr.re1.csv"))

ELPDLOO1 <- read.csv(paste0(dir,"ELPD.LOO.lme.non.re1.csv"))
ELPDLOO2 <- read.csv(paste0(dir,"ELPD.LOO.lme.incr.re1.csv"))
ELPDLOO3 <- read.csv(paste0(dir,"ELPD.LOO.lme.add.non.re1.csv"))
ELPDLOO4 <- read.csv(paste0(dir,"ELPD.LOO.lme.add.incr.re1.csv"))

ELPDWAIC1 <- read.csv(paste0(dir,"ELPD.WAIC.lme.non.re1.csv"))
ELPDWAIC2 <- read.csv(paste0(dir,"ELPD.WAIC.lme.incr.re1.csv"))
ELPDWAIC3 <- read.csv(paste0(dir,"ELPD.WAIC.lme.add.non.re1.csv"))
ELPDWAIC4 <- read.csv(paste0(dir,"ELPD.WAIC.lme.add.incr.re1.csv"))

##################################################

SIM = 100
k1 = 6; k2 = 6 
PARAM1 = PARAM2 = array(NA, dim=c(SIM,12))
PARAM3 = PARAM4 = array(NA, dim=c(SIM,12+2+k1+k2))
 
ERROR1 = ERROR2 = ERROR3 = ERROR4 =array(NA, dim=c(SIM,3*5))
ELPDLOO1 = ELPDLOO2 = ELPDLOO3 = ELPDLOO4 = array(NA, dim=c(SIM,6))
ELPDWAIC1 = ELPDWAIC2 = ELPDWAIC3 = ELPDWAIC4 = array(NA, dim=c(SIM,6))

sec = c("_a","_b","_c")
lin_ini = c( 1, 19,  36)
lin_fin = c(18, 35, 100)

for(j in 1:3){
PARAM1j <- read.csv(paste0(dir,"PARAM.lme.non.re1",sec[j],".csv"))
PARAM2j <- read.csv(paste0(dir,"PARAM.lme.incr.re1",sec[j],".csv"))
PARAM3j <- read.csv(paste0(dir,"PARAM.lme.add.non.re1",sec[j],".csv"))
PARAM4j <- read.csv(paste0(dir,"PARAM.lme.add.incr.re1",sec[j],".csv"))

ERROR1j <- read.csv(paste0(dir,"ERROR.lme.non.re1",sec[j],".csv"))
ERROR2j <- read.csv(paste0(dir,"ERROR.lme.incr.re1",sec[j],".csv"))
ERROR3j <- read.csv(paste0(dir,"ERROR.lme.add.non.re1",sec[j],".csv"))
ERROR4j <- read.csv(paste0(dir,"ERROR.lme.add.incr.re1",sec[j],".csv"))

ELPDLOO1j <- read.csv(paste0(dir,"ELPD.LOO.lme.non.re1",sec[j],".csv"))
ELPDLOO2j <- read.csv(paste0(dir,"ELPD.LOO.lme.incr.re1",sec[j],".csv"))
ELPDLOO3j <- read.csv(paste0(dir,"ELPD.LOO.lme.add.non.re1",sec[j],".csv"))
ELPDLOO4j <- read.csv(paste0(dir,"ELPD.LOO.lme.add.incr.re1",sec[j],".csv"))

ELPDWAIC1j <- read.csv(paste0(dir,"ELPD.WAIC.lme.non.re1",sec[j],".csv"))
ELPDWAIC2j <- read.csv(paste0(dir,"ELPD.WAIC.lme.incr.re1",sec[j],".csv"))
ELPDWAIC3j <- read.csv(paste0(dir,"ELPD.WAIC.lme.add.non.re1",sec[j],".csv"))
ELPDWAIC4j <- read.csv(paste0(dir,"ELPD.WAIC.lme.add.incr.re1",sec[j],".csv"))


idx = lin_ini[j]:lin_fin[j]


PARAM1[idx,] = as.matrix(PARAM1j[idx,-1])
PARAM2[idx,] = as.matrix(PARAM2j[idx,-1])
PARAM3[idx,] = as.matrix(PARAM3j[idx,-1])
PARAM4[idx,] = as.matrix(PARAM4j[idx,-1])

ERROR1[idx,] <- as.matrix(ERROR1j[idx,-1])
ERROR2[idx,] <- as.matrix(ERROR2j[idx,-1])
ERROR3[idx,] <- as.matrix(ERROR3j[idx,-1])
ERROR4[idx,] <- as.matrix(ERROR4j[idx,-1])

ELPDLOO1[idx,] <- as.matrix(ELPDLOO1j[idx,-1])
ELPDLOO2[idx,] <- as.matrix(ELPDLOO2j[idx,-1])
ELPDLOO3[idx,] <- as.matrix(ELPDLOO3j[idx,-1])
ELPDLOO4[idx,] <- as.matrix(ELPDLOO4j[idx,-1])

ELPDWAIC1[idx,] <- as.matrix(ELPDWAIC1j[idx,-1])
ELPDWAIC2[idx,] <- as.matrix(ELPDWAIC2j[idx,-1])
ELPDWAIC3[idx,] <- as.matrix(ELPDWAIC3j[idx,-1])
ELPDWAIC4[idx,] <- as.matrix(ELPDWAIC4j[idx,-1])

}


nombres1 = c("b0","b1","b2", "invtau2","tau2","tau", 
             "invsig2[1]","invsig2[2]","sig2[1]","sig2[2]", "sigma[1]","sigma[2]")
nombres2 = c("b0","b1","b2", "invtau2","tau2","tau", 
             "invsig2[1]","invsig2[2]","sig2[1]","sig2[2]", "sigma[1]","sigma[2]")
nombres3 = c("b0","b1[1]","b1[2]","b1[3]","b1[4]","b1[5]","b1[6]", 
             "b2[1]","b2[2]","b2[3]","b2[4]","b2[5]","b2[6]", 
             "invtau2","tau2","tau", 
             "lambda[1]","lambda[2]","rho[1]","rho[2]",
             "invsig2[1]","invsig2[2]","sig2[1]","sig2[2]", "sigma[1]","sigma[2]")
nombres4 = c("b0","b1[1]","b1[2]","b1[3]","b1[4]","b1[5]","b1[6]", 
             "b2[1]","b2[2]","b2[3]","b2[4]","b2[5]","b2[6]", 
             "invtau2","tau2","tau", 
             "lambda[1]","lambda[2]","rho[1]","rho[2]",
             "invsig2[1]","invsig2[2]","sig2[1]","sig2[2]", "sigma[1]","sigma[2]")


##################################################


### PARAMETROS
colnames(PARAM1) = nombres1
param1.m = apply(PARAM1,c(2),"mean")
param1.sd = apply(PARAM1,c(2),"sd")
param1.m
param1.sd

colnames(PARAM2) = nombres2
param2.m = apply(PARAM2,c(2),"mean")
param2.sd = apply(PARAM2,c(2),"sd")
param2.m
param2.sd

colnames(PARAM3) = nombres3
param3.m = apply(PARAM3,c(2),"mean")
param3.sd = apply(PARAM3,c(2),"sd")
param3.m
param3.sd

colnames(PARAM4) = nombres4
param4.m = apply(PARAM4,c(2),"mean")
param4.sd = apply(PARAM4,c(2),"sd")
param4.m
param4.sd



### LOO-WAIC
###          Estimate        SE
### elpd_loo -357.8435 11.296562
### p_loo      32.7760  3.409914
### looic     715.6870 22.593125

apply(ELPDLOO1,c(2),"mean")
apply(ELPDLOO1,c(2),"sd")
apply(ELPDLOO2,c(2),"mean")
apply(ELPDLOO2,c(2),"sd")
apply(ELPDLOO3,c(2),"mean")
apply(ELPDLOO3,c(2),"sd")
apply(ELPDLOO4,c(2),"mean")
apply(ELPDLOO4,c(2),"sd")

apply(ELPDWAIC1,c(2),"mean")
apply(ELPDWAIC1,c(2),"sd")
apply(ELPDWAIC2,c(2),"mean")
apply(ELPDWAIC2,c(2),"sd")
apply(ELPDWAIC3,c(2),"mean")
apply(ELPDWAIC3,c(2),"sd")
apply(ELPDWAIC4,c(2),"mean")
apply(ELPDWAIC4,c(2),"sd")


### ERROR
### mae   mse   rmse  
apply(ERROR1,c(2),"mean")
apply(ERROR1,c(2),"sd")
apply(ERROR2,c(2),"mean")
apply(ERROR2,c(2),"sd")
apply(ERROR3,c(2),"mean")
apply(ERROR3,c(2),"sd")
apply(ERROR4,c(2),"mean")
apply(ERROR4,c(2),"sd")


##################################################






