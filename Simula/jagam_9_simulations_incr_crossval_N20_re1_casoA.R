##################################################
### jagam_9_simulations_incr_crossval_N20_re1_casoB
##################################################
library(VGAM)
library(mgcv)
library(rstan)
library(rjags)
load.module("glm")
library(plot3D)
library(dplyr)
library(ggplot2)
library(corrplot)
library(splines)
library(splines2)
library(cgam)
library(MASS)
library(msm)
library(nlme)

library(dlookr)
library(gtools)
library(gdata)
library(ggeffects)
library(sjmisc)
library(mnormt)
library(loo)  
### http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/

dir <- "~/Documents/Sabatico/Examples/Simula_f1increasing 231028/Resultados N20 re1 casoB/"


### Specifications to simulate data 
set.seed(12345)
SIM = 100
TT <- 8   # number of time for repetitions 
N <- 20   # number of subjects 

k1 = 6
k2 = 6 
cuantiles1 = c(0.2,0.4, 0.6,0.8)
cuantiles2 =  c(0.2,0.4, 0.6,0.8)

### Subjects and time
id <- rep(1:N, each=TT)   # identify subjects
time <- rep(1:TT, times=N)  # time

f1 <- function(x) 10 + 0.1*(((x-5))^3)    ### Si es monotona 
f1error <- function(x) 10 + 0.1*(((x-5))^3) + rnorm(length(x),0,0.3)   ### Si es monotona, mas un error  
f2 <- function(x) x*((x-1)^2)*((x^2)*(10^2) - (3^4)*((x-1)^2))   ### cualquier otra forma


### Fixed effects
beta0 <- 6 

### Random effects
sigma2 <- c(0.7,0.05)

### Variance 
tau2 <- 1.7   # variance of the latent classes 


### Funciones 
### Right Normal Truncada I[Z < trb]
rnormright <- function(trb,mu,sig){
  rp <- pnorm(trb, mean=mu, sd=sig)
  u <- rp*runif(1)
  q <- qnorm(u, mean=mu, sd=sig)
  if(!is.finite(q)){ q = trb }
  return(q)
}
### Left Normal Truncada I[tra < Z]
rnormleft <- function(tra,mu,sig){
  rp <- pnorm(tra, mean=mu, sd=sig)
  u <- rp + (1-rp)*runif(1)
  q <- qnorm(u, mean=mu, sd=sig)
  if(!is.finite(q)){ q = tra }
  return(q)
}

ERRORES <- function(obs,est){
  error = obs-est
  error2 = error^2
  mae = mean(abs(error))
  mse = mean(error2)
  rmse = sqrt(mse)
  return(c("mae"=mae, "mse"=mse, "rmse"=rmse))
}


### Subjects and time
time1 <- time ###+ jitter(rep(0,N*TT),factor=10)
time2 <- time^2 
time3 <- time^3 
time4 <- time^4 

### Indicator for the beginning of observations for each subject for long tables format 
offset <- c(1,cumsum(table(id))+1)
n = N*TT
Ni = c(0,cumsum(table(id)))+1 ### offset


# 4. Definir la penalización $S1$ y $S2$ 
#Este es el código que produce la matriz de diferenciación.  
#No es el óptimo, pero funciona.  
#“k” es el número de b-splines y 
#“d” el orden de la diferenciación. 
#Adjunto el artículo donde discutimos esto (página 7).

diffMatrix = function(k, d = 2){
  if( (d<1) || (d %% 1 != 0) )stop("d must be a positive integer value"); 
  if( (k<1) || (k %% 1 != 0) )stop("k must be a positive integer value"); 
  if(d >= k)stop("d must be lower than k"); 
  out = diag(k); 
  for(i in 1:d){ 
    out = diff(out); 
  } 
  return(out) 
} 
(D1 = diffMatrix(k=k1, d=2))
(D2 = diffMatrix(k=k2, d=2))
(S1 = t(D1)%*%D1 + diag(1,k1)*10e-4) 
(S2 = t(D2)%*%D2 + diag(1,k2)*10e-4) 



### http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/

param.lme = c("b0","b1","b2", "invtau2","tau2","tau", "invsig2","sig2","sigma")

param.lme.add = c("b0","b1","b2", "invtau2","tau2","tau", "lambda","rho", "invsig2","sig2","sigma")

inits.lme.add.non.re1 <- function(){	list( 
  "b0" = rnorm(1,0,0.1) ,
  "b1" = rnorm(k1,0,0.1) ,
  "b2" = rnorm(k2,0,0.1) ,
  "invtau2" = rgamma(1,1,1) ,
  "lambda" = rgamma(2,1,1) ,
  "invsig2" = rgamma(2,1,1) 
)	} 

inits.lme.add.incr.re1 <- function(){	list( 
  "b0" = rnorm(1,0,0.1) ,
  "b1" = abs(rnorm(k1,0,0.1)) ,
  "b2" = rnorm(k2,0,0.1) ,
  "invtau2" = rgamma(1,1,1) ,
  "lambda" = rgamma(2,1,1) ,
  "invsig2" = rgamma(2,1,1) , 
  "bre1" = abs(rnorm(N,0,0.1))
)	} 


ERROR = array(NA, dim=c(SIM,3*5,4))
ELPDLOO = array(NA, dim=c(SIM,6,4))
ELPDWAIC = array(NA, dim=c(SIM,6,4))

PARAM1 = array(NA, dim=c(SIM,12,4))
PARAM3 = array(NA, dim=c(SIM,12+2+k1+k2,4))

### http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/ 

for(sim in 1:SIM){
  
  ### SIMULAR DATOS  
  
  x1 <- rnorm(N*TT, time, 0.2)
  x2 <- runif(N*TT, 0, 1)
  
  ### Random effects
  gama0 <- rnorm(N,0,sqrt(sigma2[1]))  # random effects 
  gama1 <- abs(rnorm(N,0,sqrt(sigma2[2])))  # random effects 
  
  ### Linear predictor and response 
  
  f1_obs = f1(x1)
  f2_obs = f2(x2)
  f1_con_error = f1error(x1)
  f1diff = f1_con_error-f1(x1)
  
  
  ### Linear predictor and response 
  eta <- rep(0,N*TT)   # linear predictor 
  Ytrue <- rep(NA,N*TT)   # response variable 
  etaerror = Yerror = rep(0,N*TT)   # error
  for(i in 1:N){		
    for(t in offset[i]){	
      eta[t] <- beta0 + f1(x1[t]) + f2(x2[t]) + gama0[id[t]] + gama1[id[t]]*time[t]   ### random slope
      Ytrue[t] <- rnorm(1, eta[t], sqrt(tau2))
      etaerror[t] <- eta[t] + f1diff[t]    ### random intercept 
      Yerror[t] <- Ytrue[t] + f1diff[t] 
    }
    for(t in (offset[i]+1):(offset[i+1]-1)){	
      eta[t] <- beta0 + f1(x1[t]) + f2(x2[t]) + gama0[id[t]] + gama1[id[t]]*time[t]   ### random slope
      Ytrue[t] <- rnorm(1, eta[t], sqrt(tau2))
      etaerror[t] <- eta[t] + f1diff[t]    ### random intercept 
      Yerror[t] <- Ytrue[t] + f1diff[t] 
    }
  }
  
  
  # 2. Generar la matriz diseño $X$ para los B-splines
  # Generate a basis matrix for Natural Cubic Splines 
  knots1 = quantile(x1, cuantiles1)
  knots2 = quantile(x2, cuantiles2) 
  X1 <- ns(x = x1, knots = knots1, intercept = TRUE) 
  X2 <- ns(x = x2, knots = knots2, intercept = TRUE) 
  
  # 3. Generar la matriz diseño $XI1$ para los I-splines
  ### ibs: integrated basis splines 
  ### degree = 3 cubic splines
  XI1 <- ibs(x1, knots = knots1, degree = 1, intercept = TRUE) 
  XI2 <- ibs(x2, knots = knots2, degree = 1, intercept = TRUE) 
  
  
  ### ESTIMAR 
  
  ## 5.2 LME: Lineal fit without constraints:
  
  datos.lme <- list( y = Yerror , 
                     n = length(Ytrue) , N = N , Ni = Ni, 
                     x1 = x1 , x2 = x2 , time = time , id = id )  
  
  
  fit.lme.non.re1 <- stan("jagam_9_lme_non_reslope.stan",
                          data=datos.lme,
                          chains=3, warmup=1000, iter=2000, thin=2, cores=4 )  
  
  mu = get_posterior_mean(fit.lme.non.re1,"mu")
  error.lme.non.re1 = ERRORES(Ytrue,mu[,"mean-all chains"]) 
  mufix = get_posterior_mean(fit.lme.non.re1,"mufix")
  errorfix.lme.non.re1 = ERRORES(Ytrue,mufix[,"mean-all chains"]) 
  mu1 = get_posterior_mean(fit.lme.non.re1,"mu1")
  error1.lme.non.re1 = ERRORES(f1_obs,mu1[,"mean-all chains"]) 
  error1error.lme.non.re1 = ERRORES(f1_con_error,mu1[,"mean-all chains"]) 
  mu2 = get_posterior_mean(fit.lme.non.re1,"mu2")
  error2.lme.non.re1 = ERRORES(f2_obs,mu2[,"mean-all chains"]) 
  
  
  estima1 = get_posterior_mean(fit.lme.non.re1, param.lme)[,"mean-all chains"]
  
  
  ## 6.2 LME: Lineal con restricciones creciente  :
  
  fit.lme.incr.re1 <- stan("jagam_9_lme_incr_reslope.stan",
                           data=datos.lme,
                           chains=3,warmup=1000,iter=2000,thin=2,cores=4 ) 
  
  mu = get_posterior_mean(fit.lme.incr.re1,"mu")
  error.lme.incr.re1 = ERRORES(Ytrue,mu[,"mean-all chains"]) 
  mufix = get_posterior_mean(fit.lme.incr.re1,"mufix")
  errorfix.lme.incr.re1 = ERRORES(Ytrue,mufix[,"mean-all chains"]) 
  mu1 = get_posterior_mean(fit.lme.incr.re1,"mu1")
  error1.lme.incr.re1 = ERRORES(f1_obs,mu1[,"mean-all chains"]) 
  error1error.lme.incr.re1 = ERRORES(f1_con_error,mu1[,"mean-all chains"]) 
  mu2 = get_posterior_mean(fit.lme.incr.re1,"mu2")
  error2.lme.incr.re1 = ERRORES(f2_obs,mu2[,"mean-all chains"]) 
  
  
  estima2 = get_posterior_mean(fit.lme.incr.re1, param.lme)[,"mean-all chains"]
  
  
  ## 7.2 LME: For a spline-based fit without constraints:
  
  datos.lme.add <- list( y = Yerror , 
                         id = id ,
                         n = length(Ytrue) , 
                         N = N , Ni = Ni,
                         k1=k1, k2=k2, 
                         XI1 = XI1, X2 = X2,  
                         time = time, 
                         zero = rep(0,1+k1+k2), 
                         S1=S1 , S2=S2 )  
  
  fit.lme.add.non.re1 <- stan("jagam_9_lme_add_non_reslope.stan",
                              data=datos.lme.add,
                              chains=3,warmup=1000,iter=2000,thin=2,cores=4,
                              init= inits.lme.add.non.re1) 
  
  mu = get_posterior_mean(fit.lme.add.non.re1,"mu")
  error.lme.add.non.re1 = ERRORES(Ytrue,mu[,"mean-all chains"]) 
  mufix = get_posterior_mean(fit.lme.add.non.re1,"mufix")
  errorfix.lme.add.non.re1 = ERRORES(Ytrue,mufix[,"mean-all chains"]) 
  mu1 = get_posterior_mean(fit.lme.add.non.re1,"mu1")
  error1.lme.add.non.re1 = ERRORES(f1_obs,mu1[,"mean-all chains"]) 
  error1error.lme.add.non.re1 = ERRORES(f1_con_error,mu1[,"mean-all chains"]) 
  mu2 = get_posterior_mean(fit.lme.add.non.re1,"mu2")
  error2.lme.add.non.re1 = ERRORES(f2_obs,mu2[,"mean-all chains"]) 
  
  
  estima3 = get_posterior_mean(fit.lme.add.non.re1, param.lme.add)[,"mean-all chains"]
  
  
  ## 8.2. LME: Spline con restricciones creciente
  
  fit.lme.add.incr.re1 <- stan("jagam_9_lme_add_incr_reslope.stan",
                               data=datos.lme.add,
                               chains=3,warmup=1000,iter=2000,thin=2,cores=4,
                               init= inits.lme.add.incr.re1) 
  
  mu = get_posterior_mean(fit.lme.add.incr.re1,"mu")
  error.lme.add.incr.re1 = ERRORES(Ytrue,mu[,"mean-all chains"]) 
  mufix = get_posterior_mean(fit.lme.add.incr.re1,"mufix")
  errorfix.lme.add.incr.re1 = ERRORES(Ytrue,mufix[,"mean-all chains"]) 
  mu1 = get_posterior_mean(fit.lme.add.incr.re1,"mu1")
  error1.lme.add.incr.re1 = ERRORES(f1_obs,mu1[,"mean-all chains"]) 
  error1error.lme.add.incr.re1 = ERRORES(f1_con_error,mu1[,"mean-all chains"]) 
  mu2 = get_posterior_mean(fit.lme.add.incr.re1,"mu2")
  error2.lme.add.incr.re1 = ERRORES(f2_obs,mu2[,"mean-all chains"]) 
  
  estima4 = get_posterior_mean(fit.lme.add.incr.re1, param.lme.add)[,"mean-all chains"]
  
  
  ### COMPARAR
  
  loo_sample_non_re1 = fit.lme.non.re1 # reslope
  loo_sample_incr_re1 = fit.lme.incr.re1 # reslope
  loo_sample_add_non_re1 = fit.lme.add.non.re1 # reslope
  loo_sample_add_incr_re1 = fit.lme.add.incr.re1 # reslope
  
  ### we have to extract those log-likelihood terms that we so carefully had Stan calculate for us:
  log_lik_non_re1 = extract_log_lik(loo_sample_non_re1, merge_chains = F)
  log_lik_incr_re1 = extract_log_lik(loo_sample_incr_re1, merge_chains = F)
  log_lik_add_non_re1 = extract_log_lik(loo_sample_add_non_re1, merge_chains = F)
  log_lik_add_incr_re1 = extract_log_lik(loo_sample_add_incr_re1, merge_chains = F)
  
  r_eff_non_re1 = relative_eff(log_lik_non_re1)
  r_eff_incr_re1 = relative_eff(log_lik_incr_re1)
  r_eff_add_non_re1 = relative_eff(log_lik_add_non_re1)
  r_eff_add_incr_re1 = relative_eff(log_lik_add_incr_re1)
  
  ###  look at the results for each model, first the one with mu estimated:
  loo_non_re1 <- loo(log_lik_non_re1, r_eff=r_eff_non_re1) 
  loo_incr_re1 <- loo(log_lik_incr_re1, r_eff=r_eff_incr_re1) 
  loo_add_non_re1 <- loo(log_lik_add_non_re1, r_eff=r_eff_add_non_re1) 
  loo_add_incr_re1 <- loo(log_lik_add_incr_re1, r_eff=r_eff_add_incr_re1) 
  
  waic_non_re1 <- waic(log_lik_non_re1, r_eff=r_eff_non_re1) 
  waic_incr_re1 <- waic(log_lik_incr_re1, r_eff=r_eff_incr_re1) 
  waic_add_non_re1 <- waic(log_lik_add_non_re1, r_eff=r_eff_add_non_re1) 
  waic_add_incr_re1 <- waic(log_lik_add_incr_re1, r_eff=r_eff_add_incr_re1) 
  
  
  PARAM1[sim,,1] <- estima1
  PARAM1[sim,,2] <- estima2
  PARAM3[sim,,3] <- estima3
  PARAM3[sim,,4] <- estima4
  
  ELPDLOO[sim,,1] <- as.vector(loo_non_re1$estimates)   
  ELPDLOO[sim,,2] <- as.vector(loo_incr_re1$estimates)   
  ELPDLOO[sim,,3] <- as.vector(loo_add_non_re1$estimates)   
  ELPDLOO[sim,,4] <- as.vector(loo_add_incr_re1$estimates)   
  
  ELPDWAIC[sim,,1] <- as.vector(waic_non_re1$estimates)   
  ELPDWAIC[sim,,2] <- as.vector(waic_incr_re1$estimates)   
  ELPDWAIC[sim,,3] <- as.vector(waic_add_non_re1$estimates)   
  ELPDWAIC[sim,,4] <- as.vector(waic_add_incr_re1$estimates)   
  
  ERROR[sim,,1] <- c(error.lme.non.re1, errorfix.lme.non.re1, error1.lme.non.re1, error2.lme.non.re1, error1error.lme.non.re1) 
  ERROR[sim,,2] <- c(error.lme.incr.re1, errorfix.lme.incr.re1, error1.lme.incr.re1, error2.lme.incr.re1, error1error.lme.incr.re1) 
  ERROR[sim,,3] <- c(error.lme.add.non.re1, errorfix.lme.add.non.re1, error1.lme.add.non.re1, error2.lme.add.non.re1, error1error.lme.add.non.re1) 
  ERROR[sim,,4] <- c(error.lme.add.incr.re1, errorfix.lme.add.incr.re1, error1.lme.add.incr.re1, error2.lme.add.incr.re1, error1error.lme.add.incr.re1) 
  
  write.csv(PARAM1[,,1], paste0(dir,"PARAM.lme.non.re1.csv"))
  write.csv(PARAM1[,,2], paste0(dir,"PARAM.lme.incr.re1.csv"))
  write.csv(PARAM3[,,3], paste0(dir,"PARAM.lme.add.non.re1.csv"))
  write.csv(PARAM3[,,4], paste0(dir,"PARAM.lme.add.incr.re1.csv"))
  
  write.csv(ERROR[,,1], paste0(dir,"ERROR.lme.non.re1.csv"))
  write.csv(ERROR[,,2], paste0(dir,"ERROR.lme.incr.re1.csv"))
  write.csv(ERROR[,,3], paste0(dir,"ERROR.lme.add.non.re1.csv"))
  write.csv(ERROR[,,4], paste0(dir,"ERROR.lme.add.incr.re1.csv"))
  
  write.csv(ELPDLOO[,,1], paste0(dir,"ELPD.LOO.lme.non.re1.csv"))
  write.csv(ELPDLOO[,,2], paste0(dir,"ELPD.LOO.lme.incr.re1.csv"))
  write.csv(ELPDLOO[,,3], paste0(dir,"ELPD.LOO.lme.add.non.re1.csv"))
  write.csv(ELPDLOO[,,4], paste0(dir,"ELPD.LOO.lme.add.incr.re1.csv"))
  
  write.csv(ELPDWAIC[,,1], paste0(dir,"ELPD.WAIC.lme.non.re1.csv"))
  write.csv(ELPDWAIC[,,2], paste0(dir,"ELPD.WAIC.lme.incr.re1.csv"))
  write.csv(ELPDWAIC[,,3], paste0(dir,"ELPD.WAIC.lme.add.non.re1.csv"))
  write.csv(ELPDWAIC[,,4], paste0(dir,"ELPD.WAIC.lme.add.incr.re1.csv"))
  
}




pairs(cbind(time, x1,"f1"=f1(x1), x2,"f2"=f2(x2), "eta"=eta,Ytrue))

pairs(cbind(time, x1,"f1"=f1_con_error, x2,"f2"=f2(x2), "eta"=etaerror,Yerror))

plot(x1,f1(x1)) 
plot(x2,f2(x2))
plot(x1,f1_con_error)
plot(x1,f1diff)

# 2. Generar la matriz diseño $X$ para los B-splines
# Generate a basis matrix for Natural Cubic Splines 
matplot(x1, X1) 
matplot(x2, X2) 
# 3. Generar la matriz diseño $XI1$ para los I-splines
### ibs: integrated basis splines 
### degree = 3 cubic splines
matplot(x1, XI1) 
abline(v = knots1, h = knots1, lty = 2, col = "gray") 
matplot(x2, XI2) 
abline(v = knots2, h = knots2, lty = 2, col = "gray") 


### dev, portrait, 5.4 x 3.7
### Graphics 
data1 <- data.frame("id"=id, "Ytrue"=Ytrue, "eta"=eta,
                    "time"=time, "x1"=x1, "x2"=x2, 
                    "f1_con_error"=f1_con_error, "f1_obs"=f1_obs, 
                    "Yerror"=Yerror)
str(data1)

p0 <- ggplot(data1, aes(x=time, y=Ytrue, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=time, y=Ytrue, group=1),
              method="loess", se=TRUE, level=0.95, 
              lwd=1, alpha=0.5, color="red") + 
  ggtitle("Scenario 1") + 
  theme(plot.title = element_text(size=17), 
        axis.title.x=element_text(size=15),  
        axis.title.y=element_text(size=15), 
        axis.text=element_text(size=13)) +  
  theme(legend.position="bottom", 
        legend.text=element_text(size=13)) +  
  xlab(expression(paste("Time ",t))) + 
  ylab(expression(paste("Response ", y)))
p0


p1 <- ggplot(data1, aes(x=x1, y=Ytrue, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=x1, y=Ytrue, group=1),
              method="loess", se=TRUE, level=0.95, 
              lwd=1, alpha=0.5, color="red") + 
  ggtitle("Scenario 1") + 
  theme(plot.title = element_text(size=17), 
        axis.title.x=element_text(size=15),  
        axis.title.y=element_text(size=15), 
        axis.text=element_text(size=13)) +  
  theme(legend.position="bottom", 
        legend.text=element_text(size=13)) +  
  xlab(expression(paste("Variable ",x[1]))) + 
  ylab(expression(paste("Response ", y)))
p1

p2 <- ggplot(data1, aes(x=x2, y=Ytrue, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=x2, y=Ytrue, group=1),
              method="loess", se=TRUE, level=0.95,
              lwd=1, alpha=0.5, color="red") + 
  ggtitle("Scenario 1") + 
  theme(plot.title = element_text(size=17), 
        axis.title.x=element_text(size=15),  
        axis.title.y=element_text(size=15), 
        axis.text=element_text(size=13)) +  
  theme(legend.position="bottom", 
        legend.text=element_text(size=13)) +  
  xlab(expression(paste("Variable ",x[2]))) + 
  ylab(expression(paste("Response ", y)))
p2

p3 <- ggplot(data1, aes(x=time, y=eta, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=time, y=eta, group=1),
              method="loess", se=TRUE, level=0.95, 
              lwd=1, alpha=0.5, color="red") + 
  ggtitle("Scenario 1") + 
  theme(plot.title = element_text(size=17), 
        axis.title.x=element_text(size=15),  
        axis.title.y=element_text(size=15), 
        axis.text=element_text(size=13)) +  
  theme(legend.position="bottom", 
        legend.text=element_text(size=13)) +  
  xlab("time t") + 
  ylab(expression(paste("Linear predictor ", eta)))
p3


p4 <- ggplot(data1, aes(x=x1, y=f1_con_error, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=x1, y=f1_obs, group=1),
              method="loess", se=TRUE, level=0.100, 
              lwd=1, alpha=0.5, color="red") + 
  ggtitle("Scenario 1") + 
  theme(plot.title = element_text(size=17), 
        axis.title.x=element_text(size=15),  
        axis.title.y=element_text(size=15), 
        axis.text=element_text(size=13)) +  
  theme(legend.position="bottom", 
        legend.text=element_text(size=13)) +  
  xlab(expression(paste("Variable ",x[1]))) + 
  ylab(expression(paste(f[omega](x[1]))))
p4


p5 <- ggplot(data1, aes(x=x1, y=Yerror, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=x1, y=beta0+f1_obs+f2_obs, group=1),
              method="loess", se=TRUE, level=0.95, 
              lwd=1, alpha=0.5, color="red") + 
  ggtitle("Scenario 1") + 
  theme(plot.title = element_text(size=17), 
        axis.title.x=element_text(size=15),  
        axis.title.y=element_text(size=15), 
        axis.text=element_text(size=13)) +  
  theme(legend.position="bottom", 
        legend.text=element_text(size=13)) +
  xlab(expression(paste("Variable ",x[1]))) + 
  ylab(expression(paste("Response ", y)))
p5


p6 <- ggplot(data1, aes(x=x2, y=Yerror, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=x2, y=beta0+f1_obs+f2_obs, group=1),
              method="loess", se=TRUE, level=0.95, 
              lwd=1, alpha=0.5, color="red") + 
  ggtitle("Scenario 1") + 
  theme(plot.title = element_text(size=17), 
        axis.title.x=element_text(size=15),  
        axis.title.y=element_text(size=15), 
        axis.text=element_text(size=13)) +  
  theme(legend.position="bottom", 
        legend.text=element_text(size=13)) + 
  xlab(expression(paste("Variable ",x[2]))) + 
  ylab(expression(paste("Response ", y)))
p6




### PARAMETROS
param1.m = apply(PARAM1,c(2,3),"mean")
rownames(param1.m) = names(estima1)
param1.sd = apply(PARAM1,c(2,3),"sd")
rownames(param1.sd) = names(estima1)
param1.m
param1.sd

param3.m = apply(PARAM3,c(2,3),"mean")
rownames(param3.m) = names(estima3)
param3.sd = apply(PARAM3,c(2,3),"sd")
rownames(param3.sd) = names(estima3)
param3.m
param3.sd


### LOO-WAIC
###          Estimate        SE
### elpd_loo -357.8435 11.296562
### p_loo      32.7760  3.409914
### looic     715.6870 22.593125

print(loo_add_incr_re1$estimates)

apply(ELPDLOO,c(2,3),"mean")
apply(ELPDLOO,c(2,3),"sd")

apply(ELPDWAIC,c(2,3),"mean")
apply(ELPDWAIC,c(2,3),"sd")


### ERROR
### mae   mse   rmse  
print(error.lme.add.incr.re1) 
apply(ERROR,c(2,3),"mean")
apply(ERROR,c(2,3),"sd")









