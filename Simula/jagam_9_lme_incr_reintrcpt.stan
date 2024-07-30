/*
Bayesian 
Linear regression 
Increasing constraints
*/

data{
  int<lower=0> N;
  int<lower=0> n;
  int<lower=1> Ni[N+1];
  real y[n]; 
  real x1[n]; 
  real x2[n]; 
  int<lower=1> id[n]; 
}

parameters{
  real b0; 
  real<lower=0> b1; 
  real b2; 
  real bre0[N];   // random effects
  real<lower=0> invtau2; 
  real<lower=0> invsig2; 
}

transformed parameters{
  vector[n] mu;
  real<lower=0> tau2;
  real<lower=0> tau;
  real<lower=0> sig2;
  real<lower=0> sigma;
  for(i in 1:N){
    for(t in Ni[i]:(Ni[i+1]-1)){
      mu[t] = b0 + x1[t]*b1 + x2[t]*b2 + bre0[id[t]] ; // expected response
    }
  }
  tau2 = 1/invtau2; // convert invtau2 to standard GLM tau2
  tau = pow(tau2,0.5);
  sig2 = 1/invsig2; // convert to standard GLM tau2
  sigma = pow(sig2,0.5);
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior 
  invsig2 ~ gamma(.05,.005); // precision parameter prior 
  b0 ~ normal(0,9.9e+06); 
  b1 ~ normal(0,9.9e+06); 
  b2 ~ normal(0,9.9e+06); 
  for(i in 1:N){
    bre0[i] ~ normal(0, sigma);
    for(t in Ni[i]:(Ni[i+1]-1)){
      y[t] ~ normal(mu[t], tau);   // response    
    }
  }
}

generated quantities {
  vector[n] log_lik;
  vector[n] mu1;  
  vector[n] mu2;  
  vector[n] mufix;  
  vector[n] mumix;
  for(i in 1:N){
    for(t in Ni[i]:(Ni[i+1]-1)){
      log_lik[t] = normal_lpdf(y[t] | mu[t], tau);
      mu1[t] = x1[t]*b1 ; 
      mu2[t] = x2[t]*b2 ; 
      mufix[t] = b0 + x1[t]*b1 + x2[t]*b2 ; // 
      mumix[t] = mufix[t] + bre0[id[t]] ; // 
    }
  }
}

