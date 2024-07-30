/*
Bayesian 
Linear regression 
Increasing constraints
*/

data{
  int<lower=0> N;
  int<lower=0> n;
  int<lower=1> Ni[N+1];
  int<lower=1,upper=4> y[n]; 
  real x1[n]; 
  int<lower=1> id[n]; 
}

parameters{
  ordered[3] kappa; 
  real<lower=0> b1; 
  real<lower=0> bre1[N];   // random effects
  real<lower=0> invsig2; 
}

transformed parameters{
  vector[n] mu;
  real<lower=0> sig2;
  real<lower=0> sigma;
  for(i in 1:N){
    for(t in Ni[i]:(Ni[i+1]-1)){
      mu[t] = x1[t]*b1 + x1[t]*bre1[i] ; // expected response
    }
  }
  sig2 = 1/invsig2; // convert to standard GLM tau2
  sigma = pow(sig2,0.5);
}

model {
  invsig2 ~ gamma(.05,.005); // precision parameter prior 
  b1 ~ normal(0,9.9e+06); 
  for(i in 1:N){
    bre1[i] ~ normal(0, sigma);
    for(t in Ni[i]:(Ni[i+1]-1)){
      y[t] ~ ordered_logistic(mu[t], kappa);   // response    
    }
  }
}


generated quantities {
  vector[n] log_lik;
  for(i in 1:N){
    for(t in Ni[i]:(Ni[i+1]-1)){
      log_lik[t] = ordered_logistic_lpmf(y[t] | mu[t], kappa);
    }
  }
}



