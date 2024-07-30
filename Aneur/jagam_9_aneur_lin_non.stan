/*
Bayesian 
Linear regression 
Non constraints
*/

data{
  int<lower=0> n;
  real y[n]; 
  real x1[n]; 
}

parameters{
  real b0; 
  real b1; 
  real<lower=0> invtau2; 
}

transformed parameters{
  vector[n] mu;
  real<lower=0> tau2;
  real<lower=0> tau;
  for(i in 1:n){
    mu[i] = b0 + x1[i]*b1 ; // expected response
  }
  tau2 = 1/invtau2; // convert invtau2 to standard GLM tau2
  tau = pow(tau2,0.5);
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior 
  b0 ~ normal(0,9.9e+06); 
  b1 ~ normal(0,9.9e+06);  
  for (i in 1:n) { 
    y[i] ~ normal(mu[i], tau);   // response    
  }
}

generated quantities {
  vector[n] log_lik;
  for (i in 1:n) {
    log_lik[i] = normal_lpdf(y[i] | mu[i], tau);
  }
}
