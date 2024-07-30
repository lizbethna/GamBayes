/*
Bayesian 
Ordinal Logistic regression 
Non constraints
*/

data{
  int<lower=0> n;
  int<lower=1,upper=4> y[n]; 
  real x1[n]; 
}

parameters{
  ordered[3] kappa; 
  real b1; 
}

transformed parameters{
  vector[n] mu;
  for(i in 1:n){
    mu[i] = x1[i]*b1 ; // expected response
  }
}

model {
  b1 ~ normal(0,9.9e+06);  
  for (i in 1:n) { 
    y[i] ~ ordered_logistic(mu[i], kappa);   // response    
  }
}

generated quantities {
  vector[n] log_lik;
  for (i in 1:n) {
    log_lik[i] = ordered_logistic_lpmf(y[i] | mu[i], kappa);
  }
}
