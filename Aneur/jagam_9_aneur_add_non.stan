/*
Bayesian 
*/

data{
  int<lower=0> n;
  int<lower=0> k1;
  real y[n]; 
  matrix[n,k1] XI1; 
  vector[1+k1] zero;
  matrix[k1,k1] S1;
}

parameters{
  real b0; 
  vector[k1] b1; 
  real<lower=0> invtau2; 
  real<lower=0> lambda;
}

transformed parameters{
  vector[n] mu;
  real<lower=0> tau2;
  real<lower=0> tau;
  real rho; 
  matrix[k1,k1] K1;  
  for(i in 1:n){
    mu[i] =  b0 + dot_product(XI1[i,] , b1) ; // expected response
  }
  tau2 = 1/invtau2; // convert invtau2 to standard GLM tau2
  tau = pow(tau2,0.5);
  K1 = S1 * lambda ; 
  rho = log(lambda);
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior 
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate! 
  // prior for s(x1) and s(x2) 
  b0 ~ normal(0,9.9e+06); 
  b1 ~ multi_normal(zero[(1+1):(1+k1)],K1); 
  // smoothing parameter priors CHECK...
    lambda ~ gamma(.05,.005);
  for (i in 1:n) { 
    y[i] ~ normal(mu[i], tau);   // response    
  }
}

generated quantities {
  vector[n] log_lik;
  vector[n] mu1; 
  for(i in 1:n){
    log_lik[i] = normal_lpdf(y[i] | mu[i], tau);
    mu1[i] = dot_product(XI1[i,] , b1);  
  }
}


