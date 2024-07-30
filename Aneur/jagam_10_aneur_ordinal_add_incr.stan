/*
Bayesian 
Ordinal Logistic regression 
*/

data{
  int<lower=0> n;
  int<lower=0> k1;
  int<lower=1,upper=4> y[n]; 
  matrix[n,k1] XI1; 
  vector[1+k1] zero;
  matrix[k1,k1] S1;
}

parameters{
  ordered[3] kappa; 
  vector<lower=0>[k1] b1;  
  real<lower=0> lambda;
}

transformed parameters{
  vector[n] mu;
  real rho; 
  matrix[k1,k1] K1;   
  for(i in 1:n){
    mu[i] =  dot_product(XI1[i,] , b1) ; // expected response
  }
  K1 = S1 * lambda; 
  rho = log(lambda);
}

model {
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate! 
  // prior for s(x1) and s(x2)... 
  b1 ~ multi_normal(zero[(1+1):(1+k1)],K1); 
  // smoothing parameter priors CHECK...
    lambda ~ gamma(.05,.005);
  for (i in 1:n) { 
    y[i] ~ ordered_logistic(mu[i], kappa);   // response    
  }
}

generated quantities {
  vector[n] log_lik;
  vector[n] mu1;
  for(i in 1:n){
    log_lik[i] = ordered_logistic_lpmf(y[i] | mu[i], kappa);
    mu1[i] = dot_product(XI1[i,] , b1);   
  }
}


