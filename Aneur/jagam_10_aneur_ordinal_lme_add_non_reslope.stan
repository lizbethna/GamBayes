/*
Bayesian 
*/

data{
  int<lower=0> N;
  int<lower=0> n;
  int<lower=0> Ni[N+1];
  int<lower=1,upper=4> y[n]; 
  int<lower=1> id[n]; 
  int<lower=0> k1;
  matrix[n,k1] XI1;  
  real x1[n];  
  vector[1+k1] zero;
  matrix[k1,k1] S1;
}

parameters{
  ordered[3] kappa; 
  vector[k1] b1; 
  real bre1[N];   // random effects
  real<lower=0> lambda;
  real<lower=0> invsig2; 
}

transformed parameters{
  vector[n] mu;
  real rho; 
  real<lower=0> sig2;
  real<lower=0> sigma;
  matrix[k1,k1] K1; 
  for(i in 1:N){
    for(t in Ni[i]:(Ni[i+1]-1)){
      mu[t] = dot_product(XI1[t,] , b1) + x1[t]*bre1[i] ; 
 // expected response
    }
  }
  K1 = S1 * lambda ;
  rho = log(lambda);
  sig2 = 1/invsig2; // convert to standard GLM tau2
  sigma = pow(sig2,0.5); 
}

model {
  invsig2 ~ gamma(.5,.05); // precision parameter prior 
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate! 
  // prior for s(x1) ... 
  b1 ~ multi_normal(zero[(1+1):(1+k1)],K1); 
  // smoothing parameter priors CHECK...
    lambda ~ gamma(.05,.005); 
  for(i in 1:N){
    bre1[i] ~ normal(0, sigma);
    for(t in Ni[i]:(Ni[i+1]-1)){
      y[t] ~ ordered_logistic(mu[t], kappa);   // response    
    }
  }
}

generated quantities {
  vector[n] log_lik;
  vector[n] mu1;
  vector[n] mu2;

  for(i in 1:N){
    for(t in Ni[i]:(Ni[i+1]-1)){
      log_lik[t] = ordered_logistic_lpmf(y[t] | mu[t], kappa);
      mu1[t] = dot_product(XI1[t,] , b1) ; // expected response
      mu2[t] = dot_product(XI1[t,] , b1) + x1[t]*bre1[i] ; 
    }
  }

}


