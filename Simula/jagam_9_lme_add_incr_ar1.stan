/*
Bayesian 
*/

data{
  int<lower=0> N;
  int<lower=0> n;
  int<lower=0> Ni[N+1];
  vector[n] y; 
  int<lower=1> id[n]; 
  int<lower=1> k1;
  int<lower=1> k2;
  matrix[n,k1] XI1; 
  matrix[n,k2] X2;
  vector[1+k1+k2] zero;
  matrix[k1,k1] S1;
  matrix[k2,k2] S2;
}

parameters{
  real b0; 
  vector[k1] b1; 
  vector[k2] b2; 
  real bre0[N];   // random effects
  real<lower=0> bre1[N];   // random effects
  real<lower=0> invtau2; 
  vector<lower=0>[2] lambda;
  real<lower=0> invsig2[2]; 

}

transformed parameters{
  vector[n] mu;
  real<lower=0> tau2;
  real<lower=0> tau;
  vector[2] rho; 
  real<lower=0> sig2[2];
  real<lower=0> sigma[2];
  matrix[k1,k1] K1; 
  matrix[k2,k2] K2;  

  for(i in 1:N){
    for(t in Ni[i]:Ni[i]){
      mu[t] = b0 + dot_product(XI1[t,],b1) + dot_product(X2[t,],b2) + bre0[id[t]] ; // expected response  
    }
    for(t in (Ni[i]+1):(Ni[i+1]-1)){
      mu[t] = b0 + dot_product(XI1[t,],b1) + dot_product(X2[t,],b2) + bre0[id[t]] + bre1[id[t]]*(y[t-1]-b0) ; // expected response  
    }
  }
  tau2 = 1/invtau2; // convert invtau2 to standard GLM tau2
  tau = pow(tau2,0.5);
  K1 = S1 * lambda[1] ; 
  K2 = S2 * lambda[2] ; 
  rho = log(lambda);
  for(j in 1:2){
    sig2[j] = 1/invsig2[j]; // convert to standard GLM tau2
    sigma[j] = pow(sig2[j],0.5); 
  }
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior 
  invsig2 ~ gamma(.5,.05); // precision parameter prior 
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate! 
  // prior for s(x1) and s(x2) ... 
  b0 ~ normal(0,9.9e+06); 
  b1 ~ multi_normal(zero[(1+1):(1+k1)],K1); 
  b2 ~ multi_normal(zero[(1+1):(1+k2)],K2); 
  // smoothing parameter priors CHECK...
  for(j in 1:2){
    lambda[j] ~ gamma(.05,.005);
  }
  for(i in 1:N){
    bre0[i] ~ normal(0, sigma[1]);
    bre1[i] ~ normal(0, sigma[2]);
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

      mu1[t] = dot_product(XI1[t,] , b1);  
      mu2[t] = dot_product(X2[t,] , b2);   
    } 
    for(t in Ni[i]:Ni[i]){
      mufix[t] = b0 + dot_product(XI1[t,],b1) + dot_product(X2[t,],b2) ; // 
      mumix[t] = mufix[t] + bre0[id[t]] ; // 
    }
    for(t in (Ni[i]+1):(Ni[i+1]-1)){
      mufix[t] = b0 + dot_product(XI1[t,],b1) + dot_product(X2[t,],b2)  ; // 
      mumix[t] = mufix[t] + bre0[id[t]] + bre1[id[t]]*(y[t-1]-b0) ; // 
    }
  }
}


