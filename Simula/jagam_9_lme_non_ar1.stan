/*
Bayesian 
*/

data{
  int<lower=0> N;
  int<lower=0> n;
  int<lower=1,upper=n+1> Ni[N+1];
  vector[n] y; 
  vector[n] x1; 
  vector[n] x2; 
  int<lower=1,upper=N> id[n];  
}

parameters{
  real b0; 
  real b1; 
  real b2; 
  real bre0[N];   // random effects
  real bre1[N];   // random effects
  real<lower=0> invtau2; 
  real<lower=0> invsig2[2]; 

}

transformed parameters{
  vector[n] mu;
  real<lower=0> tau2;
  real<lower=0> tau;
  real<lower=0> sig2[2];
  real<lower=0> sigma[2]; 
 
  for(i in 1:N){
    for(t in Ni[i]:Ni[i]){
      mu[t] = b0 + x1[t]*b1 + x2[t]*b2 + bre0[id[t]] ; // expected response  
    }
    for(t in (Ni[i]+1):(Ni[i+1]-1)){
      mu[t] = b0 + x1[t]*b1 + x2[t]*b2 + bre0[id[t]] + bre1[id[t]]*(y[t-1]-b0) ; // expected response  
    }
  }
  tau2 = 1/invtau2; // convert invtau2 to standard GLM tau2
  tau = pow(tau2,0.5);
  for(j in 1:2){
    sig2[j] = 1/invsig2[j]; // convert to standard GLM tau2
    sigma[j] = pow(sig2[j],0.5);
  }
}

model { 
  invtau2 ~ gamma(.05,.005); // precision parameter prior 
  for(j in 1:2){
    invsig2[j] ~ gamma(.05,.005); // precision parameter prior 
  }
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate! 
  // prior for s(x1) ... 
  b0 ~ normal(0,9.9e+06); 
  b1 ~ normal(0,9.9e+06); 
  b2 ~ normal(0,9.9e+06); 
  for(i in 1:N){
    bre0[i] ~ normal(0,sigma[1]);
    bre1[i] ~ normal(0,sigma[2]);
    for(t in Ni[i]:(Ni[i+1]-1)){
      y[t] ~ normal(mu[t], tau);   // response    
    }
  }   
}

generated quantities {
  vector[n] log_lik;
  vector[n] mufix;  
  vector[n] mumix;
  for(i in 1:N){
    for(t in Ni[i]:(Ni[i+1]-1)){
      log_lik[t] = normal_lpdf(y[t] | mu[t], tau);
    } 
    for(t in Ni[i]:Ni[i]){
      mufix[t] = b0 + x1[t]*b1 + x2[t]*b2 ; // 
      mumix[t] = mufix[t] + bre0[id[t]] ; // 
    }
    for(t in (Ni[i]+1):(Ni[i+1]-1)){
      mufix[t] = b0 + x1[t]*b1 + x2[t]*b2 ; // 
      mumix[t] = mufix[t] + bre0[id[t]] + bre1[id[t]]*(y[t-1]-b0) ; // 
    }
  }
}


