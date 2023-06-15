data{
  int <lower=1> m; // number of all records
  int <lower=1> k1; // number of NPI variables
  int <lower=1> k2; // number of vac variables
  matrix[m, k1] X1; // NPI variables
  matrix[m, k2] X2; // vac variables
  vector[m] Rt;
  vector[m] Rt_sd;
  real r0;
}
parameters{
  real<lower=0> alpha_hier[k1]; // real priors
  real beta_hier[k2]; // real priors
  //real<lower=0> k;
  //real<lower=0> u;
  real<lower=2> R0;
}
transformed parameters {
  vector[k1] alpha;
  vector[k2] beta;
  vector[m] mu;
  
  {
    for(i in 1:k1){
      alpha[i] = alpha_hier[i]- (log(1.05)/6.0);
      }
    for(i in 1:k2){
      beta[i] = beta_hier[i];
      }
    mu = -X1*alpha-X2*beta;
    //mu = -X1*alpha;
  }
  
}

model{
//////////////////  R0 baseline  ///////////////////////////////
  //k ~ normal(0,0.1);
  R0 ~ gamma(r0,0.1); // R0 - SARS-CoV-2
  //u ~ uniform(0,1);
///////////////////////////////////////////////////////
  alpha_hier ~ gamma(0.1667,1);
  beta_hier ~ normal(0,0.5);
  //eli_hier ~ normal(0,0.5);
  for(i in 1:m){
    Rt[i] ~ gamma((exp(mu[i])*R0),Rt_sd[i]);
  }
}

generated quantities {               
  vector[m] ypre;
  for(i in 1:m){
      ypre[i]=exp(mu[i])*R0;
    }
    
}


