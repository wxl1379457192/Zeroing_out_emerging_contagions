data{
  int <lower=1> m; // number of variant
  int <lower=1> lor; //length of each variant
  int <lower=1> lde; //length of each variant
  int <lower=1> lom; //length of each variant
  int <lower=1> n; // number of dataset
  int <lower=1> k1; // number of interest variables
  int <lower=1> k2; // number of control variables
  matrix[n, k1] X1; // interest variables
  matrix[n, k2] X2; // control variables
  vector[n] Rt;
  vector[n] Rt_sd;
  vector[m] r0;
}
parameters{
  real<lower=0> aic[k1,m]; //alpha for each country
  real bec[k2]; //beta for each country
  //real<lower=0> elic[m];// elipson for each country
  real<lower=2> R0_alpha;
  real<lower=2> R0_delta;
  real<lower=2> R0_omicron;
}
transformed parameters {
  matrix[k1,m] alpha;
  vector[k2] beta;
  //vector[15*n]eli;
  vector[n] mu;
  
  {
    // calculation of variants-based real R0
    // calculation of effective vaccined population for each record against local variants context
    for (j in 1:m){
      for(i in 1:k1){
        alpha[i,j] = aic[i,j]- (log(1.05)/6.0);
      }
    }
    for(i in 1:k2){
      beta[i] = bec[i];
    }
    mu[1:lor] = -X1[1:lor,] * alpha[,1] - X2[1:lor,] * beta;
    mu[(lor+1):lde] = -X1[(lor+1):lde,] * alpha[,2] - X2[(lor+1):lde,] * beta;
    mu[(lde+1):lom] = -X1[(lde+1):lom,] * alpha[,3] - X2[(lde+1):lom,] * beta;
   
  }
}

model{
/////////////////// variants based R0 ///////////////////
//////////////////  R0 baseline  ///////////////////////////////
  for (j in 1:m){
    // R0 - SARS-CoV-2
    aic[,j] ~ gamma(0.1667,1);
  }
  bec ~ normal(0,0.5);
  R0_alpha ~ gamma(r0[1],0.1);
  R0_delta ~ gamma(r0[2],0.1);
  R0_omicron ~ gamma(r0[3],0.1);
  for(i in 1:lor){
      Rt[i] ~ gamma((exp(mu[i])*R0_alpha),Rt_sd[i]);
  }
    for(i in (lor+1):lde){
      Rt[i] ~ gamma((exp(mu[i])*R0_delta),Rt_sd[i]);
  }
    for(i in (lde+1):lom){
      Rt[i] ~ gamma((exp(mu[i])*R0_omicron),Rt_sd[i]);
  }
}








