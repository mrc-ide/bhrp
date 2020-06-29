data {
  int<lower=1> N; // days of observed data
  int<lower=1> N2; // days to run exogeneous component for
  int cases[N]; // reported cases
  real SI[N]; // fixed SI using empirical data
  
}
transformed data {
  vector[N] SI_rev; // SI in reverse order
  for(i in 1:N)
  SI_rev[i] = SI[N-i+1];
  
}
parameters {
  real<lower=0> phi;
  vector[N+1] weekly_effect;
  real<lower=0, upper=1> weekly_rho;
  real<lower=0, upper=1> weekly_rho1;
  real<lower=0> weekly_sd;
  vector<lower=0>[N2] mu;
}

transformed parameters {
  vector[N] prediction=rep_vector(1e-5,N);
  vector<lower=0>[N] Rt;
  
  {
    Rt[1:N] = exp(weekly_effect[1:N]);
    prediction[1:N2]= prediction[1:N2] +  mu[1:N2];
    for (i in 2:N) {
      real convolution = dot_product(prediction[1:(i-1)], tail(SI_rev, i-1));
      prediction[i] = prediction[i] + Rt[i] * convolution;
    }
    
  }
}
model {
  weekly_sd ~ normal(0,1);
  weekly_rho ~ normal(0.8, 0.05);
  weekly_rho1 ~ normal(0.1, 0.05);
  phi ~ normal(0,5);
  mu[1:N2] ~ exponential(0.5);
  
  weekly_effect[3:(N+1)] ~ normal( weekly_effect[2:N]* weekly_rho + weekly_effect[1:(N-1)]* weekly_rho1, 
  weekly_sd *sqrt(1-pow(weekly_rho,2)-pow(weekly_rho1,2) - 2 * pow(weekly_rho,2) * weekly_rho1/(1-weekly_rho1)));
  weekly_effect[2] ~ normal(-1,weekly_sd *sqrt(1-pow(weekly_rho,2)-pow(weekly_rho1,2) - 2 * pow(weekly_rho,2) * weekly_rho1/(1-weekly_rho1)));
  weekly_effect[1] ~ normal(-1, 0.1);
  
  cases ~ neg_binomial_2(prediction,phi);
}
