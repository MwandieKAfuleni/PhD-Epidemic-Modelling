




// ////////////////// with 44 days ///////////////////////////////////////


data {
  int<lower=1> n_days;                 // total days (44)
  int<lower=1> n_obs;                  // number of valid observations (42)
  array[n_obs] int<lower=1, upper=n_days> obs_idx;  // indices of valid days
  vector[n_days] time;                 // time for all 44 days
  array[n_days] int<lower=0> cases;          // cases for all 44 days
  real<lower=0> T_change;              // change point
}

parameters {
  real<lower=0> C0;
  real<lower=0, upper=0.5> r;
  real<lower=0> r_minus;
  real<lower=0> phi_inv;
}

transformed parameters {
  real phi = 1 / phi_inv;
  vector[n_days] log_mu;
  vector[n_days] mu;

  for (t in 1:n_days) {
    if (time[t] < T_change)
      log_mu[t] = log(C0) + r * time[t];
    else
      log_mu[t] = log(C0) + r * T_change - r_minus * (time[t] - T_change);
  }

  mu = exp(log_mu);
}

model {
  // Priors
  C0 ~ lognormal(log(5), 0.8);
  r ~ normal(0.2, 0.2);
  r_minus ~ normal(0.1, 0.1);
  phi_inv ~ exponential(10);

  // Likelihood only for valid days
  for (i in 1:n_obs)
    cases[obs_idx[i]] ~ neg_binomial_2(mu[obs_idx[i]], phi);
}

generated quantities {
  array[n_days] int pred_cases;
  vector[n_days] log_lik;

  for (t in 1:n_days) {
    pred_cases[t] = neg_binomial_2_rng(mu[t], phi);
    log_lik[t] = neg_binomial_2_lpmf(cases[t] | mu[t], phi);
  }
}








