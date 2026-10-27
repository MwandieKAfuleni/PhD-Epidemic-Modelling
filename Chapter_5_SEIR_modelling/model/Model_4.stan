 

// ///////////////////    With 44 days  ////////////////////////////////////


// 2 Beta


functions {
  real get_beta(real t, int n_beta, array[] real beta_vec, array[] real t_changes) {
    int i = 1;
    while (i <= n_beta - 1 && t > t_changes[i]) {
      i += 1;
    }
    return beta_vec[i];
  }

  array[] real seir_tv_beta(real t, array[] real y, array[] real theta, array[] real x_r, array[] int x_i) {
    int N = x_i[1];
    int n_beta = x_i[2];

    real gamma = theta[1];
    real a     = theta[2];

    array[n_beta] real beta_vec;
    array[n_beta - 1] real t_changes;

    for (i in 1:n_beta)
      beta_vec[i] = theta[2 + i];
    for (i in 1:(n_beta - 1))
      t_changes[i] = theta[2 + n_beta + i];

    real beta = get_beta(t, n_beta, beta_vec, t_changes);

    real S = y[1];
    real E = y[2];
    real I = y[3];
    real R = y[4];

    real dS_dt = -beta * I * S / N;
    real dE_dt =  beta * I * S / N - a * E;
    real dI_dt = a * E - gamma * I;
    real dR_dt = gamma * I;

    return {dS_dt, dE_dt, dI_dt, dR_dt};
  }
}

data {
  int<lower=1> n_days;              // 44 calendar days
  int<lower=1> n_obs;               // 42 valid days
  array[n_obs] int<lower=1, upper=n_days> obs_idx; // indices of valid days

  real t0;
  array[n_days] real ts;
  int N;
  array[n_days] int cases;
  real<lower=1> i0;

  int<lower=1> n_beta;
  array[n_beta - 1] real<lower=0> t_changes;
}

transformed data {
  array[2] int x_i = { N, n_beta };
}

parameters {
  vector<lower=0>[n_beta] beta_vec;
  real<lower=0> a_inv;
  real<lower=0> gamma_inv;
  real<lower=1> e0;
  real<lower=0, upper=1> p_reported;
  real<lower=0> phi_inv;
}

transformed parameters {
  real gamma = 1 / gamma_inv;
  real a     = 1 / a_inv;
  array[n_days, 4] real y;
  array[n_days] real incidence;
  real phi = 1. / phi_inv;

  array[2 + n_beta + (n_beta - 1)] real theta;

  theta[1] = gamma;
  theta[2] = a;
  for (i in 1:n_beta)
    theta[2 + i] = beta_vec[i];
  for (i in 1:(n_beta - 1))
    theta[2 + n_beta + i] = t_changes[i];

  array[4] real y0 = {N - i0 - e0, e0, i0, 0};
  y = integrate_ode_rk45(seir_tv_beta, y0, t0, ts, theta, rep_array(0.0, 0), x_i);

  for (i in 1:n_days) {
    incidence[i] = a * y[i, 2] * p_reported;
  }
}

model {
  beta_vec ~ lognormal(-0.6094, 0.2);
  a_inv ~ lognormal(1.6, 0.13);
  gamma_inv ~ lognormal(1.6, 0.13);
  e0 ~ lognormal(log(10), 0.8);
  p_reported ~ beta(0.28, 100);
  phi_inv ~ exponential(10);

  // Likelihood only on valid days
  for (i in 1:n_obs)
    cases[obs_idx[i]] ~ neg_binomial_2(incidence[obs_idx[i]], phi);
}

generated quantities {
  array[n_beta] real R0;
  array[n_days] real Reff;
  array[n_days] real pred_cases;
  real recovery_time = gamma_inv;
  real incubation_time = a_inv;

  for (i in 1:n_beta)
    R0[i] = beta_vec[i] / gamma;

  for (i in 1:n_days)
    pred_cases[i] = neg_binomial_2_rng(incidence[i], phi);

  for (i in 1:n_days) {
    int interval = 1;
    while (interval <= n_beta - 1 && ts[i] > t_changes[interval]) {
      interval += 1;
    }
    Reff[i] = R0[interval];
  }
}






