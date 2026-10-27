// ///////////////////// 44 observations ////////////////////////////////

functions {
  array[] real seir_ode(
    real t,
    array[] real y,
    array[] real theta,
    array[] real x_r,
    array[] int x_i
  ) {
    real beta  = theta[1];
    real alpha = theta[2];
    real gamma = theta[3];

    real S = y[1];
    real E = y[2];
    real I = y[3];
    real R = y[4];

    real N = S + E + I + R;

    real dS_dt = -beta * S * I / N;
    real dE_dt =  beta * S * I / N - alpha * E;
    real dI_dt =  alpha * E - gamma * I;
    real dR_dt =  gamma * I;

    return { dS_dt, dE_dt, dI_dt, dR_dt };
  }
}

data {
  int<lower=1> n_days;
  int<lower=1> n_total_days;
  int<lower=1> burn_in_days;

  array[n_days] int<lower=1> obs_days;
  array[n_total_days + burn_in_days] real ts;

  real t0;
  array[n_days] int<lower=0> cases;
  real<lower=0> N;
}

transformed data {
  array[0] real x_r;
  array[0] int x_i;
}

parameters {
  real<lower=0> R_0;
  real<lower=0> alpha_inv;
  real<lower=0> gamma_inv;
  real<lower=0> phi_inv;
  real<lower=0, upper=1> p_reported;
  real<lower=1, upper=1e5> E0;
  real<lower=1, upper=1e5> I0;
}

transformed parameters {
  real alpha = 1.0 / alpha_inv;
  real gamma = 1.0 / gamma_inv;
  real beta  = R_0 * gamma;

  real S0 = N - E0 - I0;
  array[n_total_days + burn_in_days, 4] real y;

  // Interval incidence
  array[n_total_days] real incidence;

  real phi = 1.0 / phi_inv;

  if (S0 < 0)
    S0 = 1e-6;

  {
    array[4] real y0 = { S0, E0, I0, 0 };
    array[3] real theta = { beta, alpha, gamma };

    y = integrate_ode_bdf(
      seir_ode, y0, t0, ts, theta, x_r, x_i
    );
  }

  // Compute incidence from the exposed compartment
  for (t in 1:n_total_days) {
    real E_t = y[t + burn_in_days, 2];
    incidence[t] = fmax(1e-6, p_reported * alpha * E_t);
  }
}

model {
  R_0 ~ lognormal(1.0, 0.2);
  alpha_inv ~ lognormal(1.6, 0.13);
  gamma_inv ~ lognormal(1.6, 0.13);

  I0 ~ lognormal(log(5), 0.8);
  E0 ~ lognormal(log(10), 0.8);

  p_reported ~ beta(0.28, 100);
  phi_inv ~ exponential(10);

  // Likelihood on observed interval incidence
  for (i in 1:n_days) {
    int t = obs_days[i];

    if (t > 1)
      cases[i] ~ neg_binomial_2(incidence[t], phi);
  }
}

generated quantities {
  array[n_total_days] real Reff;
  array[n_total_days] int pred_cases;

  for (t in 1:n_total_days) {
    pred_cases[t] = neg_binomial_2_rng(incidence[t], phi);
    Reff[t] = R_0;
  }
}

