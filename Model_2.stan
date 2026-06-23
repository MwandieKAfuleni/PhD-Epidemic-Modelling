// /////////////////////////////////////////////////////


// functions {
//   real[] seir_ode(
//     real t,
//     real[] y,
//     real[] theta,
//     real[] x_r,
//     int[] x_i
//   ) {
//     real beta  = theta[1];
//     real alpha = theta[2];
//     real gamma = theta[3];
//     real N     = x_r[1];
// 
//     real S = y[1];
//     real E = y[2];
//     real I = y[3];
//     real R = y[4];
// 
//     real dS_dt = -beta * S * I / N;
//     real dE_dt =  beta * S * I / N - alpha * E;
//     real dI_dt =  alpha * E - gamma * I;
//     real dR_dt =  gamma * I;
// 
//     return { dS_dt, dE_dt, dI_dt, dR_dt };
//   }
// }
// 
// data {
//   int<lower=1> n_days;
//   int<lower=0> burn_in_days;
//   real<lower=0> t0;
//   real ts[n_days + burn_in_days];
//   int cases[n_days];
// 
//   real<lower=0> N;
// 
//   // fixed epidemiological parameters
//   real<lower=0> alpha;
//   real<lower=0> gamma;
// }
// 
// transformed data {
//   real x_r[1] = { N };
//   int x_i[0];
// }
// 
// parameters {
//   real<lower=0> R_0;
//   real<lower=0> I0;
//   real<lower=0> epsilon;
//   real<lower=0, upper=1> p_reported;
//   real<lower=0> phi_inv;
// }
// 
// 
// transformed parameters {
//   real beta = R_0 * gamma;
// 
//   /*
//     Eigenvector-based centre for epsilon = I0 / E0.
//   */
//   real epsilon_center =
//     2.0 * alpha /
//     (
//       -alpha + gamma
//       + sqrt(square(alpha - gamma) + 4.0 * alpha * beta)
//     );
// 
//   real E0 = I0 / epsilon;
//   real S0 = N - E0 - I0;
// 
//   real y[n_days + burn_in_days, 4];
//   real incidence[n_days];
//   real phi = 1.0 / phi_inv;
// 
//   {
//     real y0[4] = { S0, E0, I0, 0 };
//     real theta[3] = { beta, alpha, gamma };
// 
//     y = integrate_ode_bdf(
//       seir_ode, y0, t0, ts,
//       theta, x_r, x_i
//     );
// 
//     for (t in 1:n_days) {
//       int idx = t + burn_in_days;
//       incidence[t] =
//         fmax(1e-9, p_reported * alpha * y[idx, 2]);
//     }
//   }
// }
// 
// model {
//   // Priors
//   R_0     ~ normal(2.5, 0.5);
//   I0      ~ lognormal(log(50), 0.8);
// 
//   // Eigenvector‑based prior
//   epsilon ~ lognormal(log(epsilon_center), 0.5);
// 
//   p_reported ~ beta(2, 38);    //(8, 20);
//   phi_inv ~ exponential(5);
// 
//   // Likelihood
//   cases ~ neg_binomial_2(incidence, phi);
// }
// 
// generated quantities {
//   vector[n_days] log_lik;
//   int pred_cases[n_days];
// 
//   for (t in 1:n_days) {
//     log_lik[t] =
//       neg_binomial_2_lpmf(cases[t] | incidence[t], phi);
//     pred_cases[t] =
//       neg_binomial_2_rng(incidence[t], phi);
//   }
// }
// 
// 


// ///////////////   42 observations   ///////////////////////////////////////////
// functions {
//   real[] seir_ode(
//     real t,
//     real[] y,
//     real[] theta,
//     real[] x_r,
//     int[] x_i
//   ) {
// 
//     real beta  = theta[1];
//     real alpha = theta[2];
//     real gamma = theta[3];
//     real N     = x_r[1];
// 
//     real S = y[1];
//     real E = y[2];
//     real I = y[3];
//     real R = y[4];
// 
//     real dS_dt = -beta * S * I / N;
//     real dE_dt =  beta * S * I / N - alpha * E;
//     real dI_dt =  alpha * E - gamma * I;
//     real dR_dt =  gamma * I;
// 
//     return {dS_dt, dE_dt, dI_dt, dR_dt};
//   }
// }
// 
// data {
// 
//   int<lower=1> n_days;
//   int<lower=0> burn_in_days;
// 
//   real<lower=0> t0;
// 
//   real ts[n_days + burn_in_days];
// 
//   int cases[n_days];
// 
//   real<lower=0> N;
// 
//   // fixed epidemiological parameters
//   real<lower=0> alpha;
//   real<lower=0> gamma;
// }
// 
// transformed data {
// 
//   real x_r[1] = {N};
// 
//   int x_i[0];
// }
// 
// parameters {
// 
//   real<lower=0> R_0;
// 
//   real<lower=0> I0;
// 
//   real<lower=0, upper=1> p_reported;
// 
//   real<lower=0> phi_inv;
// }
// 
// transformed parameters {
// 
//   real beta = R_0 * gamma;
// 
//   // Eigenvector relationship
//   real epsilon =
//     2.0 * alpha /
//     (
//       -alpha + gamma
//       + sqrt(square(alpha - gamma)
//              + 4.0 * alpha * beta)
//     );
// 
//   real E0 = I0 / epsilon;
// 
//   real S0 = N - E0 - I0;
// 
//   real y[n_days + burn_in_days, 4];
// 
//   real incidence[n_days];
// 
//   real phi = 1.0 / phi_inv;
// 
//   {
//     real y0[4] = {S0, E0, I0, 0};
// 
//     real theta[3] = {beta, alpha, gamma};
// 
//     y = integrate_ode_bdf(seir_ode,y0,t0,ts,theta,x_r,x_i);
// 
//     for (t in 1:n_days) {
// 
//       int idx = t + burn_in_days;
// 
//       incidence[t] = fmax(1e-9, p_reported * alpha * y[idx, 2]);
//     }
//   }
// }
// 
// model {
// 
//   // Priors
//   R_0 ~ lognormal(1.0, 0.2); //normal(2.5, 0.5);
// 
//   I0 ~ lognormal(log(5), 0.8); //lognormal(log(2), 0.8);
// 
//   p_reported ~  beta(0.28, 100); //beta(2, 38);
// 
//   phi_inv ~ exponential(10);
// 
//   // Likelihood
//   cases ~ neg_binomial_2(incidence, phi);
// }
// 
// generated quantities {
// 
//   real epsilon_out =
//     2.0 * alpha /
//     (
//       -alpha + gamma
//       + sqrt(square(alpha - gamma)
//              + 4.0 * alpha * beta)
//     );
// 
//   real E0_out = I0 / epsilon_out;
// 
//   //vector[n_days] log_lik;
// 
//   int pred_cases[n_days];
// 
//   for (t in 1:n_days) {
// 
//     //log_lik[t] =neg_binomial_2_lpmf(cases[t] |incidence[t],phi);
// 
//     pred_cases[t] = neg_binomial_2_rng(incidence[t],phi);
//   }
// }





//    ////////////////////// 44 observations  ///////////////////
functions {
  real[] seir_ode(
    real t,
    real[] y,
    real[] theta,
    real[] x_r,
    int[] x_i
  ) {

    real beta  = theta[1];
    real alpha = theta[2];
    real gamma = theta[3];
    real N     = x_r[1];

    real S = y[1];
    real E = y[2];
    real I = y[3];
    real R = y[4];

    real dS_dt = -beta * S * I / N;
    real dE_dt =  beta * S * I / N - alpha * E;
    real dI_dt =  alpha * E - gamma * I;
    real dR_dt =  gamma * I;

    return {dS_dt, dE_dt, dI_dt, dR_dt};
  }
}

data {

  int<lower=1> n_days;
  int<lower=0> burn_in_days;

  real<lower=0> t0;

  real ts[n_days + burn_in_days];

  int cases[n_days];

  real<lower=0> N;

  // fixed epidemiological parameters
  real<lower=0> alpha;
  real<lower=0> gamma;

  // NEW: outlier handling
  int<lower=1> n_obs;
  int<lower=1, upper=n_days> obs_idx[n_obs];
}

transformed data {

  real x_r[1] = {N};

  int x_i[0];
}

parameters {

  real<lower=0> R_0;

  real<lower=0> I0;

  real<lower=0, upper=1> p_reported;

  real<lower=0> phi_inv;
}

transformed parameters {

  real beta = R_0 * gamma;

  // Eigenvector relationship
  real epsilon =
    2.0 * alpha /
    (
      -alpha + gamma
      + sqrt(square(alpha - gamma)
             + 4.0 * alpha * beta)
    );

  real E0 = I0 / epsilon;

  real S0 = N - E0 - I0;

  real y[n_days + burn_in_days, 4];

  real incidence[n_days];

  real phi = 1.0 / phi_inv;

  {
    real y0[4] = {S0, E0, I0, 0};

    real theta[3] = {beta, alpha, gamma};

    y = integrate_ode_bdf(seir_ode,y0,t0,ts,theta,x_r,x_i);

    for (t in 1:n_days) {

      int idx = t + burn_in_days;

      incidence[t] = fmax(1e-9, p_reported * alpha * y[idx, 2]);
    }
  }
}

model {

  // Priors
  R_0 ~ lognormal(1.0, 0.2);

  I0 ~ lognormal(log(5), 0.8);

  p_reported ~  beta(0.28, 100);

  phi_inv ~ exponential(10);

  // Likelihood: ONLY non-outlier days
  for (j in 1:n_obs) {
    int k = obs_idx[j];
    cases[k] ~ neg_binomial_2(incidence[k], phi);
  }
}

generated quantities {

  real epsilon_out =
    2.0 * alpha /
    (
      -alpha + gamma
      + sqrt(square(alpha - gamma)
             + 4.0 * alpha * beta)
    );

  real E0_out = I0 / epsilon_out;

  int pred_cases[n_days];

  for (t in 1:n_days) {
    pred_cases[t] = neg_binomial_2_rng(incidence[t],phi);
  }
}

