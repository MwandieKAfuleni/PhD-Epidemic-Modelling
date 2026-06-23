// /////////////////////////////////////////////////////


// Manuscript code
// 
// 
// 
// functions {
//   real[] seir_ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
//     real beta = theta[1];
//     real alpha = theta[2];
//     real gamma = theta[3];
// 
//     real S = y[1];
//     real E = y[2];
//     real I = y[3];
//     real R = y[4];
// 
//     real N = S + E + I + R;
// 
//     real dS_dt = -beta * S * I / N;
//     real dE_dt =  beta * S * I / N - alpha * E;
//     real dI_dt = alpha * E - gamma * I;
//     real dR_dt = gamma * I;
// 
//     return {dS_dt, dE_dt, dI_dt, dR_dt};
//   }
// }
// 
// data {
//   int<lower=1> n_days;
//   int<lower=1> burn_in_days;
//   real ts[n_days + burn_in_days];
//   real t0;
//   int cases[n_days];
//   real<lower=0> N;
// }
// 
// 
// transformed data {
//   real x_r[0];
//   int x_i[0];
// }
// 
// parameters {
//   real<lower=0> R_0;
//   real<lower=0> alpha;
//   real<lower=0> gamma;
//   real<lower=0> phi_inv;
//   real<lower=0, upper=1> p_reported;
//   real<lower=1> E0;
//   real<lower=1> I0;
// }
// 
// transformed parameters {
//   real beta = R_0 * gamma;
//   real S0 = N - E0 - I0;
//   real y[n_days + burn_in_days, 4];
//   real incidence[n_days - 1];
//   real phi = 1. / phi_inv;
//   real y0[4] = {S0, E0, I0, 0};
//   real theta[3] = {beta, alpha, gamma};
// 
//   y = integrate_ode_bdf(seir_ode, y0, t0, ts, theta, x_r, x_i);
// 
//   for (t in 1:(n_days - 1)) {
//     //incidence[t] = alpha * y[t, 2] * p_reported;
//     incidence[t] = fmax(1e-6, alpha * y[t + burn_in_days - 1, 2] * p_reported);
//   }
// }
// 
// model {
//   R_0 ~ lognormal(1.009619, 0.2);
//   alpha ~ lognormal(-1.68, 0.2);
//   gamma ~ lognormal(-1.986, 0.2);
//   phi_inv ~ exponential(10);
//   p_reported ~ beta(0.28, 99.72);
//   I0 ~ lognormal(7.699, 0.5);  // mean ≈ 2500
//   E0 ~ lognormal(8.169, 0.5);  // mean ≈ 4000
// 
//   cases[1:(n_days - 1)] ~ neg_binomial_2(incidence, phi);
// }
// 
// 
// 
// generated quantities {
//   real Reff[n_days];
//   int pred_cases[n_days - 1];
//   vector[n_days - 1] log_lik;
// 
//   for (t in 1:(n_days - 1)) {
// 
//     pred_cases[t] =
//       neg_binomial_2_rng(incidence[t], phi);
// 
//     log_lik[t] =
//       neg_binomial_2_lpmf(cases[t] | incidence[t], phi);
//   }
// 
//   for (t in 1:n_days) {
//     Reff[t] = beta / gamma;
//   }
// }
//









// ////// shorter epidemic window ///////////////////



// functions {
//   real[] seir_ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
//     real beta = theta[1];
//     real alpha = theta[2];
//     real gamma = theta[3];
// 
//     real S = y[1];
//     real E = y[2];
//     real I = y[3];
//     real R = y[4];
// 
//     real N = S + E + I + R;
// 
//     real dS_dt = -beta * S * I / N;
//     real dE_dt =  beta * S * I / N - alpha * E;
//     real dI_dt = alpha * E - gamma * I;
//     real dR_dt = gamma * I;
// 
//     return {dS_dt, dE_dt, dI_dt, dR_dt};
//   }
// }
// 
// data {
//   int<lower=1> n_days;
//   int<lower=1> burn_in_days;
//   real ts[n_days + burn_in_days];
//   real t0;
//   int cases[n_days];
//   real<lower=0> N;
// }
// 
// 
// transformed data {
//   real x_r[0];
//   int x_i[0];
// }
// 
// parameters {
//   real<lower=0> R_0;
//   real<lower=0> alpha_inv;
//   real<lower=0> gamma_inv;
//   real<lower=0> phi_inv;
//   real<lower=0, upper=1> p_reported;
//   real<lower=1> E0;
//   real<lower=1> I0;
// }
// 
// transformed parameters {
//   
//   real alpha = 1.0/alpha_inv;
//   real gamma = 1.0/gamma_inv;
//   real beta = R_0 * gamma;
//   
//   real S0 = N - E0 - I0;
//   real y[n_days + burn_in_days, 4];
//   real incidence[n_days];
//   real phi = 1. / phi_inv;
//   
//   for (t in 1:(n_days)) {
//     //incidence[t] = alpha * y[t, 2] * p_reported;
//     incidence[t] = fmax(1e-6, p_reported * alpha * y[t + burn_in_days, 2]);
//   }
//   real y0[4] = {S0, E0, I0, 0};
//   real theta[3] = {beta, alpha, gamma};
// 
//   y = integrate_ode_bdf(seir_ode, y0, t0, ts, theta, x_r, x_i);
// 
//   
// }
// 
// 
// 
// model {
//   R_0 ~ lognormal(1.0, 0.2);
//   alpha_inv ~ lognormal(1.6, 0.13); //lognormal(log(4.5), 0.2);
//   gamma_inv ~ lognormal(1.6, 0.13); //lognormal(log(5), 0.2);
//   // I0 ~ lognormal(3.4, 0.02); //30
//   // E0 ~ lognormal(4, 0.02); //50
//   // I0 ~ lognormal(2.9951, 0.0353); //20
//   // E0 ~ lognormal(3.4, 0.02); //30
//   I0 ~ lognormal(log(5), 0.8);
//   E0 ~ lognormal(log(10), 0.8);
//   
//   p_reported ~ beta(0.28, 100);
//   phi_inv ~ exponential(10);
//   
//   
//   //cases[1:(n_days - 1)] ~ neg_binomial_2(incidence, phi);
//   cases ~ neg_binomial_2(incidence, phi);
// }
// 
// 
// 
// generated quantities {
//   real Reff[n_days];
//   int pred_cases[n_days];
//   vector[n_days - 1] log_lik;
// 
//   for (t in 1:(n_days - 1)) {
// 
//     pred_cases[t] =
//       neg_binomial_2_rng(incidence[t], phi);
// 
//     log_lik[t] =
//       neg_binomial_2_lpmf(cases[t] | incidence[t], phi);
//   }
// 
//   for (t in 1:n_days) {
//     Reff[t] = beta / gamma;
//   }
// }
// 






// ///////////////// 42 days, does not consider removed data points////////////////////////////////
// functions {
//   real[] seir_ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
//     real beta  = theta[1];
//     real alpha = theta[2];
//     real gamma = theta[3];
// 
//     real S = y[1];
//     real E = y[2];
//     real I = y[3];
//     real R = y[4];
// 
//     real N = S + E + I + R;
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
//   int<lower=1> burn_in_days;
//   real ts[n_days + burn_in_days];
//   real t0;
//   int cases[n_days];
//   real<lower=0> N;
// }
// 
// transformed data {
//   real x_r[0];
//   int x_i[0];
// }
// 
// parameters {
//   real<lower=0> R_0;
//   real<lower=0> alpha_inv;
//   real<lower=0> gamma_inv;
//   real<lower=0> phi_inv;
//   real<lower=0, upper=1> p_reported;
//   real<lower=1, upper=1e5> E0;
//   real<lower=1, upper=1e5> I0;
// }
// 
// transformed parameters {
//   real alpha = 1.0 / alpha_inv;
//   real gamma = 1.0 / gamma_inv;
//   real beta  = R_0 * gamma;
// 
//   real S0 = N - E0 - I0;
//   real y[n_days + burn_in_days, 4];
//   real incidence[n_days];
//   real phi = 1.0 / phi_inv;
// 
//   // enforce non-negative S0 (softly)
//   // if you want to be strict, constrain E0 + I0 < N in parameters instead
//   if (S0 < 0)
//     S0 = 1e-6;
// 
//   {
//     real y0[4]    = { S0, E0, I0, 0 };
//     real theta[3] = { beta, alpha, gamma };
// 
//     y = integrate_ode_bdf(seir_ode, y0, t0, ts, theta, x_r, x_i);
//   }
// 
//   // Incidence: p_reported * alpha * E(t)
//   for (t in 1:n_days) {
//     incidence[t] = fmax(1e-6, p_reported * alpha * y[t + burn_in_days, 2]);
//   }
// }
// 
// model {
//   R_0 ~ lognormal(1.0, 0.2);
//   alpha_inv ~ lognormal(1.6, 0.13); //lognormal(log(4.5), 0.2);
//   gamma_inv ~ lognormal(1.6, 0.13); //lognormal(log(5), 0.2);
//   // I0 ~ lognormal(3.4, 0.02); //30
//   // E0 ~ lognormal(4, 0.02); //50
//   // I0 ~ lognormal(2.9951, 0.0353); //20
//   // E0 ~ lognormal(3.4, 0.02); //30
//   I0 ~ lognormal(log(5), 0.8);
//   E0 ~ lognormal(log(10), 0.8);
// 
//   p_reported ~ beta(0.28, 100);
//   phi_inv ~ exponential(10);
// 
//   cases ~ neg_binomial_2(incidence, 1.0 / phi_inv);
// }
// 
// generated quantities {
//   real Reff[n_days];
//   int pred_cases[n_days];
//   vector[n_days] log_lik;
//   //real gamma = 1.0 / gamma_inv;
//   //real beta  = R_0 * gamma;
//   //real phi   = 1.0 / phi_inv;
// 
//   for (t in 1:n_days) {
//     pred_cases[t] = neg_binomial_2_rng(incidence[t], phi);
//     log_lik[t]    = neg_binomial_2_lpmf(cases[t] | incidence[t], phi);
//     Reff[t]       = beta / gamma;
//   }
// }
// 




//      Instantaneous incidence
// //////////////////  44 days considers removed data points ///////////

// functions {
//   real[] seir_ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
//     real beta  = theta[1];
//     real alpha = theta[2];
//     real gamma = theta[3];
// 
//     real S = y[1];
//     real E = y[2];
//     real I = y[3];
//     real R = y[4];
// 
//     real N = S + E + I + R;
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
//   int<lower=1> n_days;              // number of observed days
//   int<lower=1> n_total_days;        // full 44 days
//   int<lower=1> burn_in_days;
//   int<lower=1> obs_days[n_days];    // indices of non-outlier days
//   real ts[n_total_days + burn_in_days];
//   real t0;
//   int cases[n_days];                // observed cases only
//   real<lower=0> N;
// }
// 
// transformed data {
//   real x_r[0];
//   int x_i[0];
// }
// 
// parameters {
//   real<lower=0> R_0;
//   real<lower=0> alpha_inv;
//   real<lower=0> gamma_inv;
//   real<lower=0> phi_inv;
//   real<lower=0, upper=1> p_reported;
//   real<lower=1, upper=1e5> E0;
//   real<lower=1, upper=1e5> I0;
// }
// 
// transformed parameters {
//   real alpha = 1.0 / alpha_inv;
//   real gamma = 1.0 / gamma_inv;
//   real beta  = R_0 * gamma;
// 
//   real S0 = N - E0 - I0;
//   real y[n_total_days + burn_in_days, 4];
//   real incidence[n_total_days];
//   real phi = 1.0 / phi_inv;
// 
//   if (S0 < 0)
//     S0 = 1e-6;
// 
//   {
//     real y0[4]    = { S0, E0, I0, 0 };
//     real theta[3] = { beta, alpha, gamma };
// 
//     y = integrate_ode_bdf(seir_ode, y0, t0, ts, theta, x_r, x_i);
//   }
// 
//   for (t in 1:n_total_days) {
//     incidence[t] = fmax(1e-6, p_reported * alpha * y[t + burn_in_days, 2]);
//   }
// }
// 
// model {
//   R_0 ~ lognormal(1.0, 0.2);
//   alpha_inv ~ lognormal(1.6, 0.13);
//   gamma_inv ~ lognormal(1.6, 0.13);
// 
//   I0 ~ lognormal(log(5), 0.8);
//   E0 ~ lognormal(log(10), 0.8);
// 
//   p_reported ~ beta(0.28, 100);
//   phi_inv ~ exponential(10);
// 
//   // likelihood only on non-outlier days
//   for (i in 1:n_days) {
//     cases[i] ~ neg_binomial_2(incidence[obs_days[i]], 1.0 / phi_inv);
//   }
// }
// 
// generated quantities {
//   real Reff[n_total_days];
//   int pred_cases[n_total_days];
//   
// 
//   for (t in 1:n_total_days) {
//     pred_cases[t] = neg_binomial_2_rng(incidence[t], 1.0 / phi_inv);
//     Reff[t]       = R_0;
//   }
// 
//   
// }







// Interval incidence
// //////////////////  44 days considers removed data points ///////////

functions {
  real[] seir_ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
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
  int<lower=1> obs_days[n_days];    
  real ts[n_total_days + burn_in_days];
  real t0;
  int cases[n_days];                
  real<lower=0> N;
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real<lower=0> R_0;
  real<lower=0> alpha_inv;
  real<lower=0> gamma_inv;
  real<lower=0> phi_inv;
  real<lower=0, upper=1> p_reported;
  real<lower=1, upper=1e5> E0;
  //real<lower=1, upper=1e5> I0;
}

transformed parameters {
  real alpha = 1.0 / alpha_inv;
  real gamma = 1.0 / gamma_inv;
  real beta  = R_0 * gamma;
  real I0 = 5;    // fixed value

  real S0 = N - E0 - I0;
  real y[n_total_days + burn_in_days, 4];

  // INTERVAL incidence: length = n_total_days - 1
  real incidence[n_total_days];

  real phi = 1.0 / phi_inv;

  if (S0 < 0)
    S0 = 1e-6;

  {
    real y0[4]    = { S0, E0, I0, 0 };
    real theta[3] = { beta, alpha, gamma };

    y = integrate_ode_bdf(seir_ode, y0, t0, ts, theta, x_r, x_i);
  }

  // Compute interval incidence: infections between t and t+1
  for (t in 1:(n_total_days)) {
    real E_t = y[t + burn_in_days, 2];
    incidence[t] = fmax(1e-6, p_reported * alpha * E_t);
  }
}

model {
  R_0 ~ lognormal(1.0, 0.2);
  alpha_inv ~ lognormal(1.6, 0.13);
  gamma_inv ~ lognormal(1.6, 0.13);

  //I0 ~ lognormal(log(5), 0.8);
  E0 ~ lognormal(log(10), 0.8);

  p_reported ~ beta(0.28, 100);
  phi_inv ~ exponential(10);

  // likelihood on interval incidence
  for (i in 1:n_days) {
    int t = obs_days[i];
    if (t > 1)
      cases[i] ~ neg_binomial_2(incidence[t], phi);
  }
}

 
generated quantities {
  real Reff[n_total_days];
  int pred_cases[n_total_days];

  for (t in 1:(n_total_days)) {
    pred_cases[t] = neg_binomial_2_rng(incidence[t], phi);
    Reff[t]       = R_0;
  }
}

