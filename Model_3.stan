// // /////////   42 days ///////////////////////////////////////////
// //   
//   functions {
//     real switch_eta(real t, real t1, real eta, real nu, real xi) {
//       return eta + (1 - eta) / (1 + exp(xi * (t - t1 - nu)));
//     }
//     real[] sir(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
//       int N = x_i[1];
//       real tswitch = x_r[1];
// 
//       real beta = theta[1];
//       real gamma = 1 / theta[2]; // gamma_inv
//       real a = 1 / theta[3]; // a_inv
//       real eta = theta[4];
//       real nu = theta[5];
//       real xi = theta[6];
//       real i0 = theta[7];
//       real e0 = theta[8];
//       real omega = theta[9]; // prior immunity
//       real forcing_function = switch_eta(t, tswitch, eta, nu, xi);
//       real beta_eff = beta * forcing_function;
// 
//       real init[4] = {(1 - omega) * N - i0 - e0, e0, i0, omega * N};
//       real S = y[1] + init[1];
//       real E = y[2] + init[2];
//       real I = y[3] + init[3];
//       real R = y[4] + init[4];
// 
//       real dS_dt = -beta_eff * I * S / N;
//       real dE_dt =  beta_eff * I * S / N - a * E;
//       real dI_dt = a * E - gamma * I;
//       real dR_dt = gamma * I;
// 
//       return {dS_dt, dE_dt, dI_dt, dR_dt};
//     }
//   }
// data {
//   int<lower=1> n_days;
//   real t0;
//   real tswitch;
//   real ts[n_days];
//   int N;
//   int cases[n_days];
// }
// transformed data {
//   int x_i[1] = { N };
//   real x_r[1] = {tswitch};
// }
// parameters {
//   real<lower=0> beta;
//   real<lower=0> a_inv;
//   real<lower=0> gamma_inv;
//   real<lower=0> i0;
//   real<lower=0> e0;
//   real<lower=0, upper=1> p_reported;
//   real<lower=0> phi_inv;
//   real<lower=0,upper=1> eta;
//   real<lower=0> nu;
//   real<lower=0,upper=1> xi_raw;
//   real<lower=0,upper=1> omega;
// }
// transformed parameters {
//   real y[n_days, 4];
//   real incidence[n_days - 1];
//   real phi = 1. / phi_inv;
//   real xi = xi_raw + 0.5;
//   real theta[9] = {beta, gamma_inv, a_inv, eta, nu, xi, i0, e0, omega};
// 
//   y = integrate_ode_rk45(sir, rep_array(0.0, 4), t0, ts, theta, x_r, x_i);
//   for (i in 1:n_days-1) {
//     incidence[i] = -(y[i+1, 2] - y[i, 2] + y[i+1, 1] - y[i, 1]) * p_reported;
// 
//   }
// }
// model {
//   beta ~ lognormal(-0.6094, 0.2); //gamma(7, 20);
//   a_inv ~ lognormal(1.6, 0.13); //lognormal(1.5, 0.5);  n
//   gamma_inv ~ lognormal(1.6, 0.13); //lognormal(1.5, 0.5);
//   i0 ~ lognormal(log(5), 0.8); //normal(50, 5);
//   e0 ~ lognormal(log(10), 0.8); //normal(100, 5);
//   p_reported ~ beta(0.28, 100); //beta(2, 38);
//   phi_inv ~ exponential(10); //exponential(5);
//   eta ~ beta(2.5, 4);
//   nu ~ exponential(1. / 5); //exponential(1. / 5);
//   xi_raw ~ beta(2, 5);
//   omega ~ beta(2, 5); // prior on immunity
// 
//   cases[1:(n_days-1)] ~ neg_binomial_2(incidence, phi);
// }
// generated quantities {
//   real gamma = 1 / gamma_inv;
//   real a = 1 / a_inv;
//   real R0 = beta / gamma;
//   real Reff[n_days];
//   real R_infty;
//   real recovery_time = gamma_inv;
//   real incubation_time = a_inv;
//   real pred_cases[n_days-1];
// 
//   pred_cases = neg_binomial_2_rng(incidence, phi);
//   for (i in 1:n_days)
//     Reff[i] = switch_eta(i, tswitch, eta, nu, xi) * beta / gamma;
// 
//   R_infty = sum(cases) / (p_reported * (1 - omega) * N);  //final attack rate approximation
//   //R_infty = (1 - omega) * (1 - exp(-R0)); // final attack rate approximation
// }
// 











// // /////////   44 days ///////////////////////////////////////////
// // incidence[i] = -(y[i+1, 2] - y[i, 2] + y[i+1, 1] - y[i, 1]) * p_reported;  
// 
// functions {
//   real switch_eta(real t, real t1, real eta, real nu, real xi) {
//     return eta + (1 - eta) / (1 + exp(xi * (t - t1 - nu)));
//   }
// 
//   real[] sir(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
//     int N = x_i[1];
//     real tswitch = x_r[1];
// 
//     real beta = theta[1];
//     real gamma = 1 / theta[2];
//     real a = 1 / theta[3];
//     real eta = theta[4];
//     real nu = theta[5];
//     real xi = theta[6];
//     real i0 = theta[7];
//     real e0 = theta[8];
//     real omega = theta[9];
// 
//     real forcing_function = switch_eta(t, tswitch, eta, nu, xi);
//     real beta_eff = beta * forcing_function;
// 
//     // Initial conditions added inside ODE
//     real init[4] = {(1 - omega) * N - i0 - e0, e0, i0, omega * N};
//     real S = y[1] + init[1];
//     real E = y[2] + init[2];
//     real I = y[3] + init[3];
//     real R = y[4] + init[4];
// 
//     real dS_dt = -beta_eff * I * S / N;
//     real dE_dt =  beta_eff * I * S / N - a * E;
//     real dI_dt = a * E - gamma * I;
//     real dR_dt = gamma * I;
// 
//     return {dS_dt, dE_dt, dI_dt, dR_dt};
//   }
// }
// 
// data {
//   int<lower=1> n_days;       
//   int<lower=1> n_obs;        
//   int<lower=1> obs_idx[n_obs]; 
//   real t0;
//   real tswitch;
//   real ts[n_days];
//   int N;
//   int cases[n_days];
// }
// 
// transformed data {
//   int x_i[1] = { N };
//   real x_r[1] = { tswitch };
// }
// 
// parameters {
//   real<lower=0> beta;
//   real<lower=0> a_inv;
//   real<lower=0> gamma_inv;
//   //real<lower=0> i0;
//   real<lower=0> e0;
//   real<lower=0, upper=1> p_reported;
//   real<lower=1e-6> phi_inv;
//   real<lower=0,upper=1> eta;
//   real<lower=0> nu;
//   real<lower=0,upper=1> xi_raw;
//   real<lower=0,upper=1> omega;
// }
// 
// transformed parameters {
//   real y[n_days, 4];
//   real incidence[n_days];
//   real phi = 1. / phi_inv;
//   real i0 = 5;    // fixed value
//   real xi = xi_raw + 0.5;
//   real theta[9] = {beta, gamma_inv, a_inv, eta, nu, xi, i0, e0, omega};
// 
//   y = integrate_ode_rk45(sir, rep_array(0.0, 4), t0, ts, theta, x_r, x_i);
// 
//   for (i in 1:(n_days)) {
//     real raw_inc = -(y[i+1, 2] - y[i, 2] + y[i+1, 1] - y[i, 1]) * p_reported;
//     incidence[i] = fmax(raw_inc, 1e-9);
//   }
// }
// 
// model {
//   // Hard constraint to avoid impossible initial conditions
//   if (i0 + e0 >= (1 - omega) * N)
//     target += negative_infinity();
// 
//   beta ~ lognormal(-0.6094, 0.2);
//   a_inv ~ lognormal(1.6, 0.13);
//   gamma_inv ~ lognormal(1.6, 0.13);
//   //i0 ~ lognormal(log(5), 0.8);
//   e0 ~ lognormal(log(10), 0.8);
//   p_reported ~ beta(0.28, 100);
//   phi_inv ~ exponential(10);
//   eta ~ beta(2.5, 4);
//   nu ~ exponential(1. / 5);
//   xi_raw ~ beta(2, 5);
//   omega ~ beta(2, 5);
// 
//   // Likelihood only for non-outlier days (and only where incidence is defined)
//   for (j in 1:n_obs) {
//     int k = obs_idx[j];          // k ∈ {1, ..., n_days - 1}
//     cases[k] ~ neg_binomial_2(incidence[k], phi);
//   }
// }
// 
// generated quantities {
//   real gamma = 1 / gamma_inv;
//   real a = 1 / a_inv;
//   real R0 = beta / gamma;
//   real Reff[n_days];
//   real R_infty;
//   real recovery_time = gamma_inv;
//   real incubation_time = a_inv;
//   real pred_cases[n_days];
// 
//   pred_cases = neg_binomial_2_rng(incidence, phi);
// 
//   for (i in 1:n_days)
//     Reff[i] = switch_eta(i, tswitch, eta, nu, xi) * beta / gamma;
// 
//   R_infty = sum(cases) / (p_reported * (1 - omega) * N);
// }




// /////////   44 days ///////////////////////////////////////////
// incidence[i] = a * E * p_reported;  
// 
functions {
  real switch_eta(real t, real t1, real eta, real nu, real xi) {
    return eta + (1 - eta) / (1 + exp(xi * (t - t1 - nu)));
  }

  real[] sir(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    int N = x_i[1];
    real tswitch = x_r[1];

    real beta = theta[1];
    real gamma = 1 / theta[2];
    real a = 1 / theta[3];
    real eta = theta[4];
    real nu = theta[5];
    real xi = theta[6];
    real i0 = theta[7];
    real e0 = theta[8];
    real omega = theta[9];

    real forcing_function = switch_eta(t, tswitch, eta, nu, xi);
    real beta_eff = beta * forcing_function;

    // Initial conditions embedded inside ODE
    real init[4] = {(1 - omega) * N - i0 - e0, e0, i0, omega * N};
    real S = y[1] + init[1];
    real E = y[2] + init[2];
    real I = y[3] + init[3];
    real R = y[4] + init[4];

    real dS_dt = -beta_eff * I * S / N;
    real dE_dt =  beta_eff * I * S / N - a * E;
    real dI_dt = a * E - gamma * I;
    real dR_dt = gamma * I;

    return {dS_dt, dE_dt, dI_dt, dR_dt};
  }
}

data {
  int<lower=1> n_days;
  int<lower=1> n_obs;
  int<lower=1> obs_idx[n_obs];
  real t0;
  real tswitch;
  real ts[n_days];
  int N;
  int cases[n_days];
}

transformed data {
  int x_i[1] = { N };
  real x_r[1] = { tswitch };
}

parameters {
  real<lower=0> beta;
  real<lower=0> a_inv;
  real<lower=0> gamma_inv;
  real<lower=0> e0;
  real<lower=0, upper=1> p_reported;
  real<lower=1e-6> phi_inv;
  real<lower=0,upper=1> eta;
  real<lower=0> nu;
  real<lower=0,upper=1> xi_raw;
  real<lower=0,upper=1> omega;
}

transformed parameters {
  real y[n_days, 4];
  real incidence[n_days];
  real phi = 1.0 / phi_inv;
  real i0 = 5;                     // fixed
  real xi = xi_raw + 0.5;
  real a = 1.0 / a_inv;
  real theta[9] = {beta, gamma_inv, a_inv, eta, nu, xi, i0, e0, omega};

  y = integrate_ode_rk45(sir, rep_array(0.0, 4), t0, ts, theta, x_r, x_i);

  // NEW: incidence defined for ALL n_days
  for (i in 1:n_days) {
    real E_i = y[i, 2] + e0;       // E(t_i)
    incidence[i] = fmax(a * E_i * p_reported, 1e-9);
  }
}

model {
  if (i0 + e0 >= (1 - omega) * N)
    target += negative_infinity();

  beta ~ lognormal(-0.6094, 0.2);
  a_inv ~ lognormal(1.6, 0.13);
  gamma_inv ~ lognormal(1.6, 0.13);
  e0 ~ lognormal(log(10), 0.8);
  p_reported ~ beta(0.28, 100);
  phi_inv ~ exponential(10);
  eta ~ beta(2.5, 4);
  nu ~ exponential(1.0 / 5);
  xi_raw ~ beta(2, 5);
  omega ~ beta(2, 5);

  // Likelihood for non-outlier days
  for (j in 1:n_obs) {
    int k = obs_idx[j];
    cases[k] ~ neg_binomial_2(incidence[k], phi);
  }
}

generated quantities {
  real gamma = 1 / gamma_inv;
  //real a = 1 / a_inv;
  real R0 = beta / gamma;
  real Reff[n_days];
  real R_infty;
  real recovery_time = gamma_inv;
  real incubation_time = a_inv;
  real pred_cases[n_days];

  for (i in 1:n_days)
    pred_cases[i] = neg_binomial_2_rng(incidence[i], phi);

  for (i in 1:n_days)
    Reff[i] = switch_eta(i, tswitch, eta, nu, xi) * beta / gamma;

  R_infty = sum(cases) / (p_reported * (1 - omega) * N);
}


















