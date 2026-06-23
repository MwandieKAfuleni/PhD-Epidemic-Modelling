################################################################################

# Model 1

# Full SEIR

################################################################################

library(tidyverse)
library(ggplot2)
library(deSolve)
library(tidyr)
library(lubridate)

library(rstan)

set.seed(42)
# Load Data
file_path <- file.path("data", "MalawiBetaCases.csv")

BetaVariant <- read.csv(file_path) %>%
  mutate(date = lubridate::dmy(date))   

BetaVariant <- BetaVariant%>%
  filter(date >= as.Date("2021-01-02") &
           date <= as.Date("2021-02-14"))



cases <- BetaVariant$new_cases    
N <- 17563749   

all_days <- length(cases)        

obs_days <- setdiff(1:all_days, c(11, 14))

cases_obs <- cases[obs_days]
n_days <- length(cases_obs)       

burn_in_days <- 32
ts <- seq_len(all_days + burn_in_days)

stan_data_1 <- list(
  n_days = n_days,
  n_total_days = all_days,
  burn_in_days = burn_in_days,
  obs_days = obs_days,
  t0 = 0,
  ts = ts,
  cases = cases_obs,
  N = N
)



#==============================================
# Compile and fit
#==============================================
model1 <- stan_model("model", "Model_1.stan")
 
seir_fit_1 <- sampling(model1,
                  data = stan_data_1,
                  chains = 4,
                  cores = 4,
                  iter = 4000,
                  warmup = 2000,
                  control = list(adapt_delta = 0.9999, max_treedepth = 20),
                  seed = 123)

# Save the results
saveRDS(seir_fit_1, file = "FinalModel_1.rds")

# Load the results
#seir_fit_1 <- readRDS(------------------/FinalModel_1.RDS")


#==============================================
# Calculate \beta  
#==============================================
post <- extract(seir_fit_1)
beta <- post$R_0 * post$gamma    
mean_beta <- mean(beta)
crI_beta <- quantile(beta, c(0.025, 0.975))
mean_beta
crI_beta


#==============================================
# Calculate \pi to 4 decimal places
#==============================================
post <- rstan::extract(seir_fit_1)
p <- post$p_reported
mean_p <- mean(p)
crI_p <- quantile(p, c(0.025, 0.975))
sprintf("%.4f (%.4f, %.4f)",
        mean_p, crI_p[1], crI_p[2])


#==============================================
# Traceplot
#==============================================
stan_trace(seir_fit_1, pars = c("R_0", "alpha_inv", "gamma_inv", "E0", "p_reported", "phi_inv")) 


#==============================================
# Posterior summary
#==============================================
print(seir_fit_1, pars = c("R_0", "alpha_inv", "gamma_inv", "E0", "p_reported", "phi_inv"))


#==============================================
# Density plot
#==============================================
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

draws <- as_draws_df(seir_fit_1)

draws_long <- draws %>%
  select(.chain, R_0, alpha_inv, gamma_inv, E0, p_reported, phi_inv) %>%
  pivot_longer(
    cols = c(R_0, alpha_inv, gamma_inv, E0, p_reported, phi_inv),
    names_to = "parameter",
    values_to = "value"
  )

ggplot(draws_long,
       aes(x = value, colour = factor(.chain))) +
  geom_density(linewidth = 1.5) +
  facet_wrap(
    ~ parameter,
    scales = "free",
    labeller = as_labeller(
      c(
        R_0 = "R[0]",
        alpha_inv = "T[L]",
        gamma_inv = "T[I]",
        #I0 = "I[0]",
        E0 = "E[0]",
        p_reported = "pi",
        phi_inv = "phi^{-1}"
      ),
      label_parsed
    )
  ) +
  labs(colour = "Chain", fill = "Chain") +
  theme_bw()


#==============================================
# Pairwiseplot
#==============================================
pairs(seir_fit_1,
      pars = c("R_0", "alpha_inv", "gamma_inv", "phi_inv", "p_reported", "I0", "E0"),
      labels = c(
        expression(R[0]),
        expression(T[L]),
        expression(T[I]),
        #expression(I[0]),
        expression(E[0]),
        expression(pi),
        expression(phi^{-1})
      )
)$summary


#==============================================
## Posterior predictive check
#==============================================

outliers <- c(11, 14)

cases <- BetaVariant$new_cases    

 
obs_days <- setdiff(1:length(cases), outliers)

cases_obs <- cases[obs_days]

 
pred_summary <- as.data.frame(
  summary(seir_fit_1, pars = "pred_cases")$summary
)

 
t_pred <- 1:(length(cases))

 
interval_idx <- obs_days
pred_summary_obs <- pred_summary[interval_idx, ]

Pred1 <- cbind(
  pred_summary_obs,
  day = obs_days,       
  cases = cases_obs
)

colnames(Pred1) <- make.names(colnames(Pred1))

ggplot(Pred1, aes(x = day)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "darkred", alpha = 0.35) +
  geom_line(aes(y = X50., color = "Median prediction")) +
  geom_point(aes(y = cases, color = "Observed cases")) +
  scale_color_manual(values = c("Median prediction" = "red", "Observed cases" = "black")) +
  labs(
    #title = "Posterior Predictive Check",
    x = "Day",
    y = "Incidence (interval)",
    color = "Legend"
  ) +
  theme_minimal()

 

#==============================================
# Prior versus posterior distribution
#==============================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

# Extract posterior samples
post <- as.data.frame(
  rstan::extract(
    seir_fit_1,
    pars = c(
      "R_0",
      "alpha_inv",
      "gamma_inv",
      "phi_inv",
      "p_reported",
      #"I0",
      "E0"
    )
  )
)


# Prior density functions
priors <- list(
  R_0        = function(x) dlnorm(x, 1.009619, 0.2),
  alpha_inv      = function(x) dlnorm(x, log(4.5), 0.5),
  gamma_inv      = function(x) dlnorm(x, log(5), 0.5),
  phi_inv    = function(x) dexp(x, 10),
  p_reported = function(x) dbeta(x, 0.28, 100),
  #E0         = function(x) dlnorm(x, log(10), 0.5)
  E0         = function(x) dlnorm(10^x, log(10), 0.5)
)

# Facet labels
param_labels <- c(
  R_0        = "R[0]",
  alpha_inv      = "T[L]",
  gamma_inv      = "T[I]",
  phi_inv    = "phi^{-1}",
  p_reported = "pi",
  E0         = "E[0]"
)

# Prior centres for vertical lines
prior_centres <- data.frame(
  
  param = factor(
    unname(param_labels),
    levels = unname(param_labels)
  ),
  
  xintercept = c(
    exp(1.009619 + 0.2^2/2),   #R_0
    exp(log(5) + 0.13^2/2),   # T_L     
    exp(log(5) + 0.13^2/2),  #T_I    
    1/10,         # phi_inv    
    0.0028,       # p_reported
    log10(exp(log(10) + 0.8^2/2)) # E0       
  )
)

# Build plotting data
plot_data <- lapply(names(post), function(p) {
  #posterior <- density(post[[p]])
  posterior <- if (p == "E0") {
    density(log10(post[[p]]))
  } else {
    density(post[[p]])
  }
  prior_x <- switch(p,
    p_reported = seq(0, 0.005, length.out = 300),
    phi_inv    = seq(0, 0.5, length.out = 300),
    alpha_inv      = seq(0, 12, length.out = 300),
    gamma_inv      = seq(0, 6, length.out = 300),
    R_0        = seq(0, 6, length.out = 300),
    #E0         = seq(0, 300, length.out = 300),
    E0         = seq(log10(1), log10(4000), length.out = 300),
    
    seq(
      min(posterior$x),
      max(posterior$x),
      length.out = 300
    )
  )
  
  #prior_y <- priors[[p]](prior_x)
  prior_y <- if (p == "E0"){
    dlnorm(10^prior_x, meanlog=log(10), sdlog = 0.8)
  } else {
    priors[[p]](prior_x)
  }
  
  data.frame(
    
    param = rep(unname(param_labels[p]),
                length(c(prior_x, posterior$x))),
    
    #x = c(prior_x, posterior$x),
    x = c(
      if (p == "E0") prior_x else prior_x,
      if (p == "E0") posterior$x else posterior$x
    ),
    
    density = c(
      prior_y / max(prior_y),
      posterior$y / max(posterior$y)
    ),
    
    type = c(
      rep("Prior", length(prior_x)),
      rep("Posterior", length(posterior$x))
    )
  )
  
}) %>% bind_rows()

# Factor ordering
plot_data$param <- factor(plot_data$param,levels = unname(param_labels))

# Plot
ggplot(plot_data, aes(x = x, y = density, colour = type)) +
  
  # prior centre line
  geom_vline(
    data = prior_centres,
    aes(xintercept = xintercept),
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +
  
# density curves
geom_line(linewidth = 1.1) +
  
  # facets
  facet_wrap(
    ~param,
    scales = "free",
    ncol = 2,
    labeller = label_parsed
  ) +
  
  labs(
    x = "Parameter value",
    y = "Normalised density",
    colour = NULL,
    #title = "Prior vs Posterior Distributions"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    legend.position = "top",
    strip.text = element_text(size = 13),
    plot.title = element_text(face = "bold")
  )






























################################################################################

# Model 2

# SEIR with Eigenvector informed initialisation

################################################################################

library(tidyverse)
library(ggplot2)
library(deSolve)
library(tidyr)
library(lubridate)
library(rstan)
 

set.seed(42)

# Load Data
file_path <- file.path("data","MalawiBetaCases.csv")

BetaVariant <- read.csv(file_path) %>%
  mutate(date = lubridate::dmy(date))

BetaVariant <- BetaVariant%>%
  filter(date >= as.Date("2021-01-02") &
           date <= as.Date("2021-02-14"))

cases <- BetaVariant$new_cases
N <- 17563749
n_days <- length(cases)

# Fixed epidemiological parameters  
alpha_fix <- 1 / 5
gamma_fix <- 1 / 5

burn_in_days <- 120      
t0 <- 0
ts <- seq_len(n_days + burn_in_days)

outliers <- c(11, 14)

obs_idx <- setdiff(1:n_days, outliers)
n_obs <- length(obs_idx)


stan_data_2 <- list(
  n_days       = n_days,
  burn_in_days = burn_in_days,
  t0           = t0,
  ts           = ts,
  cases        = cases,
  N            = N,
  alpha        = alpha_fix,
  gamma        = gamma_fix,
  n_obs        = n_obs,
  obs_idx      = obs_idx
)

 
#==============================================
# Compile and fit
#==============================================
model2 <- stan_model("model", "Model_2.stan")

seir_fit_2 <- sampling(
  model2,
  data = stan_data_2,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.9999, max_treedepth = 20),
  seed = 123
)


# Save the results
saveRDS(seir_fit_2, file = "FinalModel_2.rds")

# Load the results
#seir_fit_2 <- readRDS("------------------/FinalModel_2.RDS")

#==============================================
# calculate \beta
#==============================================
post <- extract(seir_fit_2)
gamma <- 1 / post$gamma_inv
beta <- post$R_0 * gamma 
mean_beta <- mean(beta)
crI_beta <- quantile(beta, c(0.025, 0.975))
mean_beta
crI_beta


#==============================================
# Calculate \pi 
#==============================================
post <- rstan::extract(seir_fit_2)
p <- post$p_reported
mean_p <- mean(p)
crI_p <- quantile(p, c(0.025, 0.975))
sprintf("%.4f (%.4f, %.4f)",
        mean_p, crI_p[1], crI_p[2])





#==============================================
# Traceplot
#==============================================
stan_trace(seir_fit_2, pars = c("R_0","I0", "p_reported", "phi_inv")) 

#==============================================
# Summary estimates
#==============================================
print(seir_fit_2, pars = c("R_0","I0", "p_reported", "phi_inv"))


#==============================================
# Density plot
#==============================================
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

draws <- as_draws_df(seir_fit_2)

draws_long <- draws %>%
  select(.chain, R_0, I0, p_reported, phi_inv) %>%
  pivot_longer(
    cols = c(R_0, I0, p_reported, phi_inv),
    names_to = "parameter",
    values_to = "value"
  )

ggplot(draws_long,
       aes(x = value, colour = factor(.chain))) +
  geom_density(linewidth = 1.5) +
  facet_wrap(
    ~ parameter,
    scales = "free",
    labeller = as_labeller(
      c(
        R_0 = "R[0]",
        I0 = "I[0]",
        p_reported = "pi",
        phi_inv = "phi^{-1}"
      ),
      label_parsed
    )
  ) +
  labs(colour = "Chain", fill = "Chain") +
  theme_bw()




#==============================================
# Pairwise plot
#==============================================
pairs(seir_fit_2,
      pars = c("R_0","I0", "p_reported", "phi_inv"),
      labels = c(
        expression(R[0]),
        expression(I[0]),
        expression(pi),
        expression(phi^{-1})
        
        
      )
)$summary
 


#==============================================
# PPC
#==============================================
outliers <- c(11, 14)

cases_plot <- cases
cases_plot[outliers] <- NA

pred_summary <- as.data.frame(summary(seir_fit_2, pars = "pred_cases")$summary)

Pred2 <- cbind(
  pred_summary,
  t = 1:n_days,          # correct time index for your Stan model
  cases = cases_plot
)

colnames(Pred2) <- make.names(colnames(Pred2))

ggplot(Pred2, aes(x = t)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.),
              fill = "darkred", alpha = 0.35) +
  geom_line(aes(y = X50., color = "Median prediction")) +
  geom_point(aes(y = cases, color = "Observed cases")) +
  scale_color_manual(values = c(
    "Median prediction" = "red",
    "Observed cases" = "black"
  )) +
  labs(
    #title = "Posterior Predictive Check",
    x = "Day",
    y = "Incidence",
    color = "Legend"
  ) +
  theme_minimal()


#==============================================
# Prior vs posterior plot
#==============================================
 
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

# Extract posterior samples
post <- as.data.frame(
  rstan::extract(
    seir_fit_2,
    pars = c(
      "R_0",
      "I0",
      "p_reported",
      "phi_inv"
      
    )
  )
)


priors <- list(
  R_0        = function(x) dlnorm(x, meanlog = 1.0, sdlog = 0.2),
  I0         = function(x) dlnorm(x, meanlog = log(5), sdlog = 0.8),
  p_reported = function(x) dbeta(x, shape1 = 0.28, shape2 = 100),
  phi_inv    = function(x) dexp(x, rate = 10)
)

# ---------------------------------------------------
# Facet labels
# ---------------------------------------------------

param_labels <- c(
  R_0        = "R[0]",
  I0      = "I[0]",
  p_reported = "pi",
  phi_inv    = "phi^{-1}"
)

# ---------------------------------------------------
# Prior centres for vertical lines
# ---------------------------------------------------

prior_centres <- data.frame(
  param = factor(
    unname(param_labels),
    levels = unname(param_labels)
  ),
  xintercept = c(
    exp(1 + 0.2^2/2),          # R_0 mean ≈ 2.773
    exp(log(5) + 0.8^2/2),     # I0 mean ≈ 6.885
    0.28 / (0.28 + 100),       # p_reported mean ≈ 0.00279
    1/10                       # phi_inv mean = 0.1
  )
)
# ---------------------------------------------------
# Build plotting data
# ---------------------------------------------------

plot_data <- lapply(names(post), function(p) {
  
  posterior <- density(post[[p]])
  
  prior_x <- switch(
    p,
    
    R_0        = seq(0, 4, length.out = 300),
    I0         = seq(0, 20, length.out = 300),
    p_reported = seq(0, 0.005, length.out = 300),
    phi_inv    = seq(0, 0.5, length.out = 300),
    
    seq(
      min(posterior$x),
      max(posterior$x),
      length.out = 300
    )
  )
  
  prior_y <- priors[[p]](prior_x)
  
  data.frame(
    
    param = rep(unname(param_labels[p]),
                length(c(prior_x, posterior$x))),
    
    x = c(prior_x, posterior$x),
    
    density = c(
      prior_y / max(prior_y),
      posterior$y / max(posterior$y)
    ),
    
    type = c(
      rep("Prior", length(prior_x)),
      rep("Posterior", length(posterior$x))
    )
  )
  
}) %>% bind_rows()

# ---------------------------------------------------
# Factor ordering
# ---------------------------------------------------

plot_data$param <- factor(
  plot_data$param,
  levels = unname(param_labels)
)

# ---------------------------------------------------
# Plot
# ---------------------------------------------------

ggplot(plot_data,
       aes(x = x,
           y = density,
           colour = type)) +
  
  # prior centre line
  geom_vline(
    data = prior_centres,
    aes(xintercept = xintercept),
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +
  
  # density curves
  geom_line(linewidth = 1.1) +
  
  # facets
  facet_wrap(
    ~param,
    scales = "free",
    ncol = 2,
    labeller = label_parsed
  ) +
  
  labs(
    x = "Parameter value",
    y = "Normalised density",
    colour = NULL
    #title = "Prior vs Posterior Distributions"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    legend.position = "top",
    strip.text = element_text(size = 13),
    plot.title = element_text(face = "bold")
  )






























################################################################################

# Model 3

# SEIR Model - Transmission modelled with a logistic function

################################################################################

library(tidyverse)
library(ggplot2)
library(deSolve)
library(tidyr)
library(lubridate)
library(rstan)

 
file_path <- file.path("data" , "MalawiBetaCasesDeaths.csv")  
BetaVariant <- read.csv(file_path) %>% mutate(date = dmy(date))
BetaVariant <- BetaVariant%>%
  filter(date >= as.Date("2021-01-02") &
           date <= as.Date("2021-02-14"))

cases <- BetaVariant$new_cases
N <- 17563749  # Total population size
n_days <- length(cases)
 

date_switch <- "2021-01-17"
tswitch <- BetaVariant %>% filter(date < date_switch) %>% nrow() + 1

outliers <- c(11, 14)

obs_idx <- setdiff(1:(n_days), outliers)
n_obs <- length(obs_idx)

ts <- 1:n_days

data_seir3 <- list(
  n_days = n_days,
  n_obs = n_obs,
  obs_idx = obs_idx,
  t0 = 0,
  ts = ts,
  cases = cases,
  N = N,
  tswitch = tswitch
)

#==============================================
## Compile and fit 
#==============================================
 model3 <- stan_model("model", "Model_3.stan")

seir_fit_3 <- sampling(model3, 
                       data = data_seir3,  
                       iter = 4000, 
                       warmup = 2000, 
                       chains = 4, 
                       cores = 4,
                       seed = 420, 
                       control = list(adapt_delta = 0.999, max_treedepth = 20)
)

#  Save the results
saveRDS(seir_fit_3, file = "FinalModel_3.rds")

# Load the results
#seir_fit_3 <- readRDS("------------------/FinalModel_3.RDS")




#==============================================
# Traceplot
#==============================================
stan_trace(seir_fit_3, pars = c("beta", "a_inv", "gamma_inv", "e0", 
                                "p_reported", "phi_inv", "eta", "nu", 
                                "xi_raw", "omega"))

#==============================================
# Summary estimates
#==============================================
summary(seir_fit_3, pars = c("beta", "a_inv", "gamma_inv", "e0", 
                             "p_reported", "phi_inv", "eta", "nu", 
                             "xi_raw", "omega", "R_infty", "R0"))$summary



#==============================================
# density plot
#==============================================
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

draws <- as_draws_df(seir_fit_3)

draws_long <- draws %>%
  select(.chain, beta, a_inv, gamma_inv, e0, p_reported, 
         phi_inv, eta, nu, xi_raw, omega) %>%
  pivot_longer(
    cols = c(beta, a_inv, gamma_inv, e0, p_reported, 
             phi_inv, eta, nu, xi_raw, omega),
    names_to = "parameter",
    values_to = "value"
  )

ggplot(draws_long,
       aes(x = value, colour = factor(.chain))) +
  geom_density(linewidth = 1.5) +
  facet_wrap(
    ~ parameter,
    scales = "free",
    labeller = as_labeller(
      c(
        beta = "beta",
        a_inv = "T[L]",
        gamma_inv = "T[I]",
        #i0 = "I[0]",
        e0 = "E[0]",
        p_reported = "pi",
        phi_inv = "phi^{-1}",
        eta = "eta", 
        nu = "nu",
        xi_raw = "xi_raw", 
        omega = "omega"
        
      ),
      label_parsed
    )
  ) +
  labs(colour = "Chain", fill = "Chain") +
  theme_bw()




#==============================================
#Pairwise plot
#==============================================
pairs(seir_fit_3,
      pars = c("beta", "a_inv", "gamma_inv", "e0",
               "p_reported", "phi_inv", "eta", "nu",
               "xi_raw", "omega"),
      labels = c(
        expression(beta),
        expression(T[L]),
        expression(T[I]),
        #expression(I[0]),
        expression(E[0]),
        expression(pi),
        expression(phi^{-1}),
        expression(eta),
        expression(nu),
        expression(xi[raw]),
        expression(omega)
      )
)$summary


#==============================================
# PPC
#==============================================
outliers <- c(11, 14)

cases_plot <- cases
cases_plot[outliers] <- NA

pred_summary <- as.data.frame(summary(seir_fit_3, pars = "pred_cases")$summary)

Pred3 <- cbind(
  pred_summary,
  t = ts,            
  cases = cases_plot
)

colnames(Pred3) <- make.names(colnames(Pred3))

ggplot(Pred3, aes(x = t)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "darkred", alpha = 0.35) +
  geom_line(aes(y = X50., color = "Median prediction")) +
  geom_point(aes(y = cases, color = "Observed cases")) +
  scale_color_manual(values = c("Median prediction" = "red", "Observed cases" = "black")) +
  labs(
    #title = "Posterior Predictive Check",
    x = "Day",
    y = "Incidence",
    color = "Legend"
  ) +
  theme_minimal()



#============================================== 
# Prior vs posterior plot
#==============================================
 
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

 
post <- as.data.frame(
  rstan::extract(
    seir_fit_3,
    pars = c(
      "beta",
      "a_inv",
      "gamma_inv",
      "e0",
      "p_reported",
      "phi_inv",
      "eta",
      "nu",
      "xi_raw",
      "omega"
      
    )
  )
)


# ---------------------------------------------------
# Prior density functions
# ---------------------------------------------------

priors <- list(
  
  beta        = function(x) dlnorm(x, meanlog = -0.6094, sdlog = 0.2),
  a_inv       = function(x) dlnorm(x, meanlog = 1.6, sdlog = 0.13),
  gamma_inv   = function(x) dlnorm(x, meanlog = 1.6, sdlog = 0.13),
  #e0          = function(x) dlnorm(x, meanlog = log(10), sdlog = 0.8),
  e0          = function(x) dlnorm(10^x, meanlog = log(10), sdlog = 0.8),
  p_reported  = function(x) dbeta(x, shape1 = 0.28, shape2 = 100),
  phi_inv     = function(x) dexp(x, rate = 10),
  eta         = function(x) dbeta(x, shape1 = 2.5, shape2 = 4),
  nu          = function(x) dexp(x, rate = 1/5),
  xi_raw      = function(x) dbeta(x, shape1 = 2, shape2 = 5),
  omega       = function(x) dbeta(x, shape1 = 2, shape2 = 5)
)

# ---------------------------------------------------
# Facet labels
# ---------------------------------------------------

param_labels <- c(
   
  beta       = "beta",    #expression(beta),
  a_inv      = "T[L]",    #expression(T[L]),         
  gamma_inv  = "T[I]",    #expression(T[I]),         
  e0         = "E[0]",      #expression(E0),
  p_reported = "pi",      #expression(pi),
  phi_inv    = "phi^{-1}",#expression(phi^{-1}),
  eta        = "eta",     #expression(eta),
  nu         = "nu",      #expression(nu),
  xi_raw     = "xi[raw]", #expression(xi[raw]),
  omega      = "omega"    #expression(omega)
)

# ---------------------------------------------------
# Prior centres for vertical lines (means)
# ---------------------------------------------------

prior_centres <- data.frame(
  
  param = factor(
    unname(param_labels),
    levels = unname(param_labels)
  ),
  
  xintercept = c(
    
    exp(-0.6094 + 0.2^2/2),            # beta mean
    exp(1.6 + 0.13^2/2),               # a_inv mean
    exp(1.6 + 0.13^2/2),               # gamma_inv mean
    #exp(log(10) + 0.8^2/2),            # e0 mean
    log10(exp(log(10) + 0.8^2/2)),            # e0 mean
    0.28/(0.28 + 100),                 # p_reported mean
    1/10,                              # phi_inv mean
    2.5/(2.5 + 4),                     # eta mean
    5,                                 # nu mean
    2/(2 + 5),                         # xi_raw mean
    2/(2 + 5)                          # omega mean
  )
)
# ---------------------------------------------------
# Build plotting data
# ---------------------------------------------------

plot_data <- lapply(names(post), function(p) {
  
  #posterior <- density(post[[p]])
  posterior <- if (p == "e0"){
    density(log10(post[[p]]))
  } else {
    density(post[[p]])
  }
  
  prior_x <- switch(
    p,
    
    beta        = seq(0, 1, length.out = 300),
    a_inv       = seq(0, 8, length.out = 300),
    gamma_inv   = seq(0, 8, length.out = 300),
    #e0          = seq(0, 30, length.out = 300),
    e0          = seq(log10(1), log10(10000), length.out = 300),
    p_reported  = seq(0, 0.15, length.out = 300),
    phi_inv     = seq(0, 1, length.out = 300),
    eta         = seq(0, 0.75, length.out = 300),
    nu          = seq(0, 7, length.out = 300),
    xi_raw      = seq(0, 0.5, length.out = 300),
    omega       = seq(0, 0.5, length.out = 300),
    
    seq(
      min(posterior$x),
      max(posterior$x),
      length.out = 300
    )
  )
  
  prior_y <- priors[[p]](prior_x)
  
  data.frame(
    
    param = rep(unname(param_labels[p]),
                length(c(prior_x, posterior$x))),
    
    x = c(prior_x, posterior$x),
    
    density = c(
      prior_y / max(prior_y),
      posterior$y / max(posterior$y)
    ),
    
    type = c(
      rep("Prior", length(prior_x)),
      rep("Posterior", length(posterior$x))
    )
  )
  
}) %>% bind_rows()

# ---------------------------------------------------
# Factor ordering
# ---------------------------------------------------

plot_data$param <- factor(
  plot_data$param,
  levels = unname(param_labels)
)

# ---------------------------------------------------
# Plot
# ---------------------------------------------------

ggplot(plot_data,
       aes(x = x,
           y = density,
           colour = type)) +
  
  # prior centre line
  geom_vline(
    data = prior_centres,
    aes(xintercept = xintercept),
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +
  
  # density curves
  geom_line(linewidth = 1.1) +
  
  # facets
  facet_wrap(
    ~param,
    scales = "free",
    ncol = 2,
    labeller = label_parsed
  ) +
  
  labs(
    x = "Parameter value",
    y = "Normalised density",
    colour = NULL
    #title = "Prior vs Posterior Distributions"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    legend.position = "top",
    strip.text = element_text(size = 13),
    plot.title = element_text(face = "bold")
  )






























################################################################################

# Model 4

# Piecewise constant transmission

################################################################################
 
library(tidyverse)
library(ggplot2)
library(deSolve)
library(tidyr)
library(lubridate)
library(rstan)

 
file_path <- file.path("dta", "MalawiBetaCases.csv")

BetaVariant <- read.csv(file_path) %>%
  mutate(date = lubridate::dmy(date))  

BetaVariant <- BetaVariant%>%
  filter(date >= as.Date("2021-01-02") &
           date <= as.Date("2021-02-14"))


cases <- BetaVariant$new_cases   
date <- BetaVariant$date 


#==============================================
# GAM function to identify turning point(s)
#==============================================
library(mgcv)
library(dplyr)
library(ggplot2)
library(gratia)

 
df <- data.frame(day = date, cases = cases)

df <- df %>%
  arrange(day) %>%
  mutate(day = as.numeric(day - min(day)) + 1)

 
outliers <- c(11, 14)
df$cases[outliers] <- NA


gam_fit <- gam(cases ~ s(day, k = 20), data = df, family = poisson())

# Derivatives of the smooth
deriv <- derivatives(gam_fit, select = "s(day)", order = 1)

# Identify sign changes in derivative
turning_points_idx <- which(diff(sign(deriv$.derivative)) != 0)

# All change points (unfiltered)
t_changes_all <- deriv$day[turning_points_idx]

# Filtered change points (>7 days apart)
t_changes_filtered <- t_changes_all[diff(c(0, t_changes_all)) > 7]

# Add fitted values (must use newdata to get 44 values)
df$smoothed <- predict(gam_fit, newdata = df, type = "response")

# Plot
ggplot(df, aes(x = day)) +
  geom_line(aes(y = cases), color = "steelblue", alpha = 0.4) +
  geom_line(aes(y = smoothed), color = "darkred", linewidth = 1) +
  
  # All change points (grey dashed)
  geom_vline(xintercept = t_changes_all,
             linetype = "dashed", color = "grey50") +
  
  # Filtered change point (purple, thicker)
  geom_vline(xintercept = t_changes_filtered,
             linetype = "dashed", color = "purple", linewidth = 1.2) +
  
  labs(
    title = "Smoothed cases with change points (GAM Derivatives)",
    y = "Cases",
    x = "Day"
  ) +
  theme_minimal()



#==============================================
# Model
#==============================================

library(mgcv)
library(dplyr)
library(ggplot2)
library(gratia)   


df <- data.frame(day = date, cases = cases)

df <- df %>%
  arrange(day) %>%
  mutate(day = as.numeric(day - min(day)) + 1)

 
N <- 17563749   

outliers <- c(11, 14)

n_days <- 44

obs_idx <- setdiff(1:n_days, outliers)
n_obs <- length(obs_idx)   # 42

ts <- 1:(length(cases))

t_changes <- as.array(t_changes_filtered)
n_beta <- length(t_changes) + 1

data_seir4 <- list(
  n_days = n_days,
  n_obs = n_obs,
  obs_idx = obs_idx,
  t0 = 0,
  ts = ts,
  cases = cases,      
  N = N,
  i0 = 5,
  n_beta = n_beta,
  t_changes = t_changes
)

#==============================================
# Compile and fit
#==============================================
model4 <- stan_model("model", "Model_4.stan")

seir_fit_4 <- sampling(model4,
                      data_seir4,
                      iter=4000,
                      warmup = 2000,
                      chains = 4,
                      cores = 4,
                      seed = 42,
                      control = list(adapt_delta = 0.9999, max_treedepth = 20)
                       
)


# Save the results
saveRDS(seir_fit_4, file = "FinalModel_4.rds")

# Load the results
#seir_fit_4 <- readRDS("------------------/"FinalModel_4.RDS")



#==============================================
# Traceplot
#==============================================
stan_trace(seir_fit_4, pars = c("beta_vec[1]", "beta_vec[2]", "a_inv", "gamma_inv", 
                                "e0", "p_reported", "phi_inv"))

#==============================================
#Summary estimates
#==============================================
print(seir_fit_4, pars = c("beta_vec[1]", "beta_vec[2]", "a_inv", "gamma_inv",
                           "e0", "p_reported", "phi_inv", "R0[1]", "R0[2]"))

 
#==============================================
# Density plot
#==============================================
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

draws <- as_draws_df(seir_fit_4)

draws_long <- draws %>%
  select(.chain, `beta_vec[1]`, `beta_vec[2]`, a_inv, gamma_inv, e0, p_reported, 
         phi_inv) %>%
  pivot_longer(
    cols = c(`beta_vec[1]`, `beta_vec[2]`, a_inv, gamma_inv, e0, p_reported, 
             phi_inv),
    names_to = "parameter",
    values_to = "value"
  )

ggplot(draws_long,
       aes(x = value, colour = factor(.chain))) +
  geom_density(linewidth = 1.5) +
  facet_wrap(
    ~ parameter,
    scales = "free",
    labeller = as_labeller(
      c(
        `beta_vec[1]` = "beta[1]",
        `beta_vec[2]` = "beta[2]",
        a_inv = "T[L]",
        gamma_inv = "T[I]",
        e0 = "E[0]",
        p_reported = "pi",
        phi_inv = "phi^{-1}"
        
        
      ),
      label_parsed
    )
  ) +
  labs(colour = "Chain", fill = "Chain") +
  theme_bw()



#==============================================
# Pairwise plot
#==============================================
pairs(seir_fit_4,
      pars = c("beta_vec[1]", "beta_vec[2]", "a_inv",
               "gamma_inv", "phi_inv", "e0" , "p_reported"),
      labels = c(
        expression(beta[1]),
        expression(beta[2]),
        expression(T[L]),
        expression(T[I]),
        expression(phi^{-1}),
        expression(E[0]),
        expression(pi)
      )
)$summary



#==============================================
# PPC
#==============================================
outliers <- c(11, 14)

cases_plot <- cases
cases_plot[outliers] <- NA

pred_summary <- as.data.frame(summary(seir_fit_4, pars = "pred_cases")$summary)

Pred4 <- cbind(
  pred_summary,
  t = ts,            
  cases = cases_plot
)

colnames(Pred4) <- make.names(colnames(Pred4))

ggplot(Pred4, aes(x = t)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "darkred", alpha = 0.35) +
  geom_line(aes(y = X50.), color = "red") +
  geom_point(aes(y = cases), color = "black") +
  labs(
    x = "Day",
    y = "Incidence",
    #title = "Posterior Predictive Check"
  ) +
  theme_minimal()



#==============================================
## Posterior versus posterior plot
#==============================================

library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------
# Extract posterior samples
# ---------------------------------------------------------

post <- as.data.frame(
  rstan::extract(
    seir_fit_4,
    pars = c(
      "beta_vec",
      "a_inv",
      "gamma_inv",
      "e0",
      "p_reported",
      "phi_inv"
    )
  )
)


# =========================================================
# PRIOR VS POSTERIOR:
# beta_vec[1] and beta_vec[2] on SAME AXIS
# =========================================================


beta_prior_mean <- 0.555

# ---------------------------------------------------
# Build combined beta data
# ---------------------------------------------------

beta_data <- do.call(
  rbind,
  lapply(1:2, function(i) {
    
    param_name <- paste0("beta_vec.", i)
    
    posterior_vals <- as.numeric(post[[param_name]])
    
    posterior_vals <- posterior_vals[is.finite(posterior_vals)]
    
    posterior <- density(posterior_vals)
    
    prior_x <- seq(0.01, 1.5, length.out = 300)
    
    prior_y <- dlnorm(
      prior_x,
      meanlog = -0.6094,
      sdlog = 0.2
    )
    
    data.frame(
      beta = paste0("β", i),
      
      x = c(prior_x, posterior$x),
      
      density = c(
        prior_y / max(prior_y),
        posterior$y / max(posterior$y)
      ),
      
      type = c(
        rep("Prior", length(prior_x)),
        rep("Posterior", length(posterior$x))
      )
    )
  })
)

#-----------------------------------------------
#Plot beta parameters together
#-----------------------------------------------
p_beta <- ggplot(
  beta_data,
  aes(
    x = x,
    y = density,
    colour = beta,
    linetype = type
  )
) +

  geom_vline(
    xintercept = beta_prior_mean,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +

  geom_line(linewidth = 1.1) +

  scale_linetype_manual(
    values = c(
      "Prior" = "dashed",
      "Posterior" = "solid"
    )
  ) +

  labs(
    #title = "Prior vs Posterior for Transmission Parameters",
    x = "Parameter value",
    y = "Normalised density",
    colour = NULL,
    linetype = NULL
  ) +

  theme_minimal(base_size = 14) +

  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

print(p_beta)


# =========================================================
# Remaining parameters
# =========================================================

# ---------------------------------------------------------
# Prior density functions
# ---------------------------------------------------------

priors <- list(
  a_inv       = function(x) dlnorm(x, meanlog = 1.6, sdlog = 0.13),
  gamma_inv   = function(x) dlnorm(x, meanlog = 1.6, sdlog = 0.13),
  #e0          = function(x) dlnorm(x, meanlog = log(10), sdlog = 0.8),
  e0          = function(x) dlnorm(10^x, meanlog = log(10), sdlog = 0.8),
  p_reported  = function(x) dbeta(x, shape1 = 0.28, shape2 = 100),
  phi_inv     = function(x) dexp(x, rate = 10)
)


# ---------------------------------------------------------
# Facet labels
# ---------------------------------------------------------

param_labels <- c(
  a_inv      = "T[L]",
  gamma_inv  = "T[I]",
  e0         = "E[0]",
  p_reported = "pi",
  phi_inv    = "phi^{-1}"
  
)

# ---------------------------------------------------------
# Prior centres
# ---------------------------------------------------------

prior_centres <- data.frame(
  
  param = factor(
    unname(param_labels),
    levels = unname(param_labels)
  ),
  
  xintercept = c(
    
    exp(1.6 + 0.13^2/2),               # a_inv mean
    exp(1.6 + 0.13^2/2),               # gamma_inv mean
    #exp(log(10) + 0.8^2/2),            # e0 mean
    log10(exp(log(10) + 0.8^2/2)),    # e0 mean
    0.28/(0.28 + 100),                 # p_reported mean
    1/10                              # phi_inv mean
   )
)



#-----------------------------------------------------------
# Build plotting data
# ---------------------------------------------------------

params_to_plot <- names(param_labels)

plot_data <- lapply(params_to_plot, function(p) {
  
  #posterior <- density(post[[p]])
  posterior <- if (p == "e0") {
    density(log10(post[[p]]))
  } else {
    density(post[[p]])
  }
  
  prior_x <- switch(
    p,
    
    a_inv       = seq(0, 8, length.out = 300),
    gamma_inv   = seq(0, 8, length.out = 300),
   # e0          = seq(0, 10000, length.out = 300),
    e0          = seq(log10(1), log10(30000), length.out = 300),
    p_reported  = seq(0, 0.15, length.out = 300),
    phi_inv     = seq(0, 0.2, length.out = 300),
    
    seq(
      min(posterior$x),
      max(posterior$x),
      length.out = 300
    )
  )
  
  #prior_y <- priors[[p]](prior_x)
  prior_y <- if (p == "e0") {
    dlnorm(10^prior_x, meanlog = log(10), sdlog = 0.8)
  } else {
    priors[[p]](prior_x)
  }
  
  
  data.frame(
    
    param = rep(
      unname(param_labels[p]),
      length(c(prior_x, posterior$x))
    ),
    
    #x = c(prior_x, posterior$x),
    x = c(
      if (p == "e0") prior_x else prior_x,
      if (p == "e0") posterior$x else posterior$x
    ),
    
    density = c(
      prior_y / max(prior_y),
      posterior$y / max(posterior$y)
    ),
    
    type = c(
      rep("Prior", length(prior_x)),
      rep("Posterior", length(posterior$x))
    )
  )
  
}) %>%
  bind_rows()



# ---------------------------------------------------------
# Factor ordering
# ---------------------------------------------------------

plot_data$param <- factor(
  plot_data$param,
  levels = unname(param_labels)
)

# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------

ggplot(
  plot_data,
  aes(
    x = x,
    y = density,
    colour = type
  )
) +
  
  geom_vline(
    data = prior_centres,
    aes(xintercept = xintercept),
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +
  
  geom_line(linewidth = 1.1) +
  
  facet_wrap(
    ~param,
    scales = "free",
    ncol = 2,
    labeller = label_parsed
  ) +
  
  labs(
    x = "Parameter value",
    y = "Normalised density",
    colour = NULL,
    #title = "Prior vs Posterior Distributions"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    legend.position = "top",
    strip.text = element_text(size = 13),
    plot.title = element_text(face = "bold")
  )






























################################################################################

# Model 5 

# Linearised SEIR - Piecewise exponential function

################################################################################
library(tidyverse)
library(ggplot2)
library(deSolve)
library(tidyr)
library(lubridate)

library(rstan)

set.seed(42)

file_path <- file.path("data", "MalawiBetaCases.csv") 

BetaVariant <- read.csv(file_path) %>%
  mutate(date = lubridate::dmy(date))   

BetaVariant <- BetaVariant%>%
  filter(date >= as.Date("2021-01-02") &
           date <= as.Date("2021-02-14"))


cases <- c(BetaVariant$new_cases)   

n_days <- 44

outliers <- c(11, 14)

obs_idx <- setdiff(1:n_days, outliers)

n_obs <- length(obs_idx)

time <- 0:(n_days - 1)

stan_data_NonODE <- list(
  n_days = n_days,
  n_obs = n_obs,
  obs_idx = obs_idx,
  time = time,
  cases = cases,
  T_change = 20.5   
)


#==============================================
#Compile and fit
#==============================================
 
model5 <- stan_model("model", "Model_5.stan")
 

seir_fit_5 <- sampling(
  model5,
  data = stan_data_NonODE,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123
)

# # Save the results
saveRDS(seir_fit_5, file = "FinalModel_5.rds")

# Load the results
#seir_fit_5 <- readRDS("------------------/FinalModel_5.RDS")


#==============================================
# Traceplot
#==============================================
stan_trace(seir_fit_5, pars = c("C0", "r", "r_minus", "phi_inv")) 

#==============================================
# Summary estimates
#==============================================
print(seir_fit_5, pars = c("C0", "r", "r_minus", "phi_inv"))

#==============================================
# Density plot
#==============================================
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

draws <- as_draws_df(seir_fit_5)

draws_long <- draws %>%
  select(.chain, C0, r, r_minus, phi_inv) %>%
  pivot_longer(
    cols = c(C0, r, r_minus, phi_inv),
    names_to = "parameter",
    values_to = "value"
  )

ggplot(draws_long,
       aes(x = value, colour = factor(.chain))) +
  geom_density(linewidth = 1.5) +
  facet_wrap(
    ~ parameter,
    scales = "free",
    labeller = as_labeller(
      c(
        C0 = "C[0]",
        r = "r",
        r_minus = "r[`-`]",
        phi_inv = "phi^{-1}"
      ),
      label_parsed
    )
  ) +
  labs(colour = "Chain", fill = "Chain") +
  theme_bw()




#==============================================
# Pairwise plot
#==============================================
pairs(seir_fit_5,
      pars = c("C0", "r", "r_minus", "phi_inv"),
      labels = c(
        expression(C[0]),
        expression(r),
        expression(r[`-`]),
        expression(phi^{-1})
        
        
      )
)$summary



#==============================================
# PPC
#==============================================

outliers <- c(11, 14)

cases_plot <- cases
cases_plot[outliers] <- NA    


pred_summary <- as.data.frame(summary(seir_fit_5, pars = "pred_cases")$summary)

Pred5 <- cbind(
  pred_summary,
  t = time,           
  cases = cases_plot 
)

colnames(Pred5) <- make.names(colnames(Pred5))


ggplot(Pred5, aes(x = t)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "darkred", alpha = 0.35) +
  geom_line(aes(y = X50., color = "Median prediction")) +
  geom_point(aes(y = cases, color = "Observed cases")) +
  scale_color_manual(values = c("Median prediction" = "red", "Observed cases" = "black")) +
  labs(
    #title = "Posterior Predictive Check",
    x = "Day",
    y = "Incidence",
    color = "Legend"
  ) +
  theme_minimal()


#==============================================
# Prior vs posterior plot
#==============================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

#--------------------------------------------
# Extract posterior samples
#--------------------------------------------
post <- as.data.frame(
  rstan::extract(
    seir_fit_5,
    pars = c(
      "C0",
      "r",
      "r_minus",
      "phi_inv"
      
    )
  )
)


priors <- list(
  C0        = function(x) dlnorm(x, meanlog = log(5), sdlog = 0.8),
  r         = function(x) dnorm(x, mean = 0.2, sd = 0.2),
  r_minus   = function(x) dnorm(x, mean = 0.1, sd = 0.1),
  phi_inv   = function(x) dexp(x, rate = 10)
)


# ---------------------------------------------------
# Facet labels
# ---------------------------------------------------

param_labels <- c(
  C0        = "C[0]",
  r      = "r",
  r_minus = "r[`-`]",
  phi_inv    = "phi^{-1}"
)

 

prior_centres <- data.frame(
  param = factor(
    unname(param_labels),
    levels = unname(param_labels)
  ),
  xintercept <- c(
    exp(log(5) + 0.8^2/2),   # C0 mean ≈ 6.885
    0.2,                     # r mean
    0.1,                     # r_minus mean
    1/10                     # phi_inv mean = 0.1
  )
  
)
# ---------------------------------------------------
# Build plotting data
# ---------------------------------------------------

plot_data <- lapply(names(post), function(p) {
  
  posterior <- density(post[[p]])
  
  prior_x <- switch(
    p,
    
    C0        = seq(0, 150, length.out = 300),
    r         = seq(0, 1, length.out = 300),
    r_minus = seq(0, 0.15, length.out = 300),
    phi_inv    = seq(0, 0.5, length.out = 300),
    
    seq(
      min(posterior$x),
      max(posterior$x),
      length.out = 300
    )
  )
  
  prior_y <- priors[[p]](prior_x)
  
  data.frame(
    
    param = rep(unname(param_labels[p]),
                length(c(prior_x, posterior$x))),
    
    x = c(prior_x, posterior$x),
    
    density = c(
      prior_y / max(prior_y),
      posterior$y / max(posterior$y)
    ),
    
    type = c(
      rep("Prior", length(prior_x)),
      rep("Posterior", length(posterior$x))
    )
  )
  
}) %>% bind_rows()

# ---------------------------------------------------
# Factor ordering
# ---------------------------------------------------

plot_data$param <- factor(
  plot_data$param,
  levels = unname(param_labels)
)

# ---------------------------------------------------
# Plot
# ---------------------------------------------------

ggplot(plot_data,
       aes(x = x,
           y = density,
           colour = type)) +
  
  # prior centre line
  geom_vline(
    data = prior_centres,
    aes(xintercept = xintercept),
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +
  
  # density curves
  geom_line(linewidth = 1.1) +
  
  # facets
  facet_wrap(
    ~param,
    scales = "free",
    ncol = 2,
    labeller = label_parsed
  ) +
  
  labs(
    x = "Parameter value",
    y = "Normalised density",
    colour = NULL
    #title = "Prior vs Posterior Distributions"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    legend.position = "top",
    strip.text = element_text(size = 13),
    plot.title = element_text(face = "bold")
  )



 


























################################################################################

# Model 1 with longer epidemic window

# Full SEIR

################################################################################

library(tidyverse)
library(ggplot2)
library(deSolve)
library(tidyr)
library(lubridate)

library(rstan)

set.seed(42)
# Load Data
file_path <- file.path("data", "MalawiBetaCases.csv")

BetaVariant <- read.csv(file_path) %>%
  mutate(date = lubridate::dmy(date))   



cases <- BetaVariant$new_cases    
N <- 17563749   

all_days <- length(cases)        

 
obs_days <- setdiff(1:all_days, c(41, 44))

cases_obs <- cases[obs_days]
n_days <- length(cases_obs)       

burn_in_days <- 32
ts <- seq_len(all_days + burn_in_days)

stan_data_6 <- list(
  n_days = n_days,
  n_total_days = all_days,
  burn_in_days = burn_in_days,
  obs_days = obs_days,
  t0 = 0,
  ts = ts,
  cases = cases_obs,
  N = N
)



#==============================================
# # Compile and fit
#==============================================
model6 <- stan_model("model", "Model_1.stan")

seir_fit_6 <- sampling(model6,
                  data = stan_data_6,
                  chains = 4,
                  cores = 4,
                  iter = 4000,
                  warmup = 2000,
                  control = list(adapt_delta = 0.9999, max_treedepth = 20),
                  seed = 123)

# Save the results
saveRDS(seir_fit_6, file = "FinalModel_6.rds")

# Load the results
#seir_fit_6 <- readRDS("------------------/FinalModel_6.RDS")


#==============================================
# Calculate \beta  
#==============================================
post <- extract(seir_fit_6)
beta <- post$R_0 * post$gamma    
mean_beta <- mean(beta)
crI_beta <- quantile(beta, c(0.025, 0.975))
mean_beta
crI_beta


#==============================================
# Calculate \pi to 4 decimal places
#==============================================
post <- rstan::extract(seir_fit_6)
p <- post$p_reported
mean_p <- mean(p)
crI_p <- quantile(p, c(0.025, 0.975))
sprintf("%.4f (%.4f, %.4f)",
        mean_p, crI_p[1], crI_p[2])


#==============================================
# Traceplot
#==============================================
stan_trace(seir_fit_6, pars = c("R_0", "alpha_inv", "gamma_inv", "E0", "p_reported", "phi_inv")) 


#==============================================
# Posterior summary
#==============================================
print(seir_fit_6, pars = c("R_0", "alpha_inv", "gamma_inv", "E0", "p_reported", "phi_inv"))


#==============================================
# Density plot
#==============================================
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

draws <- as_draws_df(seir_fit_6)

draws_long <- draws %>%
  select(.chain, R_0, alpha_inv, gamma_inv, E0, p_reported, phi_inv) %>%
  pivot_longer(
    cols = c(R_0, alpha_inv, gamma_inv, E0, p_reported, phi_inv),
    names_to = "parameter",
    values_to = "value"
  )

ggplot(draws_long,
       aes(x = value, colour = factor(.chain))) +
  geom_density(linewidth = 1.5) +
  facet_wrap(
    ~ parameter,
    scales = "free",
    labeller = as_labeller(
      c(
        R_0 = "R[0]",
        alpha_inv = "T[L]",
        gamma_inv = "T[I]",
        #I0 = "I[0]",
        E0 = "E[0]",
        p_reported = "pi",
        phi_inv = "phi^{-1}"
      ),
      label_parsed
    )
  ) +
  labs(colour = "Chain", fill = "Chain") +
  theme_bw()


#==============================================
# Pairwiseplot
#==============================================
pairs(seir_fit_6,
      pars = c("R_0", "alpha_inv", "gamma_inv", "phi_inv", "p_reported", "I0", "E0"),
      labels = c(
        expression(R[0]),
        expression(T[L]),
        expression(T[I]),
        #expression(I[0]),
        expression(E[0]),
        expression(pi),
        expression(phi^{-1})
      )
)$summary


#==============================================
## Posterior predictive check
#==============================================

outliers <- c(41, 44)

cases <- BetaVariant$new_cases   # length 44

# days used in likelihood (2..44, excluding outliers)
obs_days <- setdiff(1:length(cases), outliers)

cases_obs <- cases[obs_days]

# posterior predictive summary
pred_summary <- as.data.frame(
  summary(seir_fit_6, pars = "pred_cases")$summary
)

# time index for intervals: 1..(n_total_days - 1)
t_pred <- 1:(length(cases))

# match intervals to observed days: obs_days - 1
interval_idx <- obs_days
pred_summary_obs <- pred_summary[interval_idx, ]

Pred6 <- cbind(
  pred_summary_obs,
  day = obs_days,      # calendar day (2..44, no outliers)
  cases = cases_obs
)

colnames(Pred6) <- make.names(colnames(Pred6))

ggplot(Pred6, aes(x = day)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "darkred", alpha = 0.35) +
  geom_line(aes(y = X50., color = "Median prediction")) +
  geom_point(aes(y = cases, color = "Observed cases")) +
  scale_color_manual(values = c("Median prediction" = "red", "Observed cases" = "black")) +
  labs(
    #title = "Posterior Predictive Check",
    x = "Day",
    y = "Incidence (interval)",
    color = "Legend"
  ) +
  theme_minimal()


#==============================================
# Prior versus posterior distribution
#==============================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

# Extract posterior samples
post <- as.data.frame(
  rstan::extract(
    seir_fit_6,
    pars = c(
      "R_0",
      "alpha_inv",
      "gamma_inv",
      "phi_inv",
      "p_reported",
      "E0"
    )
  )
)


# Prior density functions
priors <- list(
  R_0        = function(x) dlnorm(x, 1.009619, 0.2),
  alpha_inv      = function(x) dlnorm(x, log(4.5), 0.5),
  gamma_inv      = function(x) dlnorm(x, log(5), 0.5),
  phi_inv    = function(x) dexp(x, 10),
  p_reported = function(x) dbeta(x, 0.28, 100),
  #E0         = function(x) dlnorm(x, log(10), 0.5)
  E0         = function(x) dlnorm(10^x, log(10), 0.5)
)

# Facet labels
param_labels <- c(
  R_0        = "R[0]",
  alpha_inv      = "T[L]",
  gamma_inv      = "T[I]",
  phi_inv    = "phi^{-1}",
  p_reported = "pi",
  E0         = "E[0]"
)

# Prior centres for vertical lines
prior_centres <- data.frame(
  
  param = factor(
    unname(param_labels),
    levels = unname(param_labels)
  ),
  
  xintercept = c(
    exp(1.009619 + 0.2^2/2),   #R_0
    exp(log(5) + 0.13^2/2),   # T_L     
    exp(log(5) + 0.13^2/2),  #T_I    
    1/10,         # phi_inv    
    0.0028,       # p_reported
    log10(exp(log(10) + 0.8^2/2)) # E0       
  )
)

# Build plotting data
plot_data <- lapply(names(post), function(p) {
  #posterior <- density(post[[p]])
  posterior <- if (p == "E0") {
    density(log10(post[[p]]))
  } else {
    density(post[[p]])
  }
  prior_x <- switch(p,
                    p_reported = seq(0, 0.005, length.out = 300),
                    phi_inv    = seq(0, 0.5, length.out = 300),
                    alpha_inv      = seq(0, 12, length.out = 300),
                    gamma_inv      = seq(0, 6, length.out = 300),
                    R_0        = seq(0, 6, length.out = 300),
                    #E0         = seq(0, 300, length.out = 300),
                    E0         = seq(log10(1), log10(4000), length.out = 300),
                    
                    seq(
                      min(posterior$x),
                      max(posterior$x),
                      length.out = 300
                    )
  )
  
  #prior_y <- priors[[p]](prior_x)
  prior_y <- if (p == "E0"){
    dlnorm(10^prior_x, meanlog=log(10), sdlog = 0.8)
  } else {
    priors[[p]](prior_x)
  }
  
  data.frame(
    
    param = rep(unname(param_labels[p]),
                length(c(prior_x, posterior$x))),
    
    #x = c(prior_x, posterior$x),
    x = c(
      if (p == "E0") prior_x else prior_x,
      if (p == "E0") posterior$x else posterior$x
    ),
    
    density = c(
      prior_y / max(prior_y),
      posterior$y / max(posterior$y)
    ),
    
    type = c(
      rep("Prior", length(prior_x)),
      rep("Posterior", length(posterior$x))
    )
  )
  
}) %>% bind_rows()

# Factor ordering
plot_data$param <- factor(plot_data$param,levels = unname(param_labels))

# Plot
ggplot(plot_data, aes(x = x, y = density, colour = type)) +
  
  # prior centre line
  geom_vline(
    data = prior_centres,
    aes(xintercept = xintercept),
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +
  
  # density curves
  geom_line(linewidth = 1.1) +
  
  # facets
  facet_wrap(
    ~param,
    scales = "free",
    ncol = 2,
    labeller = label_parsed
  ) +
  
  labs(
    x = "Parameter value",
    y = "Normalised density",
    colour = NULL,
    #title = "Prior vs Posterior Distributions"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    legend.position = "top",
    strip.text = element_text(size = 13),
    plot.title = element_text(face = "bold")
  )




################################################################################
################################################################################
################################################################################
################################################################################




