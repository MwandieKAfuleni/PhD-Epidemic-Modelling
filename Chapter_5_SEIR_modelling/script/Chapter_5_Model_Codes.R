################################################################################

# Model 1

# Full SEIR

# Free I0

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
model1 <- stan_model("model", "Model_1_FreeI0.stan")
 
seir_fit_1 <- sampling(model1,
                       data = stan_data_1,
                       chains = 4,
                       cores = 4,
                       iter = 4000,
                       warmup = 2000,
                       control = list(adapt_delta = 0.9999, max_treedepth = 20),
                       seed = 123)


# Save the results
saveRDS(seir_fit_1, file = "FinalModel_1_FreeI0.rds")

# Load the results
#seir_fit_1 <- readRDS(------------------/FinalModel_1_FreeI0.rds")



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








# =======================================================================
# Pairs plot with posterior histograms and prior densities on diagonals
# =======================================================================

library(rstan)
library(ggplot2)
library(GGally)
library(rlang)

# ============================================================
# Extract posterior samples
# ============================================================

posterior_pairs <- as.data.frame(
  rstan::extract(
    seir_fit_1,
    pars = c(
      "R_0",
      "alpha_inv",
      "gamma_inv",
      "phi_inv",
      "p_reported",
      "I0",
      "E0"
    )
  )
)

# ============================================================
# Prior density functions
# ============================================================

prior_density <- list(
  
  R_0 = function(x) {
    dlnorm(x, meanlog = 1.0, sdlog = 0.2)
  },
  
  alpha_inv = function(x) {
    dlnorm(x, meanlog = 1.6, sdlog = 0.13)
  },
  
  gamma_inv = function(x) {
    dlnorm(x, meanlog = 1.6, sdlog = 0.13)
  },
  
  phi_inv = function(x) {
    dexp(x, rate = 10)
  },
  
  p_reported = function(x) {
    dbeta(x, shape1 = 0.28, shape2 = 100)
  },
  
  I0 = function(x) {
    dlnorm(x, meanlog = log(5), sdlog = 0.8)
  },
  
  E0 = function(x) {
    dlnorm(x, meanlog = log(10), sdlog = 0.8)
  }
)

# ============================================================
# Custom diagonal panel
# Posterior histogram + prior density
# ============================================================

posterior_prior_diag <- function(data, mapping, ...) {
  
  # Get parameter name
  parameter <- rlang::as_name(mapping$x)
  
  # Posterior samples
  posterior_values <- data[[parameter]]
  
  # Remove NA / non-finite values
  posterior_values <- posterior_values[
    is.finite(posterior_values)
  ]
  
  # ----------------------------------------------------------
  # Set sensible x-axis limits
  # ----------------------------------------------------------
  
  x_upper <- max(posterior_values, na.rm = TRUE)
  
  x_limits <- switch(
    parameter,
    
    R_0 = c(0, max(6, x_upper)),
    
    alpha_inv = c(0, max(8, x_upper)),
    
    gamma_inv = c(0, max(8, x_upper)),
    
    phi_inv = c(0, max(0.5, x_upper)),
    
    p_reported = c(0, max(0.01, x_upper)),
    
    I0 = c(0, max(100, x_upper)),
    
    E0 = c(0, max(200, x_upper)),
    
    c(0, x_upper)
  )
  
  # ----------------------------------------------------------
  # Create grid for prior density
  # ----------------------------------------------------------
  
  prior_grid <- data.frame(
    x = seq(
      x_limits[1],
      x_limits[2],
      length.out = 500
    )
  )
  
  # Evaluate prior density
  prior_grid$density <- prior_density[[parameter]](
    prior_grid$x
  )
  
  # ----------------------------------------------------------
  # Create diagonal plot
  # ----------------------------------------------------------
  
  ggplot(data, aes(x = .data[[parameter]])) +
    
    # Posterior histogram
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 30,
      fill = "#4C78A8",
      colour = "white",
      alpha = 0.65
    ) +
    
    # Prior density
    geom_line(
      data = prior_grid,
      aes(
        x = x,
        y = density
      ),
      inherit.aes = FALSE,
      colour = "#D62728",
      linewidth = 1
    ) +
    
    coord_cartesian(
      xlim = x_limits,
      expand = FALSE
    ) +
    
    theme_minimal(base_size = 11) +
    
    theme(
      axis.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# Create pairs plot
# ============================================================

pair_plot <- GGally::ggpairs(
  
  posterior_pairs,
  
  # ----------------------------------------------------------
  # Upper triangle: posterior correlations
  # ----------------------------------------------------------
  
  upper = list(
    continuous = GGally::wrap(
      "cor",
      size = 4
    )
  ),
  
  # ----------------------------------------------------------
  # Lower triangle: posterior samples
  # ----------------------------------------------------------
  
  lower = list(
    continuous = GGally::wrap(
      "points",
      alpha = 0.25,
      size = 0.4,
      colour = "#4C78A8"
    )
  ),
  
  # ----------------------------------------------------------
  # Diagonal: posterior histogram + prior density
  # ----------------------------------------------------------
  
  diag = list(
    continuous = posterior_prior_diag
  ),
  
  # ----------------------------------------------------------
  # Parameter labels
  # ----------------------------------------------------------
  
  columnLabels = c(
    "R[0]",
    "T[L]",
    "T[I]",
    "phi^{-1}",
    "pi",
    "I[0]",
    "E[0]"
  ),
  
  labeller = "label_parsed",
  
  progress = FALSE
  
) +
  
  theme_bw(base_size = 12) +
  
  theme(
    strip.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

# ============================================================
# Display
# ============================================================

pair_plot




#==============================================
## Posterior predictive check
#==============================================

outliers <- c(11, 14)

cases <- BetaVariant$new_cases   # length 44

# days used in likelihood (2..44, excluding outliers)
obs_days <- setdiff(1:length(cases), outliers)

cases_obs <- cases[obs_days]

# posterior predictive summary
pred_summary <- as.data.frame(
  summary(seir_fit_1, pars = "pred_cases")$summary
)

# time index for intervals: 1..(n_total_days - 1)
t_pred <- 1:(length(cases))

# match intervals to observed days: obs_days - 1
interval_idx <- obs_days
pred_summary_obs <- pred_summary[interval_idx, ]

Pred1 <- cbind(
  pred_summary_obs,
  day = obs_days,      # calendar day (2..44, no outliers)
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






# =======================================================================
# Pairs plot with posterior histograms and prior densities on diagonals
# =======================================================================

library(rstan)
library(ggplot2)
library(GGally)
library(rlang)

# ============================================================
# Extract posterior samples
# ============================================================

posterior_pairs <- as.data.frame(
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

# Retain only finite posterior samples
posterior_pairs <- posterior_pairs[
  apply(posterior_pairs, 1, function(x) all(is.finite(x))),
]

# ============================================================
# Prior density functions
# ============================================================

prior_density <- list(
  
  R_0 = function(x) {
    dlnorm(
      x,
      meanlog = 1.0,
      sdlog = 0.2
    )
  },
  
  I0 = function(x) {
    dlnorm(
      x,
      meanlog = log(5),
      sdlog = 0.8
    )
  },
  
  p_reported = function(x) {
    dbeta(
      x,
      shape1 = 0.28,
      shape2 = 100
    )
  },
  
  phi_inv = function(x) {
    dexp(
      x,
      rate = 10
    )
  }
)

# ============================================================
# Custom diagonal panel
# Posterior histogram + prior density
# ============================================================

posterior_prior_diag <- function(data, mapping, ...) {
  
  # Identify the parameter represented in the panel
  parameter <- rlang::as_label(mapping$x)
  
  # Extract finite posterior samples
  posterior_values <- data[[parameter]]
  posterior_values <- posterior_values[
    is.finite(posterior_values)
  ]
  
  # ----------------------------------------------------------
  # Set parameter-specific x-axis limits
  # ----------------------------------------------------------
  
  posterior_upper <- max(posterior_values, na.rm = TRUE)
  
  x_limits <- switch(
    parameter,
    
    R_0 = c(
      0,
      max(4, posterior_upper)
    ),
    
    I0 = c(
      0,
      max(20, posterior_upper)
    ),
    
    p_reported = c(
      0,
      max(0.005, posterior_upper)
    ),
    
    phi_inv = c(
      0,
      max(0.5, posterior_upper)
    ),
    
    range(posterior_values, na.rm = TRUE)
  )
  
  # ----------------------------------------------------------
  # Construct the prior-density grid
  # ----------------------------------------------------------
  
  # Avoid evaluating beta and lognormal densities exactly at zero
  grid_lower <- if (
    parameter %in% c("R_0", "I0", "p_reported")
  ) {
    max(1e-8, x_limits[1])
  } else {
    x_limits[1]
  }
  
  prior_grid <- data.frame(
    x = seq(
      from = grid_lower,
      to = x_limits[2],
      length.out = 500
    )
  )
  
  # prior_grid$density <- prior_densityprior_grid$x
  # 
  # # Remove any non-finite density values
  # prior_grid <- prior_grid[
  #   is.finite(prior_grid$density),
  # ]
  
  
  
  
  # ----------------------------------------------------------
  # Construct the prior-density grid
  # ----------------------------------------------------------
  
  # Avoid evaluating beta and lognormal densities exactly at zero
  grid_lower <- if (
    parameter %in% c("R_0", "I0", "p_reported")
  ) {
    max(1e-8, x_limits[1])
  } else {
    x_limits[1]
  }
  
  prior_grid <- data.frame(
    x = seq(
      from = grid_lower,
      to = x_limits[2],
      length.out = 500
    )
  )
  
  # Evaluate the appropriate prior density
  prior_grid$density <- prior_density[[parameter]](
    prior_grid$x
  )
  
  # Remove any non-finite density values
  prior_grid <- prior_grid[
    is.finite(prior_grid$density),
  ]
  
  # ----------------------------------------------------------
  # Construct the diagonal plot
  # ----------------------------------------------------------
  
  ggplot(
    data,
    aes(x = .data[[parameter]])
  ) +
    
    # Posterior histogram
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 30,
      fill = "#4C78A8",
      colour = "white",
      alpha = 0.65
    ) +
    
    # Prior density
    geom_line(
      data = prior_grid,
      aes(
        x = x,
        y = density
      ),
      inherit.aes = FALSE,
      colour = "#D62728",
      linewidth = 1
    ) +
    
    coord_cartesian(
      xlim = x_limits,
      expand = FALSE
    ) +
    
    theme_minimal(base_size = 11) +
    
    theme(
      axis.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# Create pairs plot
# ============================================================

pair_plot <- GGally::ggpairs(
  
  posterior_pairs,
  
  # Upper triangle: posterior correlations
  upper = list(
    continuous = GGally::wrap(
      "cor",
      size = 4
    )
  ),
  
  # Lower triangle: posterior samples
  lower = list(
    continuous = GGally::wrap(
      "points",
      alpha = 0.25,
      size = 0.4,
      colour = "#4C78A8"
    )
  ),
  
  # Diagonal: posterior histogram and prior density
  diag = list(
    continuous = posterior_prior_diag
  ),
  
  # Mathematical parameter labels
  columnLabels = c(
    "R[0]",
    "I[0]",
    "pi",
    "phi^{-1}"
  ),
  
  labeller = "label_parsed",
  
  progress = FALSE
  
) +
  
  theme_bw(base_size = 12) +
  
  theme(
    strip.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

# ============================================================
# Display plot
# ============================================================

pair_plot


#









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














#
































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







# =======================================================================
# Pairs plot with posterior histograms and prior densities on diagonals
# Model: seir_fit_3
# =======================================================================

library(rstan)
library(ggplot2)
library(GGally)
library(rlang)

# ============================================================
# Extract posterior samples
# ============================================================

posterior_pairs <- as.data.frame(
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

# ============================================================
# Transform e0 to the log10 scale
# ============================================================

# The posterior histogram and corresponding prior density
# for e0 are both plotted on the log10 scale.

posterior_pairs$e0 <- log10(posterior_pairs$e0)

# ============================================================
# Retain only finite posterior samples
# ============================================================

posterior_pairs <- posterior_pairs[
  apply(
    posterior_pairs,
    1,
    function(x) all(is.finite(x))
  ),
]

# ============================================================
# Prior density functions
# ============================================================

prior_density <- list(
  
  beta = function(x) {
    dlnorm(
      x,
      meanlog = -0.6094,
      sdlog = 0.2
    )
  },
  
  a_inv = function(x) {
    dlnorm(
      x,
      meanlog = 1.6,
      sdlog = 0.13
    )
  },
  
  gamma_inv = function(x) {
    dlnorm(
      x,
      meanlog = 1.6,
      sdlog = 0.13
    )
  },
  
  # Prior density of log10(e0), including the Jacobian
  e0 = function(x) {
    
    e0_original <- 10^x
    
    dlnorm(
      e0_original,
      meanlog = log(10),
      sdlog = 0.8
    ) *
      e0_original *
      log(10)
  },
  
  p_reported = function(x) {
    dbeta(
      x,
      shape1 = 0.28,
      shape2 = 100
    )
  },
  
  phi_inv = function(x) {
    dexp(
      x,
      rate = 10
    )
  },
  
  eta = function(x) {
    dbeta(
      x,
      shape1 = 2.5,
      shape2 = 4
    )
  },
  
  nu = function(x) {
    dexp(
      x,
      rate = 1 / 5
    )
  },
  
  xi_raw = function(x) {
    dbeta(
      x,
      shape1 = 2,
      shape2 = 5
    )
  },
  
  omega = function(x) {
    dbeta(
      x,
      shape1 = 2,
      shape2 = 5
    )
  }
)

# ============================================================
# Parameter-specific plotting ranges
# ============================================================

prior_ranges <- list(
  
  beta       = c(0, 1),
  a_inv      = c(0, 8),
  gamma_inv  = c(0, 8),
  
  # log10 scale: e0 from 1 to 10,000
  e0         = c(log10(1), log10(10000)),
  
  p_reported = c(0, 0.15),
  phi_inv    = c(0, 1),
  eta        = c(0, 0.75),
  nu         = c(0, 7),
  xi_raw     = c(0, 0.5),
  omega      = c(0, 0.5)
)

# ============================================================
# Custom diagonal panel
# Posterior histogram + prior density
# ============================================================

posterior_prior_diag <- function(data, mapping, ...) {
  
  # Identify the parameter represented in the panel
  parameter <- rlang::as_label(mapping$x)
  
  # Extract finite posterior samples
  posterior_values <- data[[parameter]]
  
  posterior_values <- posterior_values[
    is.finite(posterior_values)
  ]
  
  # ----------------------------------------------------------
  # Set x-axis limits
  # ----------------------------------------------------------
  
  posterior_range <- range(
    posterior_values,
    na.rm = TRUE
  )
  
  specified_range <- prior_ranges[[parameter]]
  
  x_limits <- c(
    min(
      specified_range[1],
      posterior_range[1]
    ),
    max(
      specified_range[2],
      posterior_range[2]
    )
  )
  
  # Keep positive-support parameters above zero
  if (
    parameter %in%
    c(
      "beta",
      "a_inv",
      "gamma_inv",
      "p_reported",
      "phi_inv",
      "eta",
      "nu",
      "xi_raw",
      "omega"
    )
  ) {
    x_limits[1] <- max(
      0,
      x_limits[1]
    )
  }
  
  # ----------------------------------------------------------
  # Create the prior-density grid
  # ----------------------------------------------------------
  
  grid_lower <- x_limits[1]
  grid_upper <- x_limits[2]
  
  # Avoid evaluating positive-support distributions at zero
  if (
    parameter %in%
    c(
      "beta",
      "a_inv",
      "gamma_inv",
      "p_reported",
      "phi_inv",
      "eta",
      "nu",
      "xi_raw",
      "omega"
    )
  ) {
    grid_lower <- max(
      grid_lower,
      1e-8
    )
  }
  
  # Keep beta-distribution parameters strictly inside (0, 1)
  if (
    parameter %in%
    c(
      "p_reported",
      "eta",
      "xi_raw",
      "omega"
    )
  ) {
    grid_upper <- min(
      grid_upper,
      1 - 1e-8
    )
  }
  
  prior_grid <- data.frame(
    x = seq(
      from = grid_lower,
      to = grid_upper,
      length.out = 500
    )
  )
  
  # Evaluate the appropriate prior density
  prior_grid$density <- prior_density[[parameter]](
    prior_grid$x
  )
  
  # Remove non-finite prior-density values
  prior_grid <- prior_grid[
    is.finite(prior_grid$density),
  ]
  
  # ----------------------------------------------------------
  # Construct the diagonal panel
  # ----------------------------------------------------------
  
  ggplot(
    data,
    aes(x = .data[[parameter]])
  ) +
    
    # Posterior histogram
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 30,
      fill = "#4C78A8",
      colour = "white",
      alpha = 0.65
    ) +
    
    # Prior density
    geom_line(
      data = prior_grid,
      aes(
        x = x,
        y = density
      ),
      inherit.aes = FALSE,
      colour = "#D62728",
      linewidth = 1
    ) +
    
    coord_cartesian(
      xlim = x_limits,
      expand = FALSE
    ) +
    
    theme_minimal(
      base_size = 10
    ) +
    
    theme(
      axis.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# Create pairs plot
# ============================================================

pair_plot <- GGally::ggpairs(
  
  posterior_pairs,
  
  # Upper triangle: posterior correlations
  upper = list(
    continuous = GGally::wrap(
      "cor",
      size = 3
    )
  ),
  
  # Lower triangle: posterior samples
  lower = list(
    continuous = GGally::wrap(
      "points",
      alpha = 0.20,
      size = 0.30,
      colour = "#4C78A8"
    )
  ),
  
  # Diagonal: posterior histogram + prior density
  diag = list(
    continuous = posterior_prior_diag
  ),
  
  # Mathematical parameter labels
  columnLabels = c(
    "beta",
    "T[L]",
    "T[I]",
    "log[10](E[0])",
    "pi",
    "phi^{-1}",
    "eta",
    "nu",
    "xi",
    "omega"
  ),
  
  labeller = "label_parsed",
  
  progress = FALSE
  
) +
  
  theme_bw(
    base_size = 11
  ) +
  
  theme(
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# ============================================================
# Display plot
# ============================================================

pair_plot










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













# 





































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

 






# =======================================================================
# Pairs plot with posterior histograms and prior densities on diagonals
# Model: seir_fit_4
# =======================================================================

library(rstan)
library(ggplot2)
library(GGally)
library(rlang)

# ============================================================
# Extract posterior samples
# ============================================================

posterior_pairs <- as.data.frame(
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

# ============================================================
# Rename beta parameters
# ============================================================

names(posterior_pairs)[
  names(posterior_pairs) == "beta_vec.1"
] <- "beta1"

names(posterior_pairs)[
  names(posterior_pairs) == "beta_vec.2"
] <- "beta2"

# ============================================================
# Transform e0 to the log10 scale
# ============================================================

# The posterior histogram and corresponding prior density
# for e0 are both plotted on the log10 scale.

posterior_pairs$e0 <- log10(
  posterior_pairs$e0
)

# ============================================================
# Retain only finite posterior samples
# ============================================================

posterior_pairs <- posterior_pairs[
  apply(
    posterior_pairs,
    1,
    function(x) all(is.finite(x))
  ),
]

# ============================================================
# Prior density functions
# ============================================================

prior_density <- list(
  
  beta1 = function(x) {
    dlnorm(
      x,
      meanlog = -0.6094,
      sdlog = 0.2
    )
  },
  
  beta2 = function(x) {
    dlnorm(
      x,
      meanlog = -0.6094,
      sdlog = 0.2
    )
  },
  
  a_inv = function(x) {
    dlnorm(
      x,
      meanlog = 1.6,
      sdlog = 0.13
    )
  },
  
  gamma_inv = function(x) {
    dlnorm(
      x,
      meanlog = 1.6,
      sdlog = 0.13
    )
  },
  
  # Prior density of log10(e0), including the Jacobian
  e0 = function(x) {
    
    e0_original <- 10^x
    
    dlnorm(
      e0_original,
      meanlog = log(10),
      sdlog = 0.8
    ) *
      e0_original *
      log(10)
  },
  
  p_reported = function(x) {
    dbeta(
      x,
      shape1 = 0.28,
      shape2 = 100
    )
  },
  
  phi_inv = function(x) {
    dexp(
      x,
      rate = 10
    )
  }
)

# ============================================================
# Parameter-specific plotting ranges
# ============================================================

prior_ranges <- list(
  
  beta1 = c(0, 1.5),
  beta2 = c(0, 1.5),
  
  a_inv = c(0, 8),
  gamma_inv = c(0, 8),
  
  # log10 scale: e0 from 1 to 30,000
  e0 = c(
    log10(1),
    log10(30000)
  ),
  
  p_reported = c(0, 0.15),
  
  phi_inv = c(0, 0.2)
)

# ============================================================
# Custom diagonal panel
# Posterior histogram + prior density
# ============================================================

posterior_prior_diag <- function(
    data,
    mapping,
    ...
) {
  
  # Identify the parameter represented by the panel
  parameter <- rlang::as_label(
    mapping$x
  )
  
  # Extract finite posterior samples
  posterior_values <- data[[parameter]]
  
  posterior_values <- posterior_values[
    is.finite(posterior_values)
  ]
  
  # ----------------------------------------------------------
  # Set x-axis limits
  # ----------------------------------------------------------
  
  posterior_range <- range(
    posterior_values,
    na.rm = TRUE
  )
  
  specified_range <- prior_ranges[[parameter]]
  
  x_limits <- c(
    min(
      specified_range[1],
      posterior_range[1]
    ),
    max(
      specified_range[2],
      posterior_range[2]
    )
  )
  
  # ----------------------------------------------------------
  # Construct the prior-density grid
  # ----------------------------------------------------------
  
  grid_lower <- x_limits[1]
  grid_upper <- x_limits[2]
  
  # Avoid evaluating positive-support distributions at zero
  if (
    parameter %in%
    c(
      "beta1",
      "beta2",
      "a_inv",
      "gamma_inv",
      "p_reported",
      "phi_inv"
    )
  ) {
    grid_lower <- max(
      grid_lower,
      1e-8
    )
  }
  
  # Keep beta-distribution parameter strictly inside (0, 1)
  if (
    parameter == "p_reported"
  ) {
    grid_upper <- min(
      grid_upper,
      1 - 1e-8
    )
  }
  
  prior_grid <- data.frame(
    x = seq(
      from = grid_lower,
      to = grid_upper,
      length.out = 500
    )
  )
  
  # Evaluate the appropriate prior density
  prior_grid$density <- prior_density[[parameter]](
    prior_grid$x
  )
  
  # Remove non-finite prior-density values
  prior_grid <- prior_grid[
    is.finite(prior_grid$density),
  ]
  
  # ----------------------------------------------------------
  # Construct the diagonal panel
  # ----------------------------------------------------------
  
  ggplot(
    data,
    aes(
      x = .data[[parameter]]
    )
  ) +
    
    # Posterior histogram
    geom_histogram(
      aes(
        y = after_stat(density)
      ),
      bins = 30,
      fill = "#4C78A8",
      colour = "white",
      alpha = 0.65
    ) +
    
    # Prior density
    geom_line(
      data = prior_grid,
      aes(
        x = x,
        y = density
      ),
      inherit.aes = FALSE,
      colour = "#D62728",
      linewidth = 1
    ) +
    
    coord_cartesian(
      xlim = x_limits,
      expand = FALSE
    ) +
    
    theme_minimal(
      base_size = 10
    ) +
    
    theme(
      axis.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# Create pairs plot
# ============================================================

pair_plot <- GGally::ggpairs(
  
  posterior_pairs,
  
  # Upper triangle: posterior correlations
  upper = list(
    continuous = GGally::wrap(
      "cor",
      size = 3
    )
  ),
  
  # Lower triangle: posterior samples
  lower = list(
    continuous = GGally::wrap(
      "points",
      alpha = 0.20,
      size = 0.30,
      colour = "#4C78A8"
    )
  ),
  
  # Diagonal: posterior histogram + prior density
  diag = list(
    continuous = posterior_prior_diag
  ),
  
  # Mathematical parameter labels
  columnLabels = c(
    "beta[1]",
    "beta[2]",
    "T[L]",
    "T[I]",
    "log[10](E[0])",
    "pi",
    "phi^{-1}"
  ),
  
  labeller = "label_parsed",
  
  progress = FALSE
  
) +
  
  theme_bw(
    base_size = 11
  ) +
  
  theme(
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# ============================================================
# Display plot
# ============================================================

pair_plot













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









#

































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








# =======================================================================
# Pairs plot with posterior histograms and prior densities on diagonals
# Model: seir_fit_5
# =======================================================================

library(rstan)
library(ggplot2)
library(GGally)
library(rlang)

# ============================================================
# Extract posterior samples
# ============================================================

posterior_pairs <- as.data.frame(
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

# ============================================================
# Retain only finite posterior samples
# ============================================================

posterior_pairs <- posterior_pairs[
  apply(
    posterior_pairs,
    1,
    function(x) all(is.finite(x))
  ),
]

# ============================================================
# Prior density functions
# ============================================================

prior_density <- list(
  
  C0 = function(x) {
    dlnorm(
      x,
      meanlog = log(5),
      sdlog = 0.8
    )
  },
  
  r = function(x) {
    dnorm(
      x,
      mean = 0.2,
      sd = 0.2
    )
  },
  
  r_minus = function(x) {
    dnorm(
      x,
      mean = 0.1,
      sd = 0.1
    )
  },
  
  phi_inv = function(x) {
    dexp(
      x,
      rate = 10
    )
  }
)

# ============================================================
# Parameter-specific plotting ranges
# ============================================================

prior_ranges <- list(
  
  C0       = c(0, 150),
  r        = c(0, 1),
  r_minus  = c(0, 0.15),
  phi_inv  = c(0, 0.5)
  
)

# ============================================================
# Custom diagonal panel
# Posterior histogram + prior density
# ============================================================

posterior_prior_diag <- function(
    data,
    mapping,
    ...
) {
  
  # Identify the parameter represented in the panel
  parameter <- rlang::as_label(
    mapping$x
  )
  
  # Extract finite posterior samples
  posterior_values <- data[[parameter]]
  
  posterior_values <- posterior_values[
    is.finite(posterior_values)
  ]
  
  # ----------------------------------------------------------
  # Set x-axis limits
  # ----------------------------------------------------------
  
  posterior_range <- range(
    posterior_values,
    na.rm = TRUE
  )
  
  specified_range <- prior_ranges[[parameter]]
  
  x_limits <- c(
    min(
      specified_range[1],
      posterior_range[1]
    ),
    max(
      specified_range[2],
      posterior_range[2]
    )
  )
  
  # ----------------------------------------------------------
  # Construct the prior-density grid
  # ----------------------------------------------------------
  
  grid_lower <- x_limits[1]
  grid_upper <- x_limits[2]
  
  # Avoid evaluating positive-support distributions at zero
  if (
    parameter %in%
    c(
      "C0",
      "phi_inv"
    )
  ) {
    grid_lower <- max(
      grid_lower,
      1e-8
    )
  }
  
  prior_grid <- data.frame(
    x = seq(
      from = grid_lower,
      to = grid_upper,
      length.out = 500
    )
  )
  
  # Evaluate the appropriate prior density
  prior_grid$density <- prior_density[[parameter]](
    prior_grid$x
  )
  
  # Remove non-finite prior-density values
  prior_grid <- prior_grid[
    is.finite(prior_grid$density),
  ]
  
  # ----------------------------------------------------------
  # Construct the diagonal panel
  # ----------------------------------------------------------
  
  ggplot(
    data,
    aes(
      x = .data[[parameter]]
    )
  ) +
    
    # Posterior histogram
    geom_histogram(
      aes(
        y = after_stat(density)
      ),
      bins = 30,
      fill = "#4C78A8",
      colour = "white",
      alpha = 0.65
    ) +
    
    # Prior density
    geom_line(
      data = prior_grid,
      aes(
        x = x,
        y = density
      ),
      inherit.aes = FALSE,
      colour = "#D62728",
      linewidth = 1
    ) +
    
    coord_cartesian(
      xlim = x_limits,
      expand = FALSE
    ) +
    
    theme_minimal(
      base_size = 10
    ) +
    
    theme(
      axis.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# Create pairs plot
# ============================================================

pair_plot <- GGally::ggpairs(
  
  posterior_pairs,
  
  # Upper triangle: posterior correlations
  upper = list(
    continuous = GGally::wrap(
      "cor",
      size = 4
    )
  ),
  
  # Lower triangle: posterior samples
  lower = list(
    continuous = GGally::wrap(
      "points",
      alpha = 0.25,
      size = 0.4,
      colour = "#4C78A8"
    )
  ),
  
  # Diagonal: posterior histogram + prior density
  diag = list(
    continuous = posterior_prior_diag
  ),
  
  # Mathematical parameter labels
  columnLabels = c(
    "C[0]",
    "r",
    "r['-']",
    "phi^{-1}"
  ),
  
  labeller = "label_parsed",
  
  progress = FALSE
  
) +
  
  theme_bw(
    base_size = 12
  ) +
  
  theme(
    strip.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

# ============================================================
# Display plot
# ============================================================

pair_plot




#








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







#




 


























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
model6 <- stan_model("model", "Model_1_FreeI0.stan")

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









# =======================================================================
# Pairs plot with posterior histograms and prior densities on diagonals
# =======================================================================

library(rstan)
library(ggplot2)
library(GGally)
library(rlang)

# ============================================================
# Extract posterior samples
# ============================================================

posterior_pairs <- as.data.frame(
  rstan::extract(
    seir_fit_6,
    pars = c(
      "R_0",
      "alpha_inv",
      "gamma_inv",
      "phi_inv",
      "p_reported",
      "I0",
      "E0"
    )
  )
)

# ============================================================
# Prior density functions
# ============================================================

prior_density <- list(
  
  R_0 = function(x) {
    dlnorm(x, meanlog = 1.0, sdlog = 0.2)
  },
  
  alpha_inv = function(x) {
    dlnorm(x, meanlog = 1.6, sdlog = 0.13)
  },
  
  gamma_inv = function(x) {
    dlnorm(x, meanlog = 1.6, sdlog = 0.13)
  },
  
  phi_inv = function(x) {
    dexp(x, rate = 10)
  },
  
  p_reported = function(x) {
    dbeta(x, shape1 = 0.28, shape2 = 100)
  },
  
  I0 = function(x) {
    dlnorm(x, meanlog = log(5), sdlog = 0.8)
  },
  
  E0 = function(x) {
    dlnorm(x, meanlog = log(10), sdlog = 0.8)
  }
)

# ============================================================
# Custom diagonal panel
# Posterior histogram + prior density
# ============================================================

posterior_prior_diag <- function(data, mapping, ...) {
  
  # Get parameter name
  parameter <- rlang::as_name(mapping$x)
  
  # Posterior samples
  posterior_values <- data[[parameter]]
  
  # Remove NA / non-finite values
  posterior_values <- posterior_values[
    is.finite(posterior_values)
  ]
  
  # ----------------------------------------------------------
  # Set sensible x-axis limits
  # ----------------------------------------------------------
  
  x_upper <- max(posterior_values, na.rm = TRUE)
  
  x_limits <- switch(
    parameter,
    
    R_0 = c(0, max(6, x_upper)),
    
    alpha_inv = c(0, max(8, x_upper)),
    
    gamma_inv = c(0, max(8, x_upper)),
    
    phi_inv = c(0, max(0.5, x_upper)),
    
    p_reported = c(0, max(0.01, x_upper)),
    
    I0 = c(0, max(100, x_upper)),
    
    E0 = c(0, max(200, x_upper)),
    
    c(0, x_upper)
  )
  
  # ----------------------------------------------------------
  # Create grid for prior density
  # ----------------------------------------------------------
  
  prior_grid <- data.frame(
    x = seq(
      x_limits[1],
      x_limits[2],
      length.out = 500
    )
  )
  
  # Evaluate prior density
  prior_grid$density <- prior_density[[parameter]](
    prior_grid$x
  )
  
  # ----------------------------------------------------------
  # Create diagonal plot
  # ----------------------------------------------------------
  
  ggplot(data, aes(x = .data[[parameter]])) +
    
    # Posterior histogram
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 30,
      fill = "#4C78A8",
      colour = "white",
      alpha = 0.65
    ) +
    
    # Prior density
    geom_line(
      data = prior_grid,
      aes(
        x = x,
        y = density
      ),
      inherit.aes = FALSE,
      colour = "#D62728",
      linewidth = 1
    ) +
    
    coord_cartesian(
      xlim = x_limits,
      expand = FALSE
    ) +
    
    theme_minimal(base_size = 11) +
    
    theme(
      axis.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# Create pairs plot
# ============================================================

pair_plot <- GGally::ggpairs(
  
  posterior_pairs,
  
  # ----------------------------------------------------------
  # Upper triangle: posterior correlations
  # ----------------------------------------------------------
  
  upper = list(
    continuous = GGally::wrap(
      "cor",
      size = 4
    )
  ),
  
  # ----------------------------------------------------------
  # Lower triangle: posterior samples
  # ----------------------------------------------------------
  
  lower = list(
    continuous = GGally::wrap(
      "points",
      alpha = 0.25,
      size = 0.4,
      colour = "#4C78A8"
    )
  ),
  
  # ----------------------------------------------------------
  # Diagonal: posterior histogram + prior density
  # ----------------------------------------------------------
  
  diag = list(
    continuous = posterior_prior_diag
  ),
  
  # ----------------------------------------------------------
  # Parameter labels
  # ----------------------------------------------------------
  
  columnLabels = c(
    "R[0]",
    "T[L]",
    "T[I]",
    "phi^{-1}",
    "pi",
    "I[0]",
    "E[0]"
  ),
  
  labeller = "label_parsed",
  
  progress = FALSE
  
) +
  
  theme_bw(base_size = 12) +
  
  theme(
    strip.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

# ============================================================
# Display
# ============================================================

pair_plot






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












################################################################################
################################################################################
################################################################################
################################################################################




