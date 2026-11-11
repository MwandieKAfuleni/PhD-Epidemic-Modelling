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






























