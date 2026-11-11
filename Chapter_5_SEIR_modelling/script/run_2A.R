################################################################################

# Model 2A

# SEIR with Eigenvector informed initialisation - slightly different priors

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




























