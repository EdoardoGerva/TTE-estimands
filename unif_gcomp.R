################################################################################
# Uniformly Weighted Effect - G-computation
#
# This script:
#   - Simulates longitudinal data in a calendar-time design
#   - Estimates the Uniformly Weighted Effect via G-computation
#   - Compares results to a pooled linear regression benchmark
#
################################################################################

# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

library(geex)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(parallel)
library(doSNOW)

# ------------------------------------------------------------------------------
# Call function file
# ------------------------------------------------------------------------------


source("functions.R")

# ------------------------------------------------------------------------------
# Set parameters
# ------------------------------------------------------------------------------

set.seed(123)

n <- 1000
Tau <- 2
N <- 1000

par <- c(
  alpha0 = 0.5,
  alpha1 = 1,
  alpha2 = 1,
  beta0  = 0.5,
  beta1  = -1,
  beta2  = -1,
  beta3  = -1,
  gamma0 = 0.5,
  gamma1 = 1,
  gamma2 = 1,
  gamma3 = 1
)

# ------------------------------------------------------------------------------
# Simulation
# ------------------------------------------------------------------------------

cores <- parallel::detectCores()-2
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)

pb <- txtProgressBar(min=1, max=N, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

system.time(results <- foreach(i=1:N, .options.snow=opts, .packages = c('dplyr','data.table','tidyr','stats','geex'), .combine = rbind) %dopar% {
  
  data <- calendar_dgm(Tau,par)
  
  # pooled lm
  data_long <- long_calendar(data)
  naive_model <- geepack::geeglm(Y ~ A+L+Y_lag, data_long[data_long$I==1,], family = gaussian, id=id)
  
  n_old <- sum(data$I2 == 1 & data$S1 == 1)
  n_new <- sum(data$I2 == 1 & data$S1 == 0)
  w_old <- n_old / (n_old + n_new)
  w_new <- n_new / (n_old + n_new)
  
  # m-estimation  
  res <- m_estimate(
    estFUN = est_fun_unif_gcomp,
    data = data,
    root_control = setup_root_control(start = rep(0, 17)),
    inner_args = list(mu = list(w_old = w_old, w_new = w_new))
  )
  
  # Save the results
  list(gamma1_lm = naive_model$coefficients['A'], se_lm = summary(naive_model)$coefficients[2,2], gamma1_hat = coef(res)[17], se = sqrt(vcov(res)[17,17]))  
})

close(pb)
stopCluster(cl)

results_df <- as.data.frame(results,row.names = F)
results_df <- data.frame(lapply(results_df, unlist))
results <- results_df

# ------------------------------------------------------------------------------
# Results
# ------------------------------------------------------------------------------

head(results)
summary(results$se)

mean(results$se) / sd(results$gamma1_hat)

mean(results$se)
sd(results$gamma1_hat)

mean(results$se_lm)
sd(results$gamma1_lm)

effect <- 1

mean(results$gamma1_hat) - effect
mean(results$gamma1_lm) - effect


mean(ifelse(results$gamma1_hat - 1.96 * results$se <= effect &
              results$gamma1_hat + 1.96 * results$se >= effect, 1, 0))

mean(ifelse(results$gamma1_lm - 1.96 * results$se_lm <= effect &
              results$gamma1_lm + 1.96 * results$se_lm >= effect, 1, 0))

