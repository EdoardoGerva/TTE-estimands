################################################################################
# Eligibility Weighted Effect - IPW
#
# This script:
#   - Simulates longitudinal data in a calendar-time design
#   - Estimates the Eligibility Weighted Effect via IPW
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
  
  # Nuisance estimation
  ps1_m <- glm(A1 ~ L1, data[data$I1==1,],family="binomial")
  ps1 <- predict(ps1_m, data, type = "response")
  ps2 <- 1 - ps1
  ps3_m <- glm(A2 ~ L1 + L2 + Y1, data[data$I2==1&data$S1==1&data$S2==1,], family = "binomial")
  ps3 <- predict(ps3_m, data, type = "response")
  ps4_m <- glm(A2 ~ L2,  data[data$I2==1&data$S1==0&data$S2==1,], family = "binomial")
  ps4 <- predict(ps4_m, data, type = "response")
  
  mu <- mean(data[data$S2==1,'I2'])
  
  visit1 <- (1/(1+mu))*((data$A1 * data$Y1 / ps1) - ((1 - data$A1) * data$Y1 / (1 - ps1)))[data$S1==1]
  visit2_old <- (1/(1+mu))*(((data$A2 * data$Y2 / ps3) - ((1 - data$A2) * data$Y2 / (1 - ps3))) * (data$I2))[(data$S1==1)&(data$S2==1)]
  visit2_new <- (1/(1+mu))*(((data$A2 * data$Y2 / ps4) - ((1 - data$A2) * data$Y2 / (1 - ps4))))[(data$S1==0)&(data$S2==1)]
  visit2 <- c(visit2_old,visit2_new)
  psi_hat <- (mean(visit1)+mean(visit2))
  
  theta_hat <- c(ps1_m$coefficients,ps3_m$coefficients,ps4_m$coefficients,mu,psi_hat)
  
  # m-estimation
  res <- m_estimate(
    estFUN = est_fun_elig_ipw,
    data = data,
    inner_args = list(mu = list(mu=mu,w=dim(data)[1]/n)),
    compute_roots = FALSE,
    roots = theta_hat
  )
  
  # Save the results
  list(gamma1_lm = naive_model$coefficients['A'], se_lm = summary(naive_model)$coefficients[2,2], gamma1_hat = psi_hat, se = sqrt(vcov(res)[length(theta_hat),length(theta_hat)]))
  
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

