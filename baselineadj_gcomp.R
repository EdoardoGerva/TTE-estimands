################################################################################
# Baseline-Adjusted Effect - G-computation
#
# This script:
#   - Simulates longitudinal data in a visit-time design
#   - Estimates the Baseline-Adjusted Effect via G-computation
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
  
  data <- visit_dgm(Tau,par)
  
  # pooled lm
  data_long <- long_visit(data, Tau)
  naive_model <- geepack::geeglm(Y ~ A+L+Y_lag, data_long[data_long$I==1,], family = gaussian, id=patient_id)
  
  # Nuisance estimation
  m1 <- lm(Y1 ~ L1, data %>% filter(A1==1)) 
  m0 <- lm(Y1 ~ L1, data %>% filter(A1==0))
  
  y1 <- predict(m1, data, type = "response") 
  y0 <- predict(m0, data, type = "response")
  
  psi1 <- y1-y0
  
  m21 <- lm(Y2 ~ L1 + L2 + Y1, data %>% filter(A2==1) %>% filter(I2==1)) 
  m20 <- lm(Y2 ~ L1 + L2 + Y1, data %>% filter(A2==0) %>% filter(I2==1))
  
  y21 <- predict(m21, data, type = "response")
  y20 <- predict(m20, data, type = "response")
  
  data$delta <- y21-y20  
  
  m3 <- lm(delta ~ L1, data %>% filter(I2 == 1)) 
  psi2 <- predict(m3, data, type = "response")
  
  psi_hat <- (1/2) * mean(psi1+psi2)
  
  theta_hat <- c(m1$coefficients,m0$coefficients,m21$coefficients,m20$coefficients,m3$coefficients,psi_hat)
  
  # m-estimation
  res <- m_estimate(
    estFUN = est_fun_baselineadj_gcomp,
    data = data,
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

