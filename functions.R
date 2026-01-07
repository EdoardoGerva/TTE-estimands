################################################################################
#
# This script contains the main functions for the simulation study
#   
################################################################################



# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

expit <- function(x) exp(x) / (1 + exp(x))

# ------------------------------------------------------------------------------
# Visit-time Data-generating mechanism
# ------------------------------------------------------------------------------

visit_dgm <- function(Tau,par) {
  
  # Model coefficients
  alpha0 <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  beta0  <- par[4]
  beta1  <- par[5]
  beta2  <- par[6]
  beta3  <- par[7]
  gamma0 <- par[8]
  gamma1 <- par[9]
  gamma2 <- par[10]
  gamma3 <- par[11]
  
  data <- tibble(patient_id = 1:n) %>% mutate(
    L1 = rnorm(n, alpha0, 0.2),
    A1 = rbinom(n, 1, expit(beta0 + beta1 * L1)),
    Y1 = rnorm(n, gamma0 + gamma1 * A1 + gamma2 * L1, 0.5),
    I1 = 1
  )
  
  if (Tau > 1) {
    for (t in 2:Tau) {
      prev_L <- sym(paste0("L", t - 1))
      prev_A <- sym(paste0("A", t - 1))
      prev_Y <- sym(paste0("Y", t - 1))
      
      data <- data %>% mutate(
        !!sym(paste0("L", t)) := rnorm(n, alpha0 + alpha1 * !!prev_L + alpha2 * !!prev_A, 0.2),
        !!sym(paste0("A", t)) := ifelse(!!prev_A == 1, 1, rbinom(n, 1, expit(beta0 + beta1 * !!sym(paste0("L", t)) + beta2 * !!prev_A + beta3 * !!prev_Y))),
        !!sym(paste0("Y", t)) := rnorm(n, gamma0 + gamma1 * !!sym(paste0("A", t)) + gamma2 * !!sym(paste0("L", t)) + gamma3 * !!prev_Y, 0.5),
        !!sym(paste0("I", t)) := ifelse(!!prev_A == 1, 0, 1)
      )
    }
  }
  
  return(setDT(data))
}

# ------------------------------------------------------------------------------
# Calendar-time Data-generating mechanism
# ------------------------------------------------------------------------------

calendar_dgm <- function(Tau,par) {
  
  # Model coefficients
  alpha0 <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  beta0  <- par[4]
  beta1  <- par[5]
  beta2  <- par[6]
  beta3  <- par[7]
  gamma0 <- par[8]
  gamma1 <- par[9]
  gamma2 <- par[10]
  gamma3 <- par[11]
  
  retention_prob <-0.6
  
  A <- matrix(0, nrow = n, ncol = Tau)
  L <- matrix(0, nrow = n, ncol = Tau)
  Y <- matrix(0, nrow = n, ncol = Tau)
  S <- matrix(1, nrow = n, ncol = Tau)
  I <- matrix(0, nrow = n, ncol = Tau)
  
  L[, 1] <- rnorm(n, alpha0, 0.2)
  A[, 1] <- rbinom(n, 1, expit(beta0 + beta1 * L[, 1]))
  Y[, 1] <- rnorm(n, gamma0 + gamma1 * A[, 1] + gamma2 * L[, 1], 0.5)
  I[, 1] <- 1
  
  current_id <- n + 1
  
  for (t in 2:Tau) {
    
    S[, t] <- ifelse(S[, t - 1] == 1, rbinom(n, 1, retention_prob), 0)
    newleavers <- length(which(S[, t] == 0))-length(which(S[, t-1] == 0))
    
    retained <- which(S[, t] == 1) 
    if (length(retained) > 0) {
      L[retained, t] <- rnorm(length(retained), alpha0 + alpha1 * L[retained, t - 1] + alpha2 * A[retained, t - 1], 0.2)
      A[retained, t] <- ifelse(A[retained, t - 1] == 1, 1, 
                               rbinom(length(retained), 1, expit(beta0 + beta1 * L[retained, t] + beta2 * A[retained, t - 1] + beta3 * Y[retained, t - 1])))
      Y[retained, t] <- rnorm(length(retained), gamma0 + gamma1 * A[retained, t] + gamma2 * L[retained, t] + gamma3 * Y[retained, t - 1], 0.5)
      I[retained, t] <- ifelse(A[retained, t - 1] == 1, 0, 1) 
    }
    
    num_new <- newleavers  
    if (num_new > 0) {
      
      A <- rbind(A, matrix(0, nrow = num_new, ncol = Tau))
      L <- rbind(L, matrix(0, nrow = num_new, ncol = Tau))
      Y <- rbind(Y, matrix(0, nrow = num_new, ncol = Tau))
      S <- rbind(S, matrix(0, nrow = num_new, ncol = Tau))
      I <- rbind(I, matrix(0, nrow = num_new, ncol = Tau))
      
      new_ids <- current_id:(current_id + num_new - 1)
      S[new_ids, t] <- 1  
      L[new_ids, t] <- rnorm(num_new, alpha0 , 0.2)
      A[new_ids, t] <- rbinom(num_new, 1, expit(beta0 + beta1 * L[new_ids, t]))
      Y[new_ids, t] <- rnorm(num_new, gamma0 + gamma1 * A[new_ids, t] + gamma2 * L[new_ids, t], 0.5)
      I[new_ids, t] <- 1  
      
      current_id <- current_id + num_new
    }
  }
  
  colnames(S) <- paste0("S", 1:Tau)
  colnames(L) <- paste0("L", 1:Tau)
  colnames(A) <- paste0("A", 1:Tau)
  colnames(Y) <- paste0("Y", 1:Tau)
  colnames(I) <- paste0("I", 1:Tau)
  
  data <- data.frame(id = 1:(current_id - 1), S, L, A, Y, I)
  return(data)
}

# ------------------------------------------------------------------------------
# Conversion to long format for visit data
# ------------------------------------------------------------------------------

long_visit <- function(data, time_points) {
  
  long_data <- melt(data, id.vars = "patient_id", 
                    measure.vars = patterns("^L", "^A", "^Y", "^I"),
                    variable.name = "time",
                    value.name = c("L", "A", "Y", "I"))
  long_data[, time := as.integer(time)]
  
  setorder(long_data, patient_id, time)
  
  # Lagged variables
  long_data[, `:=`(
    L_lag = shift(L, fill = 0),
    A_lag = shift(A, fill = 0),
    Y_lag = shift(Y, fill = 0)
  ), by = patient_id]
  
  return(long_data)
}

# ------------------------------------------------------------------------------
# Conversion to long format for calendar data
# ------------------------------------------------------------------------------

long_calendar <- function(data) {
  data %>%
    pivot_longer(
      cols = -id,                     
      names_to = c(".value", "time"), 
      names_pattern = "([A-Z])([0-9]+)" 
    ) %>%
    filter(S == 1) %>%                
    mutate(time = as.integer(time)) %>% 
    arrange(id, time) %>%             
    group_by(id) %>%                  
    mutate(
      Y_lag = lag(Y),                 
      A_lag = lag(A),                 
      L_lag = lag(L),                 
      visit = row_number()            
    ) %>% 
    ungroup() %>%                     
    replace_na(list(Y_lag = 0))   
}

# ------------------------------------------------------------------------------
# Estimating equations (geex)
# ------------------------------------------------------------------------------

# baseline adjusted effect

est_fun_baselineadj_ipw <- function(data) {
  
  L1 <- data$L1
  L2 <- data$L2
  A1 <- data$A1
  A2 <- data$A2
  Y1 <- data$Y1
  Y2 <- data$Y2
  I1 <- data$I1
  I2 <- data$I2
  
  function(theta) {
    
    # propensity scores
    ps1 <- expit(theta[1]+theta[2]*L1)
    ps3 <- expit(theta[3]+theta[4]*L1+theta[5]*L2+theta[6]*Y1)
    psi_hat <- 1/2 *(((A1*Y1*I1)/ps1 - ((1-A1)*Y1*I1)/(1-ps1)) + ((1/(1-ps1))*((A2*Y2*I2)/ps3 - ((1-A2)*Y2*I2)/(1-ps3))))
    
    c(
      A1 - ps1,
      L1*(A1 - ps1),
      I2*(A2 - ps3),    
      L1*I2*(A2 - ps3),
      L2*I2*(A2 - ps3),
      Y1*I2*(A2 - ps3),
      psi_hat - theta[7]
    )
  }
}

est_fun_baselineadj_gcomp <- function(data) {
  
  L1 <- data$L1
  L2 <- data$L2
  A1 <- data$A1
  A2 <- data$A2
  Y1 <- data$Y1
  Y2 <- data$Y2
  I1 <- data$I1
  I2 <- data$I2
  
  function(theta) {
    
    m1 <- theta[1] + theta[2]*L1
    m0 <- theta[3] + theta[4]*L1
    m21 <- theta[5] + theta[6]*L1 + theta[7]*L2 + theta[8]*Y1
    m20 <- theta[9] + theta[10]*L1 + theta[11]*L2 + theta[12]*Y1
    m3 <- theta[13] + theta[14]*L1
    
    psi1 <- m1 - m0
    delta <- m21 - m20
    
    psi_hat <- 1/2 * (psi1 + m3)
    
    c(
      A1*(Y1 - m1),        
      L1*A1*(Y1 - m1),
      (1-A1)*(Y1 - m0),        
      L1*(1-A1)*(Y1 - m0),
      
      I2*A2*(Y2 - m21),
      L1 * I2*A2*(Y2 - m21),
      L2 * I2*A2*(Y2 - m21),
      Y1 * I2*A2*(Y2 - m21),
      
      I2*(1-A2)*(Y2 - m20),
      L1 * I2*(1-A2)*(Y2 - m20),
      L2 * I2*(1-A2)*(Y2 - m20),
      Y1 * I2*(1-A2)*(Y2 - m20),
      
      I2*(delta - m3),
      L1 * I2*(delta - m3),
      
      psi_hat - theta[15]
    )
  }
}

# uniformly weighted effect

est_fun_unif_ipw <- function(data) {
  
  L1 <- data$L1
  L2 <- data$L2
  A1 <- data$A1
  A2 <- data$A2
  Y1 <- data$Y1
  Y2 <- data$Y2
  I1 <- data$I1
  I2 <- data$I2
  S1 <- data$S1
  S2 <- data$S2
  
  function(theta, mu) {
    
    ps1      <- expit(theta[1] + theta[2]*L1)
    ps2_old  <- expit(theta[3] + theta[4]*L1 + theta[5]*L2 + theta[6]*Y1)
    ps2_new  <- expit(theta[7] + theta[8]*L2)
    
    visit1 <- ((A1 * Y1 / ps1) - ((1 - A1) * Y1 / (1 - ps1))) * (S1 == 1)
    
    visit2_old <- ((A2 * Y2 / ps2_old) - ((1 - A2) * Y2 / (1 - ps2_old))) * (I2 /(1-ps1)) * (S1 == 1 & S2 == 1)
    visit2_new <- ((A2 * Y2 / ps2_new) - ((1 - A2) * Y2 / (1 - ps2_new))) * (S1 == 0 & S2 == 1)
    
    visit_total <- 0.5 * (visit1 + visit2_old + visit2_new) * mu$w
    
    c(
      I1 * (A1 - ps1),
      I1 * L1 * (A1 - ps1),
      
      I2 * S1 * S2 * (A2 - ps2_old),
      I2 * S1 * S2 * L1 * (A2 - ps2_old),
      I2 * S1 * S2 * L2 * (A2 - ps2_old),
      I2 * S1 * S2 * Y1 * (A2 - ps2_old),
      
      I2 * (1 - S1) * S2 * (A2 - ps2_new),
      I2 * (1 - S1) * S2 * L2 * (A2 - ps2_new),
      
      I2 - theta[9],
      
      visit_total - theta[10]
    )
  }
}

est_fun_unif_gcomp <- function(data) {
  
  L1 <- data$L1
  L2 <- data$L2
  A1 <- data$A1
  A2 <- data$A2
  Y1 <- data$Y1
  Y2 <- data$Y2
  I1 <- data$I1
  I2 <- data$I2
  S1 <- data$S1
  
  function(theta, mu) {
    
    m1_t1 <- theta[1] + theta[2]*L1
    m0_t1 <- theta[3] + theta[4]*L1
    
    m1_t2_old <- theta[5] + theta[6]*L1 + theta[7]*L2 + theta[8]*Y1
    m0_t2_old <- theta[9] + theta[10]*L1 + theta[11]*L2 + theta[12]*Y1
    
    m1_t2_new <- theta[13] + theta[14]*L2
    m0_t2_new <- theta[15] + theta[16]*L2
    
    psi1_hat <- m1_t1 - m0_t1
    
    psi2_old_hat <- m1_t2_old - m0_t2_old
    psi2_new_hat <- m1_t2_new - m0_t2_new
    psi2_hat <- mu$w_old * psi2_old_hat + mu$w_new * psi2_new_hat
    
    psi_hat <- 0.5 * (psi1_hat + psi2_hat)
    
    c(
      
      A1 * I1 * (Y1 - m1_t1),
      A1 * I1 * L1 * (Y1 - m1_t1),
      (1 - A1) * I1 * (Y1 - m0_t1),
      (1 - A1) * I1 * L1 * (Y1 - m0_t1),
      
      A2 * I2 * S1 * (Y2 - m1_t2_old),
      A2 * I2 * S1 * L1 * (Y2 - m1_t2_old),
      A2 * I2 * S1 * L2 * (Y2 - m1_t2_old),
      A2 * I2 * S1 * Y1 * (Y2 - m1_t2_old),
      (1 - A2) * I2 * S1 * (Y2 - m0_t2_old),
      (1 - A2) * I2 * S1 * L1 * (Y2 - m0_t2_old),
      (1 - A2) * I2 * S1 * L2 * (Y2 - m0_t2_old),
      (1 - A2) * I2 * S1 * Y1 * (Y2 - m0_t2_old),
      
      A2 * I2 * (1 - S1) * (Y2 - m1_t2_new),
      A2 * I2 * (1 - S1) * L2 * (Y2 - m1_t2_new),
      (1 - A2) * I2 * (1 - S1) * (Y2 - m0_t2_new),
      (1 - A2) * I2 * (1 - S1) * L2 * (Y2 - m0_t2_new),
      
      psi_hat - theta[17]
    )
  }
}

# eligibility weighted effect

est_fun_elig_ipw <- function(data) {
  
  L1 <- data$L1
  L2 <- data$L2
  A1 <- data$A1
  A2 <- data$A2
  Y1 <- data$Y1
  Y2 <- data$Y2
  I1 <- data$I1
  I2 <- data$I2
  S1 <- data$S1
  S2 <- data$S2
  
  function(theta, mu) {
    
    ps1     <- expit(theta[1] + theta[2]*L1)                           
    ps2_old <- expit(theta[3] + theta[4]*L1 + theta[5]*L2 + theta[6]*Y1) 
    ps2_new <- expit(theta[7] + theta[8]*L2)                             
    
    scale <- 1 / (1 + mu$mu)
    visit1 <- scale * ((A1 * Y1 / ps1) - ((1 - A1) * Y1 / (1 - ps1))) * (S1 == 1)
    visit2_old <- scale * ((A2 * Y2 / ps2_old) - ((1 - A2) * Y2 / (1 - ps2_old))) * I2 * (S1 == 1 & S2 == 1)
    visit2_new <- scale * ((A2 * Y2 / ps2_new) - ((1 - A2) * Y2 / (1 - ps2_new))) * (S1 == 0 & S2 == 1)
    visit_total <- (visit1 + visit2_old + visit2_new) * mu$w
    
    c(
      I1 * (A1 - ps1),
      I1 * L1 * (A1 - ps1),
      
      I2 * S1 * S2 * (A2 - ps2_old),
      I2 * S1 * S2 * L1 * (A2 - ps2_old),
      I2 * S1 * S2 * L2 * (A2 - ps2_old),
      I2 * S1 * S2 * Y1 * (A2 - ps2_old),
      
      I2 * (1 - S1) * S2 * (A2 - ps2_new),
      I2 * (1 - S1) * S2 * L2 * (A2 - ps2_new),
      
      (S2 == 1) * (I2 - theta[9]),
      
      visit_total - theta[10]
    )
  }
}

est_fun_elig_gcomp <- function(data) {
  
  L1 <- data$L1
  L2 <- data$L2
  A1 <- data$A1
  A2 <- data$A2
  Y1 <- data$Y1
  Y2 <- data$Y2
  I1 <- data$I1
  I2 <- data$I2
  S1 <- data$S1
  
  function(theta, mu) {
    
    m1_t1 <- theta[1] + theta[2]*L1
    m0_t1 <- theta[3] + theta[4]*L1
    
    m1_t2_old <- theta[5] + theta[6]*L1 + theta[7]*L2 + theta[8]*Y1
    m0_t2_old <- theta[9] + theta[10]*L1 + theta[11]*L2 + theta[12]*Y1
    
    m1_t2_new <- theta[13] + theta[14]*L2
    m0_t2_new <- theta[15] + theta[16]*L2
    
    psi1_hat <- m1_t1 - m0_t1
    
    psi2_old_hat <- m1_t2_old - m0_t2_old
    psi2_new_hat <- m1_t2_new - m0_t2_new
    psi2_hat <- mu$w_old * psi2_old_hat + mu$w_new * psi2_new_hat
    
    psi_hat <- (psi1_hat/(1+mu$PI2) + psi2_hat*(mu$PI2/(1+mu$PI2)))
    
    c(
      A1 * I1 * (Y1 - m1_t1),
      A1 * I1 * L1 * (Y1 - m1_t1),
      (1 - A1) * I1 * (Y1 - m0_t1),
      (1 - A1) * I1 * L1 * (Y1 - m0_t1),
      
      A2 * I2 * S1 * (Y2 - m1_t2_old),
      A2 * I2 * S1 * L1 * (Y2 - m1_t2_old),
      A2 * I2 * S1 * L2 * (Y2 - m1_t2_old),
      A2 * I2 * S1 * Y1 * (Y2 - m1_t2_old),
      (1 - A2) * I2 * S1 * (Y2 - m0_t2_old),
      (1 - A2) * I2 * S1 * L1 * (Y2 - m0_t2_old),
      (1 - A2) * I2 * S1 * L2 * (Y2 - m0_t2_old),
      (1 - A2) * I2 * S1 * Y1 * (Y2 - m0_t2_old),
      
      A2 * I2 * (1 - S1) * (Y2 - m1_t2_new),
      A2 * I2 * (1 - S1) * L2 * (Y2 - m1_t2_new),
      (1 - A2) * I2 * (1 - S1) * (Y2 - m0_t2_new),
      (1 - A2) * I2 * (1 - S1) * L2 * (Y2 - m0_t2_new),
      
      I2 - theta[17],
      
      psi_hat - theta[18]
    )
  }
}
