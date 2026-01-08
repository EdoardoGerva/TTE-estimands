# Model-free estimands for Target Trial Emulation
Gervasoni, E.E., De Bus, L., Vansteelandt, S., Dukes, O. "On Estimands in Target Trial Emulation" ([Pre-Print](https://arxiv.org/abs/2601.03377)).

The current repository contains the .R scripts used for the simulation studies of the project.

The scripts generate and analyse data for setting 1 (see the relative Section of the paper). 

## Setting 2
For reproducing setting 2, the data-generating process has to be modified by adding a time-dependence. To do so, in the function _visit_dgm_, the generation of the outcome at the second time point has to be changed to:

`!!sym(paste0("Y", t)) := rnorm(n, gamma0 + gamma1 * t * !!sym(paste0("A", t)) + gamma2 * !!sym(paste0("L", t)) + gamma3 * !!prev_Y, 0.5)`.

For calendar-time data, in the function _calendar_dgm_, the generation of the outcome at the second time point for retained participants has been changed to:

`Y[retained, t] <- rnorm(length(retained), gamma0 + gamma1 * t * A[retained, t] + gamma2 * L[retained, t] + gamma3 * Y[retained, t - 1], 0.5)`

and the outcome for newly entered participants to:

`Y[new_ids, t] <- rnorm(num_new, gamma0 + gamma1 * t *  A[new_ids, t] + gamma2 * L[new_ids, t], 0.5)`.

The "true effect" for constructing Bias and Coverage is computed analytically.

## Setting 3
For reproducing setting 3 (Appendix C1), the data-generating process has to be changed as follows. For visit-time data, in the function _visit_dgm_, treatment and outcome at the second time point are changed to:

`!!sym(paste0("A", t)) := ifelse(!!prev_A == 1, 1, rbinom(n, 1, pnorm(beta0 + beta1 * !!sym(paste0("L", t)) + beta2 * !!prev_A + beta3 * !!prev_Y)))`
`!!sym(paste0("Y", t)) := rbinom(n, 1, pnorm(gamma0 + gamma1 * t * !!sym(paste0("A", t)) + gamma2 * !!sym(paste0("L", t)) + rnorm(n, 0, 1)))`

For calendar-time data, in the function _calendar_dgm_, the outcome at the second time point for retained participants is changes to:

`Y[retained, t] <- rbinom(length(retained), 1, expit(gamma0 + gamma1 * t * A[retained, t] + gamma2 * L[retained, t] + gamma3 * Y[retained, t - 1]))`

and the outcome for newly entered participants to:

`Y[new_ids, t] <- rbinom(num_new, 1, expit(gamma0 + gamma1 * t * A[new_ids, t] + gamma2 * L[new_ids, t]))`.

The "true effect" for constructing Bias and Coverage is computed analytically.
