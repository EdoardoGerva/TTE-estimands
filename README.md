# Model-free estimands for Target Trial Emulation.
The current repository contains the .R scripts used for the simulation studies of the project.

The scripts are generate and analyse data for setting 1 (see the relative Section of the paper). For reproducing setting 2, the data-generating process has to be changed by adding a time-dependence in the generation of the outcome at the second time point, so that:

!!sym(paste0("Y", t)) := rnorm(n, gamma0 + gamma1 * t * !!sym(paste0("A", t)) + gamma2 * !!sym(paste0("L", t)) + gamma3 * !!prev_Y, 0.5).
