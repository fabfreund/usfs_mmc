#' Small script for parameters, sourced later

## Full paths to used C++ scripts and output directories
home_dir <- paste0(getwd(),"/")
com_phi_b <- "./../Beta/MMC-Phi-LookUpTable/MMC-Phi-LookUpTable.out"
com_phi_d <- "./../Psi/MMC-Phi-LookUpTable/MMC-Phi-LookUpTable.out"
out1_short <- "temp"
out1 <- paste0(home_dir,out1_short)
com_sim_b <- "./../Beta/MMC-CoalescentSimulator/MMC-CoalescentSimulator.out"
com_sim_d <- "./../Psi/MMC-CoalescentSimulator//MMC-CoalescentSimulator.out"
out2_short <- "tempsfs"
out2 <- paste0(home_dir,out2_short)

#' Libraries

library(parallel)
library(doMC)


#' Parallel computing: how many cores?
mc1 <- 6 #to adjust

#' Simulation parameters
r1 <- c(0,0.5,1,10) #Growth rates
n_loci <- 100 #Number independent loci
coal_beta <- seq(1,2,.05) #Beta coalescent params
coal_dirac <- seq(0.05,0.95,.05)
params_beta_exp <- expand.grid(r=r1,beta=coal_beta)
params_dirac_exp <- expand.grid(r=r1,dirac=coal_dirac)
misiden <- c(0,0.01,0.05,0.1) #Mispolarisation probability
n_simul <- 500 #How many sims per parameter triple
n_mut <- 50