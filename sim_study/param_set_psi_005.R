#' Small script for parameters, sourced later

## Full paths to used C++ scripts and output directories
home_dir <- c(getwd(),"/")
com_phi_b <- "./../main_sim_inf_tool/Beta/MMC-Phi-LookUpTable/MMC-Phi-LookUpTable.out"
com_phi_d <- "./../main_sim_inf_tool/Psi/MMC-Phi-LookUpTable/MMC-Phi-LookUpTable.out"
out1_short <- "temp"
out1 <- paste0(home_dir,out1_short)
com_sim_b <- "./../main_sim_inf_tool/Beta/MMC-CoalescentSimulator/MMC-CoalescentSimulator.out"
com_sim_d <- "./../main_sim_inf_tool/Psi/MMC-CoalescentSimulator//MMC-CoalescentSimulator.out"
out2_short <- "tempsfs"
out2 <- paste0(home_dir,out2_short)

#' Libraries

library(parallel)
library(doMC)


#' Parallel computing: how many cores?
mc1 <- 6

#' Simulation parameters
r1 <- c(0.25,2.25,11.25) #Growth rates
n_loci <- 100 #Number independent loci
coal_dirac <- c(0.005)
params_dirac_exp <- expand.grid(r=r1,dirac=coal_dirac)
misiden <- c(0,0.015,0.045,0.095) #Mispolarisation probability
n_simul <- 500 #How many sims per parameter triple
n_mut <- 50