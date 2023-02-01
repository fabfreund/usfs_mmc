library(doMC)
library(parallel)
        
#what coalescent model to run 
run_psi <- TRUE
run_beta <- TRUE

mc1 <- 6#7
## Full paths to used C++ scripts and output directories
home_dir <- c(getwd(),"/")
com_phi_b <- "./..Beta/MMC-Phi-LookUpTable/MMC-Phi-LookUpTable.out"
com_phi_d <- "./../Psi/MMC-Phi-LookUpTable/MMC-Phi-LookUpTable.out"
out1_short <- "temp"
out1 <- paste0(home_dir,out1_short)
com_sim_b <- "./../Beta/MMC-CoalescentSimulator/MMC-CoalescentSimulator.out"
com_sim_d <- "./../Psi/MMC-CoalescentSimulator//MMC-CoalescentSimulator.out"
out2_short <- "tempsfs"
out2 <- paste0(home_dir,out2_short)


#' sample size
s <- 20
#' Inference parameter Phi
psi_inf <-  c("min"=0,"max"=1,
              "steps"=20)
beta_inf <- c("min"=1,"max"=2,
              "steps"=20)
#Adjust to the coarser grid in data analysis
r_inf_psi <- c("min"=0,"max"=25,"steps"=50)

setwd(home_dir)

#Psi coalescent phi table

registerDoMC(mc1)
if (run_psi){
psi_vals <- seq(from=psi_inf[1],to=psi_inf[2],length.out=(psi_inf[3]+1))   

foreach (i = seq(along=psi_vals)) %dopar% {
  com_phi <- "./../Psi/MMC-Phi-LookUpTable/MMC-Phi-LookUpTable.out -sampleSize"
  com_phi <- paste(com_phi,s,"-MinPsi",psi_vals[i],
                   "-MaxPsi",psi_vals[i]) 
  com_phi <- paste(com_phi,"-NoStepsPsi 0 -MinRho",r_inf_psi[1],"-MaxRho",r_inf_psi[2],"-NoStepsRho",r_inf_psi[3])
  com_phi <- paste0(com_phi," -OUT temp",i,"_")
  system(com_phi)
  system(paste0("mv temp",i,"_* temp",i,".txt"))
}

system(paste0("cat temp1.txt > psi.phi.table_n",
              s,add2))
for (i in 2:length(psi_vals)){
  system(paste0("cat temp",i,".txt >> psi.phi.table_n",
                s,add2))
}

system("rm temp*")}

#Beta coalescent phi table
if (run_beta){
beta_vals <- seq(from=beta_inf[1],to=beta_inf[2],
                length.out=(beta_inf[3]+1))   

foreach (i = seq(along=beta_vals)) %dopar% {
  com_phi <- "./../Beta/MMC-Phi-LookUpTable/MMC-Phi-LookUpTable.out -sampleSize"
  com_phi <- paste(com_phi,s,"-MinAlpha",beta_vals[i],
                   "-MaxAlpha",beta_vals[i]) 
  com_phi <- paste(com_phi,"-NoStepsAlpha 0 -MinRho",r_inf_psi[1],"-MaxRho",r_inf_psi[2],"-NoStepsRho",r_inf_psi[3])
  com_phi <- paste0(com_phi," -OUT temp",i,"_")
  system(com_phi)
  system(paste0("mv temp",i,"_* temp",i,".txt"))
}

system(paste0("cat temp1.txt > beta.phi.table_n",
              s))
for (i in 2:length(beta_vals)){
  system(paste0("cat temp",i,".txt >> beta.phi.table_n",
                s))
}

system("rm temp*")}