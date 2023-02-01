#' R libraries for parallelization and reading in .xlsx 
library(doMC)
library(readxl)

 
#' Parameters to be manually set

#' Phi table params, alpha presented differently because we parallelize the phi table computation
n_alphas <- 20
alpha_steps <- 1 + (0:n_alphas)/n_alphas   
rhos <- c(0,25,50)
misi <- c(0,.15,15)

#' Number of cores for parallelization  
mc1 <- 3#Adjust as fits for your system
registerDoMC(mc1)

#' Absolute path to current folder
home_dir <- getwd()


#' Extract file names of the SFS 
sfs_names <- list.files(path = "../sfs_data/real_sfs/")



#' compute_phi regulates whether phi table need to be computed 
#' (with params as above)
#' We have all phi tables precomputed, so this is by default FALSE
compute_phi <- rep(FALSE,length(sfs_names))
names(compute_phi) <- sfs_names


#' From here downwards, this runs automatically

#' Read in passport information, including 
#' Sneq: number of third allele outgroup sites 

#USFS_data_EK <- read_excel("Data_USFS-1219/USFS_data-EK.xlsx")
USFS_data_EK <- read_excel("../sfs_data/USFS_data-EK.xlsx")
#Sneq <- USFS_data_EK$`#SNP where the ancester has a different allele (3alleles)`
Sneq <- USFS_data_EK$`#SNP outgroup has a different allele (3alleles)`
names(Sneq) <- USFS_data_EK$File_name  

for (name1 in sfs_names){
sfs1 <- read.table(paste0("../sfs_data/real_sfs/",name1))
s <- max(sfs1[,1]) + 1
sfs1 <- sfs1[,2]

#' Compute phi table
if (compute_phi[name1]){  
foreach (i = 1:(n_alphas+1)) %dopar% {
  com_phi <- "./../main_sim_inf_tool/Beta/MMC-Phi-LookUpTable/MMC-Phi-LookUpTable.out -sampleSize"
  com_phi <- paste(com_phi,s,"-MinAlpha",alpha_steps[i],
                   "-MaxAlpha",alpha_steps[i]) 
  com_phi <- paste(com_phi,"-NoStepsAlpha 0 -MinRho",rhos[1],"-MaxRho",rhos[2],"-NoStepsRho",rhos[3])
  com_phi <- paste0(com_phi," -OUT temp",i,"__")
  system(com_phi)
  system(paste0("mv temp",i,"__* temp",i,".txt"))
}

system(paste0("cat temp1.txt > phi.table.",name1))
for (i in 2:(n_alphas+1)){
  system(paste0("cat temp",i,".txt >> phi.table.",name1))
}

system("rm temp*")
system(paste0("mv phi.table.",name1," phi_tables_beta/"))
}

#' Extract Sneq
Sneq1 <- unname(Sneq[name1])
#' Extract SFS, only proceed if all SFS classes are filled (lazy, one could set missing entries to 0) 
if (length(sfs1) == s-1){write.table(t(sfs1),file="tempsfs.txt",
                                     sep = "\t",row.names = FALSE,col.names = FALSE)} else {next}
  
com_grid <- "./../main_sim_inf_tool/Beta/MMC-MaxLikelihoodInference-GridSearch/"
com_grid <- paste0(com_grid,"MMC-MaxLikelihoodInference-GridSearch.out -SFS")
sfs_name <- paste0(home_dir,"/tempsfs.txt") 
com_grid_v <- paste(com_grid,sfs_name,
                    "-Phi",paste0(home_dir,"/phi_tables_beta/phi.table.",name1),
                    "-minAlpha 1 -maxAlpha 2 -noStepsAlpha",n_alphas) 
com_grid_v <- paste(com_grid_v,"-minRho",rhos[1],"-maxRho",rhos[2],"-noStepsRho",rhos[3]) 
com_grid_v <- paste(com_grid_v,"-minMisIdent",misi[1],"-maxMisIdent",misi[2],"-noStepsMisIdent",misi[3])
com_grid_v <- paste(com_grid_v,"-noIncongruentSites",Sneq1,"-printGrid")
system(com_grid_v)

system(paste0("mv tempsfs-Lump*Grid.txt likelihood_grids_beta/grid_",name1))
system(paste0("mv tempsfs-Lump*MLEstimates.txt ML_est_beta/ML_",name1))

system(paste0("rm tempsfs*")) 
}


