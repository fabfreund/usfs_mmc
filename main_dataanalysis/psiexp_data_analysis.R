#' R libraries for parallelization and reading in .xlsx 
library(doMC)
library(readxl)

 
#' Parameters to be manually set

#' Phi table params, psi presented differently because we parallelize the phi table computation
n_psis <- 20 
psi_steps <- (0:n_psis)/n_psis   
rhos <- c(0,25,50)
misi <- c(0,.15,15)
prec1 <- 1000 #Can be put lower for large data sets

#' Number of cores for parallelization  
mc1 <- 3
registerDoMC(mc1)

#' Absolute path to current folder
home_dir <- getwd()


#' Extract file names of the SFS 
sfs_names <- list.files(path = "../sfs_data/real_sfs/")



#' compute_phi regulates whether phi table need to be computed 
#' (with params as above). Not necessary by default for replicating the 
#' results, since the tables are also in the repo
compute_phi <- rep(FALSE,length(sfs_names))
names(compute_phi) <- sfs_names


#' From here downwards, this runs automatically

#' Read in passport information, including 
#' Sneq: number of third allele outgroup sites 


USFS_data_EK <- read_excel("../sfs_data/USFS_data-EK.xlsx")
#' #SNP where the outgroup has a different allele (3alleles)`
Sneq <- USFS_data_EK$`#SNP outgroup has a different allele (3alleles)`
names(Sneq) <- USFS_data_EK$File_name  

for (name1 in sfs_names){
sfs1 <- read.table(paste0("../sfs_data/real_sfs/",name1))
s <- max(sfs1[,1]) + 1
sfs1 <- sfs1[,2]

#' Compute phi table
if (compute_phi[name1]){  
foreach (i = 1:(n_psis+1)) %dopar% {
  com_phi <- "./../main_sim_inf_tool/Psi/MMC-Phi-LookUpTable/MMC-Phi-LookUpTable.out -sampleSize"
  com_phi <- paste(com_phi,s,"-MinPsi",psi_steps[i],
                   "-MaxPsi",psi_steps[i]) 
  com_phi <- paste(com_phi,"-NoStepsPsi 0 -MinRho",rhos[1],"-MaxRho",rhos[2],"-NoStepsRho",rhos[3])
  com_phi <- paste0(com_phi," -precPhi ",prec1," -OUT temp",i,"__")
  system(com_phi)
  system(paste0("mv temp",i,"__* temp",i,".txt"))
}

system(paste0("cat temp1.txt > phi.table.",name1))
for (i in 2:(n_psis+1)){
  system(paste0("cat temp",i,".txt >> phi.table.",name1))
}

system("rm temp*")
system(paste0("mv phi.table.",name1," phi_tables_psi/"))
}

#' Extract Sneq
Sneq1 <- unname(Sneq[name1])
#' Extract SFS, only proceed if all SFS classes are filled (lazy, one could set missing entries to 0) 
if (length(sfs1) == s-1){write.table(t(sfs1),file="tempsfs.txt",
                                     sep = "\t",row.names = FALSE,col.names = FALSE)} else {next}
  
com_grid <- "./../main_sim_inf_tool/Psi/MMC-MaxLikelihoodInference-GridSearch/"
com_grid <- paste0(com_grid,"MMC-MaxLikelihoodInference-GridSearch.out -SFS")
sfs_name <- paste0(home_dir,"/tempsfs.txt") 
com_grid_v <- paste(com_grid,sfs_name,
                    "-Phi",paste0(home_dir,"/phi_tables_psi/phi.table.",name1),
                    "-minPsi 0 -maxPsi 1 -noStepsPsi",n_psis) 
com_grid_v <- paste(com_grid_v,"-minRho",rhos[1],"-maxRho",rhos[2],"-noStepsRho",rhos[3]) 
com_grid_v <- paste(com_grid_v,"-minMisIdent",misi[1],"-maxMisIdent",misi[2],"-noStepsMisIdent",misi[3])
com_grid_v <- paste(com_grid_v,"-noIncongruentSites",Sneq1,"-printGrid")
system(com_grid_v)

system(paste0("mv tempsfs-Lump*Grid.txt likelihood_grids_psi/grid_",name1))
system(paste0("mv tempsfs-Lump*MLEstimates.txt ML_est_psi/ML_",name1))

system(paste0("rm tempsfs*")) 
}


