#' From the likelihood grids, we infer 
#' 1) the ML parameter combination in each setting 
#' (KM, PSI comb. w. smallPSI,Beta)
#' 2) the BFs Psi/KM, Beta/KM, Beta/Psi
#' 3) For each best model, three effect size
#'    measures
#' 4) The chosen best model (MMC, KM, Beta, Psi)
#' where we choose 
#' -KM if 
#' BFMMC:=max(BF(Psi/KM),BF(Beta/KM))<=10
#' -Beta if BFMMC >10 AND BF(Beta/Psi)>10
#' -Dirac if BFMMC>10 AND BF(Beta/Psi)<0.1 
#' -MMC if BFMMC>10 AND 0.1<=BF(Beta/Psi)<=10
#' Moreover, for plotting, we write out the 
#' best fitting SFS for each model

#' Compute the expected relative SFS w. mispolarisation
#' first three arguments are parameters, last the corresponding
#' phi table
ratio_esfs_mispol <- function(alpha1,g1,e1,phi1){
  phipos <- which.min(apply(params_coal_exp,1,
                            function(v){sum(abs(v-c(g1,alpha1)))}))
  esfs <- phi1[phipos,] 
  esfs <- ((1-e1)*esfs + e1*esfs[length(esfs):1]) 
  return(esfs)
}

#' Defining the GOF measures
#' Cramer's V, Cohen's \omega, and d:[0,1] version 
#' nsfs: observed SFS/sum(SFS), p_esfs: Approx of freqs of 
#' theoretical SFS entries
gof_index <- function(nsfs,p_esfs){
  #' Compute Cohen's \omega
  omega2 <- (nsfs-p_esfs)^2
  omega2 <- omega2/p_esfs
  omega2 <- sum(omega2)
  #' Related Cramer's d for GOF chi2
  cram_d2 <- omega2/(length(p_esfs)-1)
  q <- min(p_esfs)
  relcram_d2 <- omega2*q/(1-q)
  #relat_cramd2 <- 
  return(data.frame(cohw=sqrt(omega2),
           cramd=sqrt(cram_d2),
           relcramd2 = relcram_d2))
}

#' SFS files

sfs_names <- list.files("../sfs_data/real_sfs/")
#' Species names
library(readxl)
USFS_data_EK <- read_excel("../sfs_data/USFS_data-EK.xlsx")
files_names <- USFS_data_EK$File_name
species_names <- paste(USFS_data_EK$...1,USFS_data_EK$...2)

ml_bf_eff <- NULL
sfs_list <- vector("list",length(sfs_names))
names(sfs_list) <- sfs_names

for (name1 in sfs_names){
  sfs_list[[name1]] <- vector("list",4)
  names(sfs_list[[name1]]) <- c("obs","km","beta","psi")
  #Record name of species
  species1 <- species_names[which(files_names==name1)]

  #load observed sfs
  sfs1 <- read.delim(paste0("../sfs_data/real_sfs/",name1), 
                     header=FALSE)
  #Scale observed sfs
  sfs_classes <- sfs1$V1
  sfs1 <- sfs1$V2
  sfs1 <- sfs1/sum(sfs1)
  
  sfs_list[[name1]][["obs"]] <- sfs1
for (model1 in c("beta","psi","smallpsi")){



#' We need to extract the right line of a phi table,
#' these are the corresponding parameter values
#' HAVE TO BE AS IN ..._data_analysis.R and 
n_coalp <- 20
coal_steps <- (0:n_coalp)/n_coalp 
if (model1=="beta"){coal_steps <- coal_steps + 1} #alpha runs 1->2, Psi 0->1 
rhos <- c(0,25,50)
if (model1=="smallpsi"){
n_coalp <- 20 
coal_steps <- 0.2*(0:n_coalp)/n_coalp   
rhos <- c(0,30,60)
}
params_coal_exp <- expand.grid(seq(rhos[1],rhos[2],
                                   length.out = rhos[3]+1),
                                   coal_steps)

colnames(params_coal_exp) <- c("g",switch(model1,"beta"="alpha","psi"="psi")) 


  #Load phi table
  phi.table <- read.delim(paste0("phi_tables_",model1,"/phi.table.",name1), 
                          header=FALSE)
  phi.table <- phi.table[,-1]
  
  grid1 <- read.delim(paste0("likelihood_grids_",model1,"/grid_",name1), 
                      header=FALSE, 
                      comment.char="#",colClasses = "numeric")
  #Best fit overall
  best_params <- grid1[which.max(grid1$V5),c(2,3,4,5)]
  #Best KM+exp fit
  km_lines <- which(grid1$V2==switch(model1,"beta"=2,"psi"=0,
                                     "smallpsi"=0))
  grid_km <- grid1[km_lines,]
  best_params_km <- grid_km[which.max(grid_km$V5),c(2,3,4,5)]
  
  
  erat_best <- ratio_esfs_mispol(alpha1 = as.numeric(best_params[1]),
                                    g1 = as.numeric(best_params[2]),
                                    e1 = as.numeric(best_params[3]),
                                    phi1 =   phi.table) 
  erat_best_KM <- ratio_esfs_mispol(alpha1 = as.numeric(best_params_km[1]),
                                    g1 = as.numeric(best_params_km[2]),
                                    e1 = as.numeric(best_params_km[3]),
                                    phi1 =   phi.table)    
  
  if (model1=="beta"){ml_beta <- best_params
          sfs_list[[name1]][["beta"]] <- erat_best
          ml_km <- best_params_km
          sfs_list[[name1]][["km"]] <- erat_best_KM}

  if (model1=="psi"){
    ml_psi <- best_params
    sfs_list[[name1]][["psi"]] <- erat_best
    if (best_params_km[4]>ml_km[4]){
      ml_km <- best_params_km
      sfs_list[[name1]][["km"]] <- erat_best_KM
    }
  }
  
if (model1=="smallpsi"){
  if (best_params[4]>ml_psi[4]){
    ml_psi <- best_params
    sfs_list[[name1]][["psi"]] <- erat_best
  }
  if (best_params_km[4]>ml_km[4]){
    ml_km <- best_params_km
    sfs_list[[name1]][["km"]] <- erat_best_KM
  }
}


} 
if (max(ml_beta[4],ml_psi[4])-ml_km[4]<=log(10)){
  modelselect1 <- "KM"} else {
    modelselect1 <- "MMC"
    if (ml_beta[4]-ml_psi[4] > log(10)){
      modelselect1 <- "BETA"}
    if (ml_beta[4]-ml_psi[4] < -log(10)){
      modelselect1 <- "PSI"}}
temp1 <- data.frame(name1,species1,ml_beta,
           gof_index(sfs1,sfs_list[[name1]][["beta"]]),
           ml_psi,
           gof_index(sfs1,sfs_list[[name1]][["psi"]]),
           ml_km,
           gof_index(sfs1,sfs_list[[name1]][["km"]]),
           max(ml_beta[4],ml_psi[4])-ml_km[4],
           ml_beta[4]-ml_psi[4],modelselect1)  
ml_bf_eff <- rbind(ml_bf_eff,temp1)
}
mlcols <- c("coalp","growth","misid","logML",
            "cohenw","cramerv","relcramd")
colnames(ml_bf_eff) <- c("file","species",
                         paste0("betaml_",mlcols),
                         paste0("psiml_",mlcols),
                         paste0("kmml_",mlcols),
                         "logBF_mmcvskm",
                         "logBF_bvspsi",
                         "selected_model")

write.table(ml_bf_eff,
            file="res_summaries/mlres_BF_effectsize.txt",
            quote = FALSE,row.names = FALSE,
            col.names = TRUE,sep=";")
save(ml_bf_eff,sfs_list,file="res_summaries/mlres_BF_effectsize.RData")