#' Assess the errors of BF > 10 criterion
#' 
#' 
#' Record type I and II error
#' Test this on 2000 randomly drawn sims 
#' from all models


#' Arguments: 1=s
#args1 <- commandArgs(TRUE)
s <- 20

#' Preloadings
library(parallel)
source("inference_fun.R")

#' Threshold
thresh <- log(10)
#' Number of sims to evaluate errors
#' (for each error)

check_sims <- 2000#
mods1 <-  c("KM","BETA","PSI")


#' Extract sims
#' Functions for this

extract_sims <- function(model1,pick1,nextract){
  sims <- switch(model1,"KM"=sims_b[pick1],
                 "BETA"=sims_b[pick1],
                 "PSI"=sims_d[pick1])
  sneq <- switch(model1,"KM"=sneq_b[pick1],
                 "BETA"=sneq_b[pick1],
                 "PSI"=sneq_d[pick1])
  no_sims_per_sc <- nrow(sims[[1]])
  picks0 <- as.vector(rmultinom(1,nextract,
                                rep(1,length(sims))))
  picks1 <- sapply(picks0,function(x){
    sample(no_sims_per_sc,x,
           replace = FALSE)})
  out1 <- NULL
  out2 <- NULL
  for (i in seq(along=picks0)){
    out1 <- rbind(out1,sims[[i]][picks1[[i]],])
    out2 <- c(out2,sneq[[i]][picks1[[i]]])
  }
  return(list(sims=out1,sneq=out2))
}


#' Define function that moves the scenarios from 
#' one misid rate to the next. f1_0,f1 
#' have both 4 misid 
add_e <- function(v1,v2){c(v1,v1+v2,v1+2*v2,
                           v1+3*v2)}

#' extract grid sims


source("param_set_psi_005.R")


#' extract non-grid sims

load(paste0("sims_errortest_psi005_s",s,".RData"))
psi_pick <- list("d0.005"=which(sapply(params_dirac_exp$dirac,
                                      function(x,a){isTRUE(all.equal(x,a))},
                                      0.005)))
howmanypsi <- nrow(params_dirac_exp)

#KM_beta_e <- add_e(KM_beta,howmanybeta)

psi_pick_e <- lapply(psi_pick,add_e,howmanypsi)
picks_list <- c(psi_pick_e)

sims_list0 <- vector("list",length(psi_pick_e))
names(sims_list0) <- c(names(psi_pick_e))
mod_vec0 <- "PSI"
for (k in seq(along=sims_list0)){
  sims_list0[[k]] <- extract_sims(mod_vec0[k],picks_list[[k]],
                                  check_sims)}
mod_vec <- mod_vec0
sims_list <- sims_list0

cat("SIMS DONE \n")
#' Common inference parameters (i.e. grid parameters and misiden params)
r_inf <- c("min"=0,"max"=25,"steps"=50)
misiden_inf <- c("min"=0.001,"max"=0.201,"steps"=20)
#' Set individual inference params and phi tables
coal_inf <- list("BETA"=c("min"=1,"max"=2,"steps"=20),
                "PSI"=c("min"=0,"max"=1,"steps"=20),
                "KM"=c("min"=2,"max"=2,"steps"=0))
phibeta <- read.delim(paste0("beta.phi.table_n",s,
                             "_f1"), 
                        header=FALSE)
phipsi <- read.delim(paste0("psi.phi.table_n",s,
                            "_f1"), 
                     header=FALSE)
phibeta <- phibeta[,-1]
phipsi <- phipsi[,-1]
phi_list <- list(
               "BETA"=phibeta,
               #growth has one more step than below
               #but there's a hidden +1-1 to get to
               #start of KM rows
"KM"=phibeta[(nrow(phibeta)-r_inf["steps"]):nrow(phibeta),],
  "PSI"=phipsi)

params_phi_list <- list("BETA"=expand.grid(
  r=seq(r_inf[1],r_inf[2],length.out = r_inf[3]+1),
  coalp=seq(coal_inf[["BETA"]][1],
            coal_inf[["BETA"]][2],
      length.out = coal_inf[["BETA"]][3]+1)),
  "PSI"=expand.grid(
    r=seq(r_inf[1],r_inf[2],length.out = r_inf[3]+1),
    coalp=seq(coal_inf[["PSI"]][1],
              coal_inf[["PSI"]][2],
              length.out = coal_inf[["PSI"]][3]+1)),
  "KM"=expand.grid(
    r=seq(r_inf[1],r_inf[2],length.out = r_inf[3]+1),
    coalp=seq(coal_inf[["KM"]][1],
              coal_inf[["KM"]][2],
              length.out = coal_inf[["KM"]][3]+1)))



#' Nested list, first entry simulated model
#' second entry inferred model

#BF_list <- vector("list",length(picks_list))
#names(BF_list) <- names(picks_list)


inferBF_f <- function(sims,sneq){
  out1 <- vector("list",3)
  names(out1) <- mods1
  misid <- seq(misiden_inf[1],
            misiden_inf[2],
            length.out = misiden_inf[3]+1)
  for (mod1 in mods1){
  out1[[mod1]] <- infparams_cpslk(
                        obs_data = sims,
                        sneq_v = sneq,
              phi_tab2 = phi_list[[mod1]],
                        misiden1 = misid,
    params_phi2 = params_phi_list[[mod1]])}
  BFs <- matrix(-1,ncol=2,nrow=nrow(sims))
  colnames(BFs) <- paste0("lgBF_",
                    c("MMCvsKM","BvsPSI"))
  aux1 <- out1[["BETA"]]["logPsL",]-out1[["KM"]]["logPsL",]
  aux1 <- cbind(aux1,out1[["PSI"]]["logPsL",]-out1[["KM"]]["logPsL",])
  BFs[,1] <- apply(aux1,1,max)
  BFs[,2] <- out1[["BETA"]]["logPsL",]-out1[["PSI"]]["logPsL",]
  return(list(BF=BFs,MLestim=out1))
  }


#for (sim1 in names(BF_list)){
#BF_list[[sim1]] <- 

BF_aux <- function(sim2){out1 <- 
  inferBF_f(sims_list[[sim2]]$sims,
            sims_list[[sim2]]$sneq)
cat("BF done",sim2,"\n")
return(out1)
}

BF_list <- mclapply(names(sims_list),BF_aux,mc.cores = mc1)

names(BF_list) <- names(sims_list)

infer_model <- function(bf1,
                        thr1=log(10),
                        thr2=log(10)){
if (bf1[1]<=thr1){out1 <- "KM"} else {
  out1 <- "MMC"
  if (bf1[2] < -thr2){out1 <- "PSI"}
  if (bf1[2] > thr2){out1 <- "BETA"}
}   
return(out1)
}

param_misid <- function(i,k){
  if (mod_vec[k]=="KM"){
    out1 <- switch(sel_mods[[k]][i],"BETA"=
              BF_list[[k]][["MLestim"]][["BETA"]]["coalp",i],
          "PSI"=BF_list[[k]][["MLestim"]][["PSI"]]["coalp",i],
          "MMC"=BF_list[[k]][["MLestim"]][["BETA"]]["coalp",i])
  }
  if (mod_vec[k]=="BETA"){
    out1 <- switch(sel_mods[[k]][i],"KM"=0,
                   "PSI"=BF_list[[k]][["MLestim"]][["PSI"]]["coalp",i],
                   "MMC"=NA)
  }
  if (mod_vec[k]=="PSI"){
    out1 <- switch(sel_mods[[k]][i],"KM"=2,
                   "BETA"=BF_list[[k]][["MLestim"]][["BETA"]]["coalp",i],
                   "MMC"=NA)
  }
  return(out1)
}


sel_mods <- vector("list",length(BF_list))
names(sel_mods) <- names(BF_list)
param_misspec <- vector("list",length(BF_list))
names(param_misspec) <- names(BF_list)

for (k in seq(along=sel_mods)){
sel_mods[[k]] <- apply(BF_list[[k]]$BF,
                          1,infer_model)
wrong_mod <- which(sel_mods[[k]]!=mod_vec[k])

#misparam_sim[[sim1]] <- sims_list[[sim1]]$chosen
param_misspec[[k]] <- sapply(wrong_mod,param_misid,k)
}
save(sims_list,BF_list,sel_mods,param_misspec,
     file=paste0("crossmodBF_psi005_s",s,".RData"))
