#' Arguments: 1=s, 2="beta"/"psi"
args1 <- commandArgs(TRUE)

library(parallel)

mc1 <- 6#25

s <- as.numeric(args1[1])
coal <- args1[2]
if (coal=="beta"){
#' Inference parameter Beta
coal_inf <- c("min"=1,"max"=2,"steps"=20)
r_inf <- c("min"=0,"max"=25,"steps"=50)
misiden_inf <- c("min"=0.001,"max"=0.201,"steps"=20)
}
#' Inference parameter Phi

if (coal=="psi"){
coal_inf <-  c("min"=0,"max"=1,"steps"=20)
r_inf <- c("min"=0,"max"=25,"steps"=50)
misiden_inf <- c("min"=0.001,"max"=0.201,"steps"=20)
}

#' Phi table 

if (coal=="beta"){
  phi.table <- read.delim(paste0("phi_tables/beta.phi.table_n",
                                 s), header=FALSE)
}
if (coal=="psi"){
  phi.table <- read.delim(paste0("phi_tables/psi.phi.table_n",
                                 s), header=FALSE)
}

phi.table <- phi.table[,-1]

#' Load simulated data to then infer the parameters of
load(paste0("simulated_snpdata/sims_errortest_betapsi_s",s,"_ongrid.RData"))


if (coal=="beta"){
  sneq_c <- sneq_b; sims1 <- sims_b
}
if (args1[2]=="psi"){
  sneq_c <- sneq_d; sims1 <- sims_d
}


## Estimate misiden x via sneq 
misiden_phylofun_pointest <- function(sfs1,sneq1){sneq1/(2*sum(sfs1))}

misiden_pointest <- vector("list",length(sneq_c))
for (i in 1:length(sneq_c)){
  misiden_pointest[[i]] <- apply(cbind(sims1[[i]],sneq_c[[i]]),1,
                                 function(v){misiden_phylofun_pointest(v[-length(v)],v[length(v)])})  
}



params_phi <- expand.grid(r=seq(r_inf[1],r_inf[2],length.out = r_inf[3]+1),
                           beta=seq(coal_inf[1],coal_inf[2],length.out = coal_inf[3]+1))

##Function to estimate the pseudolikelihood given a misiden prob, omitting
##the additional term for the pseudolikelihood

pseu_logl <- function(sfs1,x,compl=FALSE,S_neq=0){
  logl_noconst <- rep(-1,nrow(phi.table))
  x_max <- rep(-1,nrow(phi.table))
  S <- sum(sfs1)
  for (j in 1:nrow(phi.table)){
    phi <-  unname(unlist(phi.table[j,]))
    temp1 <- 0
    for (l in seq(along=phi)){ 
      temp1 <-  temp1 + sfs1[l]*log(phi[l]*(1-x)+phi[length(phi)+1-l]*x)}
    if (compl){
      if (min(x)>0){
        binp_e1 <- 2*x/(1+2*x)
        binn_e1 <- S_neq + sum(sfs1)
        temp1 <- temp1 + lchoose(binn_e1,S_neq)+S_neq*log(binp_e1)+sum(sfs1)*log(1-binp_e1)
      }}
    if (!all(is.na(temp1))){
      x_max[j] <- x[which.max(temp1)]} else {
        x_max[j] <- NA
      }
    logl_noconst[j] <- max(temp1)}
  ml_pos <- which.max(logl_noconst)
  ml_val <- logl_noconst[ml_pos] - S + S*log(S) - sum(lfactorial(sfs1)) 
  ml_params <- unlist(params_phi[ml_pos,])
  ml_x <- x_max[ml_pos]
  return(c(mlpos=ml_pos,ml_params,
           est_misiden=ml_x,logPsL=ml_val))
}


infparams_pe <- function(test_set){
print(paste("pe:",test_set))  
data1 <- cbind(sims1[[test_set]],misiden_pointest[[test_set]])
est_R <- sapply(1:nrow(data1),function(i){pseu_logl(data1[i,][-ncol(data1)],
                                                     data1[i,][ncol(data1)])})
return(est_R)
}

infparams_cpslk <- function(test_set){
  print(paste("complik:",test_set))  
  est_R2 <- sapply(1:nrow(sims1[[test_set]]),function(i){pseu_logl(sims1[[test_set]][i,],
                                             x = seq(misiden_inf[1],misiden_inf[2],
                                             length.out = misiden_inf[3]+1),compl = TRUE,
                                             S_neq = sneq_c[[test_set]][i])})
  return(est_R2)}



est_R2_list <- mclapply(1:length(sims1),infparams_cpslk,mc.cores = mc1)


save(est_R2_list,file=paste0("estimR2err_",coal,"_s",s,"_ongrid.RData"))
