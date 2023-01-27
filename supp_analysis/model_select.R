#' Give kappa and e_ts and the ith row of phi_tab (+ corresponding alpha and g)

f_pseul_noconst <- function(i,e_ts,kap,phi_tab){sc_sfs <- as.vector(unname(phi_tab[i,-1]))
l_sfs <- length(sc_sfs)
e_tv <- kap*e_ts/(1+(kap-1)*e_ts)
#fact2 <- (1+kap)*e_ts
#fact2 <- fact2/(1+2*kap*e_ts)
out1 <- sum(as.integer(SFS$ts)*log((1-e_ts)*sc_sfs+e_ts*sc_sfs[l_sfs:1]))
if (e_ts>0){out1 <- out1 + Sneq$ts*log(2*kap*e_ts/(1+2*kap*e_ts))
out1 <- out1 - sum(SFS$ts)*log(1+2*kap*e_ts)}
out1 <- out1 + sum(as.integer(SFS$tv)*log((1-e_tv)*sc_sfs+e_tv*sc_sfs[l_sfs:1]))
if (e_ts>0){out1 <- out1 + Sneq$tv*log((1+kap)*e_ts/(1+2*kap*e_ts))
out1 <- out1 + sum(SFS$tv)*log((1+(kap-1)*e_ts)/(1+2*kap*e_ts))}
#Add the two additional combinatorial factors for the misid factor
if (e_ts>0){out1 <- out1 + lchoose(sum(SFS$ts)+Sneq$ts,Sneq$ts)
out1 <-out1 + lchoose(sum(SFS$tv)+Sneq$tv,Sneq$tv)}
return(out1)}



#' Read in phi table and define corresponding parameters 
#phi_names <- 
phi_tab_alpha <- read.delim(paste0("~/forschung/multiplemerger-project/",
                                   "data_analysis_current/phi_tables_beta/phi.table.",
                                   sample_sfs), 
                            header=FALSE)
phi_tab_alpha <- as.matrix(phi_tab_alpha)

n_alphas <- 20
alpha_steps <- 1 + (0:n_alphas)/n_alphas   
n_rhos <- 50
rho_steps <- (0:n_rhos)/n_rhos*25


misi <- seq(0,.15,length.out=16)
if (!allow_ezero){misi[1] <- 0.00001}
kappa_est <- (sum(SFS$tv)+Sneq$tv)/(2*(sum(SFS$ts)+Sneq$ts))

kappa <- c(1,2/3*kappa_est,kappa_est,3/2*kappa_est) # use kappa=1 for non-ts/tv version and 
# a naive kappa estimate
# given as  "# tv/# ts" (debatable)  

beta_para <- expand.grid(g=rho_steps,alpha=alpha_steps,
                         e_ts=misi,kappa=kappa)
resb <- cbind(phirow=rep(1:nrow(phi_tab_alpha),length(misi)*length(kappa)),beta_para)  




#' TODO find max
res_v <- future_sapply(1:nrow(resb),function(j){f_pseul_noconst(resb$phirow[j],
                                                         e_ts = resb$e_ts[j],
                                                         kap = resb$kappa[j],
                                                         phi_tab = phi_tab_alpha)})


#' Psi
#' Read in phi table and define corresponding parameters 

phi_tab_psi <- read.delim(paste0("~/forschung/multiplemerger-project/",
                                 "data_analysis_current/phi_tables_smallpsi/phi.table.",
                                 sample_sfs), 
                          header=FALSE)
phi_tab_psi <- as.matrix(phi_tab_psi)

#n_alphas <- 20
#alpha_steps <- 1 + (0:n_alphas)/n_alphas   
#n_rhos <- 50
#rho_steps <- (0:n_rhos)/n_rhos*25
#misi <- seq(0,.15,length.out=16)

n_psis <- 20 
psi_steps <- 0.2*(0:n_psis)/n_psis   
n_rhos <- 60
rho_steps <- (0:n_rhos)/n_rhos*30
#rhos <- c(0,30,60)
misi <- seq(0,.2,length.out=21)

if (!allow_ezero){misi[1] <- 0.00001}

kappa_est <- (sum(SFS$tv)+Sneq$tv)/(2*(sum(SFS$ts)+Sneq$ts))
kappa <- c(1,2/3*kappa_est,kappa_est,3/2*kappa_est) # use kappa=1 for non-ts/tv version and 
# a naive kappa estimate
# given as  "# tv/# ts" (debatable)  

psi_para <- expand.grid(g=rho_steps,psi=psi_steps,
                        e_ts=misi,kappa=kappa)
resp <- cbind(phirow=rep(1:nrow(phi_tab_psi),length(misi)*length(kappa)),psi_para)  




#' TODO find max
res_w <- future_sapply(1:nrow(resp),function(j){f_pseul_noconst(resp$phirow[j],
                                                         e_ts = resp$e_ts[j],
                                                         kap = resp$kappa[j],
                                                         phi_tab = phi_tab_psi)})

cat(sample_sfs," beta\n")
sameSFS <- all(SFS$ts+SFS$tv==sfs_as_paper$V2)
cat("Same SFS:",sameSFS,"\n")
if (!sameSFS){
  summary((SFS$ts+SFS$tv-sfs_as_paper$V2)/sfs_as_paper$V2)}

res_list <- list(sample_sfs,beta_al_k=NA,beta_k1=NA,beta_khat=NA,
                 psi_allk=NA,psi_k1=NA,psi_khat=NA,
                 modsel=NA,modselk1=NA,modselkhat=NA,
                 sumSFS_ts=sum(SFS$ts),sumSFS_tv=sum(SFS$tv),
                 sneqsum_ts=Sneq$ts,sneqsum_tv=Sneq$tv,
                 sameSFS=sameSFS)

whichb_kappaest <- sapply(resb$kappa,function(x){isTRUE(all.equal(x,kappa_est))})

(res_list$beta_al_k <- resb[res_v==max(res_v),])

(res_list$beta_k1 <- resb[resb$kappa==1,][which.max(res_v[resb$kappa==1]),])

res_list$beta_khat <- resb[whichb_kappaest,][which.max(res_v[whichb_kappaest]),]

top5_b <- order(res_v,decreasing = TRUE)[1:5]
resb[top5_b,]

cat(sample_sfs," psi\n")
whichp_kappaest <- sapply(resp$kappa,function(x){isTRUE(all.equal(x,kappa_est))})


(res_list$psi_allk <- resp[res_w==max(res_w),])

(res_list$psi_k1 <- resp[resp$kappa==1,][which.max(res_w[resp$kappa==1]),])

res_list$psi_khat <- resp[whichp_kappaest,][which.max(res_w[whichp_kappaest]),]

top5_p <- order(res_w,decreasing = TRUE)[1:5]
resp[top5_p,]

#Optimize for all kappa values
logml_beta <- max(res_v[resb$alpha!=2])
logml_psi <- max(res_w[resp$psi>0])
logml_km <- max(res_v[resb$alpha==2])

if (max(logml_beta,logml_psi)-logml_km > log(10)){res_list$modsel <- "mmc"
if (logml_beta-logml_psi > log(10)){res_list$modsel <- "beta"}
if (logml_beta-logml_psi < -log(10)){res_list$modsel <- "psi"}} else {
  res_list$modsel <- "km"
}

#Optimize for kappa=1
logml_beta <- max(res_v[resb$alpha!=2 & resb$kappa==1])
logml_psi <- max(res_w[resp$psi>0 & resp$kappa==1])
logml_km <- max(res_v[resb$alpha==2 & resb$kappa==1])

if (max(logml_beta,logml_psi)-logml_km>log(10)){res_list$modselk1 <- "mmc"
if (logml_beta-logml_psi > log(10)){res_list$modselk1 <- "beta"}
if (logml_beta-logml_psi < -log(10)){res_list$modselk1 <- "psi"}} else {
  res_list$modselk1 <- "km"
}

#Optimize for kappa=\hat{kappa}
logml_beta <- max(res_v[resb$alpha!=2 & whichb_kappaest])
logml_psi <- max(res_w[resp$psi>0 & whichp_kappaest])
logml_km <- max(res_v[resb$alpha==2 & whichb_kappaest])

if (max(logml_beta,logml_psi)-logml_km>log(10)){res_list$modselkhat <- "mmc"
if (logml_beta-logml_psi > log(10)){res_list$modselkhat <- "beta"}
if (logml_beta-logml_psi < -log(10)){res_list$modselkhat <- "psi"}} else {
  res_list$modselkhat <- "km"
}


save(res_list,file=paste0("results/",sample_sfs,
                          ifelse(allow_ezero,"","_no0"),".RData"))

