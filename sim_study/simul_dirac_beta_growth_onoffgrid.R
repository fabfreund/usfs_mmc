#' IMPORTANT: the following script includes the absolute path of the 
#' simulation script and has to be checked/adjusted before running
#' There are two versions - one is commented out. We marked the relevant lines 
#' with the text "2 options - ongrid and offgrid"
 
#' 2 options - ongrid and offgrid
#source("param_set_ongrid.R")
source("param_set_offgrid.R")

#' Set sample size
s <- 25

#Compute Watterson estimator for each class
setwd(home_dir)
theta_watt_est <- function(s1,n_mut1,r,cp,model="B"){
  if (model=="B"){  
    com1 <- paste(com_phi_b,"-sampleSize",s1,"-MinAlpha",cp,"-MaxAlpha",cp,"-NoStepsAlpha 0  -MinRho",r,"-MaxRho",r,"-NoStepsRho 0 -OUT",out1)}
  if (model=="D"){
    com1 <- paste(com_phi_d,"-sampleSize",s1,"-MinPsi",cp,"-MaxPsi",cp,"-NoStepsPsi 0  -MinRho",r,"-MaxRho",r,"-NoStepsRho 0 -OUT",out1) 
  }  
  system(com1)
  system(paste0("mv ",out1_short,"* bl.temp"))
  bl <- read.delim("bl.temp", header=FALSE)[1,1]
  system("rm bl.temp")
  return(unname(2*n_mut1/bl))
}

theta_b <- apply(params_beta_exp,1,function(v){theta_watt_est(s,n_mut,v[1],v[2],"B")})
theta_d <- apply(params_dirac_exp,1,function(v){theta_watt_est(s,n_mut,v[1],v[2],"D")})


#Simulation code


test_s1 <- s
test_theta1 <- theta_b[2]
test_cp1 <- params_beta_exp[2,]

sfs_sim <- function(s1,theta1,cp1=c(0,2),model="B",nsim1=25,seed1=1){
  if ((model=="B") & (cp1[2] != 2)){theta1 <- theta1/(2*10000^(cp1[2]))} #to adjust timescale, now is ^alpha
  if ((model=="B") & (cp1[2]==2)) {theta1 <- theta1/20000}
  if (model=="D") {theta1 <- (cp1[2]^2*theta1)/(2*100)}
  if (model=="B"){  
    com1 <- paste(com_sim_b,"-s",s1,"-a",cp1[2],"-r",cp1[1],"-mu",theta1,"-num",nsim1,
                  "-SFS -out",paste0(out2,seed1,"__"))}
  if (model=="D"){
    com1 <- paste(com_sim_d,"-s",s1,"-psi",cp1[2],"-r",cp1[1],"-mu",theta1,"-num",nsim1,
                  "-SFS -out",paste0(out2,seed1,"__")) 
  }  
  system(com1)
  system(paste0("mv ",out2_short,seed1,"__","* sfs.temp",seed1))
  sfs <- read.delim(paste0("sfs.temp",seed1), header=FALSE,
                    comment.char = "#")[,-(1:3)]
  system(paste0("rm sfs.temp",seed1))
  return(sfs)}  

n_loci_sfs <- function(sfs1){apply(sfs1,2,sum)}

sim_sfs_param <- function(s1,theta1,cp1=c(0,2),model="B",
                          nsim1=25,nloc=100,seed2=1){
  replicate(n=nsim1,expr = n_loci_sfs(sfs_sim(s1,theta1,cp1,model,nloc,seed2)))}

#Simulate Beta+growth, including non-polarised site sim
sims_b_t <- mclapply(1:nrow(params_beta_exp),
                     function(i){t(sim_sfs_param(s,theta_b[i],params_beta_exp[i,],"B",n_simul,n_loci,seed2=i))},
                     mc.cores = mc1)

sims_b <- sims_b_t #Used if misiden=0 allowed, NULL else
sneq_b <- vector("list",length(sims_b))
for (i in seq(along=sims_b)){
sneq_b[[i]] <- rep(0,n_simul)}


mis_SFS <- function(sfs1,x){
  corr_counts <- sapply(sfs1,function(n){rbinom(1,n,1-x)})
  n1 <- length(sfs1)+1 #original sample size
  sfsmis <- sapply(seq(along=sfs1),
                   function(ni){corr_counts[ni]+sfs1[n1-ni]-corr_counts[n1-ni]}) 
  
  return(sfsmis)
}

mis_SFS_matrix <- function(m,x){t(apply(m,1,mis_SFS,x))}

for (i in 2:length(misiden)){#if misiden=0 allowed, start at 2
  sims_b_t2 <- lapply(sims_b_t,mis_SFS_matrix,misiden[i])
  sims_b <- c(sims_b,sims_b_t2) 
  sneq2 <- lapply(sims_b_t,function(m){xtemp <- 2*misiden[i]/(1+2*misiden[i])
  Sobs <- apply(m,1,sum)
  S0hat <- floor((1+2*misiden[i])*Sobs)
  hf1 <- function(n1){rbinom(1,n1,xtemp)}
  tempsneq <- sapply(S0hat,hf1)
  return(tempsneq)
  })
  sneq_b <- c(sneq_b,sneq2)
}

#Simulate Dirac, including non-polarised base count

sims_d_t <- mclapply(1:nrow(params_dirac_exp),
                     function(i){t(sim_sfs_param(s,theta_d[i],                              params_dirac_exp[i,],"D",n_simul,n_loci,seed2=i))},
                     mc.cores = mc1)

sims_d <- sims_d_t#sims_b_t #Used if misiden=0 allowed
sneq_d <- vector("list",length(sims_d))
for (i in seq(along=sims_d)){
  sneq_d[[i]] <- rep(0,n_simul)}


mis_SFS <- function(sfs1,x){
  corr_counts <- sapply(sfs1,function(n){rbinom(1,n,1-x)})
  n1 <- length(sfs1)+1 #original sample size
  sfsmis <- sapply(seq(along=sfs1),
                   function(ni){corr_counts[ni]+sfs1[n1-ni]-corr_counts[n1-ni]}) 
  
  return(sfsmis)
}

mis_SFS_matrix <- function(m,x){t(apply(m,1,mis_SFS,x))}

for (i in 2:length(misiden)){#if misiden=0 allowed, start at 2
  sims_d_t2 <- lapply(sims_d_t,mis_SFS_matrix,misiden[i])
  sims_d <- c(sims_d,sims_d_t2) 
  sneq2_d <- lapply(sims_d_t,function(m){xtemp <- 2*misiden[i]/(1+2*misiden[i])
  Sobs <- apply(m,1,sum)
  S0hat <- floor((1+2*misiden[i])*Sobs)
  hf1 <- function(n1){rbinom(1,n1,xtemp)}
  tempsneq <- sapply(S0hat,hf1)
  return(tempsneq)
  })
  sneq_d <- c(sneq_d,sneq2_d)
}

#' 2 options - ongrid and offgrid
#save(sims_b,sims_d,sneq_b,sneq_d,file=paste0("simulated_snpdata/sims_errortest_betapsi_s",s,"_ongrid.RData"))
save(sims_b,sims_d,sneq_b,sneq_d,file=paste0("simulated_snpdata/sims_errortest_betapsi_s",s,"_offgrid.RData"))