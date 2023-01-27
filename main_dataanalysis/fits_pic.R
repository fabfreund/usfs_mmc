#' Plot the fits of the best-fitting 
#' models from KM, BETA, PSI. 
#' as well as observed SFS and GC-bias corrected SFS
#' Plots are ordered by fit measured by 
#' Cohen's d

#devanil transforms the SFS by multiplying its ith entry with i
devanil <- TRUE#FALSE

end1 <- ifelse(devanil,"tf","")


devanilla <- function(sfs1){
  seq(along=sfs1)*sfs1}

load("res_summaries/mlres_BF_effectsize.RData")

#function to change filenames to recognie GC-bias corrected SFS
rename <- function(name1){paste0("../sfs_data/gcbias_corrected_sfs/",strsplit(name1,split = ".txt"),
                                 "_neutral.txt")}

pdf(paste0("pics/fits",end1,".pdf"))
mar.default <- par()$mar
par(mar = mar.default + c(0, 1, 0, 0)) 
k <- 0
for (name1 in ml_bf_eff$file){
temp1 <- read.table(rename(name1))
neutralsfs <- temp1[,2]/sum(temp1[,2])
if (devanil){neutralsfs <- devanilla(neutralsfs)}
k <- k+1 
#if (k==9){next}
xvals <- seq(along=sfs_list[[name1]][["obs"]])
if (devanil){
for (mod1 in c("obs","km","beta","psi")){
sfs_list[[name1]][[mod1]] <- devanilla(
  sfs_list[[name1]][[mod1]])}}
range1 <- range(unlist(sfs_list[[name1]]))
range1 <- range(range1,neutralsfs)

if (!devanil){
plot(xvals,sfs_list[[name1]][["obs"]],
     ylim=range1,pch=20,
       ylab=expression(E(xi[i])/sum(E(xi))),
       xlab="SFS classes",
     main=ml_bf_eff$species[k]
     )} else {
       plot(xvals,sfs_list[[name1]][["obs"]],
            ylim=range1,pch=20,
            ylab=expression(i*E(xi[i])/sum(E(xi))),
            xlab="SFS classes",
            main=ml_bf_eff$species[k]) 
     }
ltys <- switch(ml_bf_eff$selected_model[k],
               "BETA"=c(2,1,2),
               "PSI"=c(2,2,1),
               "KM"=c(1,2,2),
               "MMC"=c(2,1,1))
effsize <- switch(ml_bf_eff$selected_model[k],
                  "BETA"=ml_bf_eff[k,8],
                  "PSI"=ml_bf_eff[k,15],
                  "KM"=ml_bf_eff[k,22],
                  "MMC"=min(ml_bf_eff[k,8],
                            ml_bf_eff[k,15]))
effsize <- round(effsize,2)
points(xvals,sfs_list[[name1]][["km"]],
       pch=20,col="black",
       type="l",lty=ltys[1])
points(xvals,sfs_list[[name1]][["beta"]],
       pch=20,col="blue",
       type="l",lty=ltys[2])
points(xvals,sfs_list[[name1]][["psi"]],
       pch=20,col="red",
       type="l",lty=ltys[3])
points(xvals,neutralsfs,
       pch=4,col="gray",
       type="p")
legend("topright",legend=c("Observed","GC-bias corrected",
paste0("ML KM: (",paste0(unname(ml_bf_eff[k,18:19]),collapse = ","),")"),
paste0("ML Beta: (",paste0(unname(ml_bf_eff[k,3:5]),collapse = ","),")"),
paste0("ML PSI: (",paste0(unname(ml_bf_eff[k,10:12]),collapse = ","),")")),
pch=c(20,4,NA,NA,NA),
col=c("black","gray","black","blue","red"),
lty=c(NA,NA,ltys))
text(max(xvals)/4,range1[2],paste("Cramer's d:",effsize))         
  }
dev.off()

