#' Plot the fits of the best-fitting 
#' models from KM, BETA, PSI. 
#' as well as observed SFS and GC-bias corrected SFS
#' Plots are ordered alphabetically

#devanil transforms the SFS by multiplying its ith entry with i
devanil <- FALSE

end1 <- ifelse(devanil,"tf","")


devanilla <- function(sfs1){
  seq(along=sfs1)*sfs1}

load("res_summaries/mlres_BF_effectsize.RData")

#' convenience function to set filenames of GC-bias corrected SFS to exactly
#' those of the original, uncorrected SFS
rename <- function(name1){paste0("../sfs_data/gcbias_corrected_sfs/",strsplit(name1,split = ".txt"),
                                 "_neutral.txt")}

pdf(paste0("pics/fits",end1,".pdf"))
mar.default <- par()$mar
par(mar = mar.default + c(0, 1, 0, 0)) 
for (k in order(ml_bf_eff$species)){
temp1 <- read.table(rename(ml_bf_eff$file[k]))
neutralsfs <- temp1[,2]/sum(temp1[,2])
if (devanil){neutralsfs <- devanilla(neutralsfs)}
xvals <- seq(along=sfs_list[[ml_bf_eff$file[k]]][["obs"]])
if (devanil){
for (mod1 in c("obs","km","beta","psi")){
sfs_list[[ml_bf_eff$file[k]]][[mod1]] <- devanilla(
  sfs_list[[ml_bf_eff$file[k]]][[mod1]])}}
range1 <- range(unlist(sfs_list[[ml_bf_eff$file[k]]]))
range1 <- range(range1,neutralsfs)

if (!devanil){
plot(xvals,sfs_list[[ml_bf_eff$file[k]]][["obs"]],
     ylim=range1,pch=20,
       ylab=expression(E(xi[i])/sum(E(xi))),
       xlab="SFS classes",
     main=ml_bf_eff$species[k]
     )} else {
       plot(xvals,sfs_list[[ml_bf_eff$file[k]]][["obs"]],
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
points(xvals,sfs_list[[ml_bf_eff$file[k]]][["km"]],
       pch=20,col="black",
       type="l",lty=ltys[1])
points(xvals,sfs_list[[ml_bf_eff$file[k]]][["beta"]],
       pch=20,col="blue",
       type="l",lty=ltys[2])
points(xvals,sfs_list[[ml_bf_eff$file[k]]][["psi"]],
       pch=20,col="red",
       type="l",lty=ltys[3])
points(xvals,neutralsfs,
       pch=4,col="darkgray",
       type="p")
legend("topright",legend=c("Observed","GC-bias corrected",
paste0("ML KM: (",paste0(unname(ml_bf_eff[k,18:19]),collapse = ","),")"),
paste0("ML Beta: (",paste0(unname(ml_bf_eff[k,3:5]),collapse = ","),")"),
paste0("ML PSI: (",paste0(unname(ml_bf_eff[k,10:12]),collapse = ","),")")),
pch=c(20,4,NA,NA,NA),
col=c("black","darkgray","black","blue","red"),
lty=c(NA,NA,ltys))
text(max(xvals)/4,range1[2],paste("Cramer's d:",effsize))         
  }
dev.off()