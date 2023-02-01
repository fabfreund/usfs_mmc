#' Libraries
library(adegenet)
library(pegas)
#library(compiler)
library(future.apply)

#' Source filter scripts
source("aux_tstv/funct_and_consts.R")

#' Specify SFS names
sfs_name <- "Aptenodytes_SFS.txt"
#' Specify fasta files
set_fasta <- "../example_fasta/Aptenodytes_concat_alignment_clean.fas"

#' Read in data
data1 <- read.dna(file=set_fasta,format="fasta")

#' Extract outgroup 
out_rows <- c(7,8) #which rows/lines in the fasta file are from outgroup?
rownames(data1)[out_rows]
rownames(data1)[-out_rows]

#' first filter step
keep_sit1 <- apply(as.character(data1),2,ccheck_awful)
data2 <- data1[,keep_sit1]
rm(data1)

#' Count bases - which ones can be oriented etc
res1 <- apply(as.character(data2),2,ccount_bases)

#' Compare SFS - from this filtering and from the main text filtering


sfs_as_paper <- read.delim(paste0("../sfs_data/real_sfs/",sfs_name), 
                           header=FALSE)


#' Extract SFS from transitions and SFS from transversions  
SFS_ts <- sfs_as_sfs(table(res1["nmut",res1["iclean",]=="1" & res1["ts",]=="1"]))
SFS_tv <- sfs_as_sfs(table(res1["nmut",res1["iclean",]=="1" & res1["ts",]=="0"]))



gc()


SFS <- list(ts=SFS_ts,tv=SFS_tv)

#' Extract non-orientable sites (3rd base in outgroup) - both for in-sample ts and tv
Sneq <- list(ts=sum(res1["ts",]=="1" & res1["pol",]=="0" & res1["iclean",]==1,
                    na.rm = TRUE),
             tv=sum(res1["ts",]=="0" & res1["pol",]=="0" & res1["iclean",]==1
                    ,na.rm = TRUE))


#' Run model selection and ML parameter optimisation
source("aux_tstv/model_select.R")
#' Save output
save(res_list,file=paste0("aptenodytes_res/tvts_",sfs_name,".RData"))