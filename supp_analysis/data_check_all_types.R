#' Libraries
library(adegenet)
library(pegas)
#library(compiler)
library(future.apply)
#' General source
source("funct_and_consts.R")

#' Specify SFS names
sfs_names <- c("Aptenodytes_SFS.txt","Emys_SFS.txt",
               "Halictus_SFS.txt","Lepus_SFS.txt","Ostrea_SFS.txt","Parus_SFS.txt",
               "Physa_SFS.txt")
#' Specify fasta files
sets_fasta <- list.files("../../../ushape_data/")


for (it1 in 1:7){
data_no <- it1 #1:7

#' Read in data
str1 <- sets_fasta[data_no]
data1 <- read.dna(file=paste0("../../../ushape_data/",str1),format="fasta")
str1


#' Extract outgroup
rownames(data1)
out_all <- list(c(7,8),c(3,4,11,12),c(1,2),c(7,8),c(9,10,13,14),c(19,20),
                c(16,17,19,22))
out_rows <- out_all[[data_no]]
rownames(data1)[out_rows]
rownames(data1)[-out_rows]


keep_sit1 <- apply(as.character(data1),2,ccheck_awful)
data2 <- data1[,keep_sit1]
rm(data1)

res1 <- apply(as.character(data2),2,ccount_bases)

#' Compare SFS


sample_sfs <- sfs_names[data_no]
sfs_as_paper <- read.delim(paste0("../Data_USFS_2020/Data_USFS/",sample_sfs), 
                           header=FALSE)


SFS_ts <- sfs_as_sfs(table(res1["nmut",res1["iclean",]=="1" & res1["ts",]=="1"]))
SFS_tv <- sfs_as_sfs(table(res1["nmut",res1["iclean",]=="1" & res1["ts",]=="0"]))



gc()


SFS <- list(ts=SFS_ts,tv=SFS_tv)
Sneq <- list(ts=sum(res1["ts",]=="1" & res1["pol",]=="0" & res1["iclean",]==1,
                    na.rm = TRUE),
             tv=sum(res1["ts",]=="0" & res1["pol",]=="0" & res1["iclean",]==1
                    ,na.rm = TRUE))



source("model_select.R")

}
