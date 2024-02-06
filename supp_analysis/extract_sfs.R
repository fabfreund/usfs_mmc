#' Libraries
library(adegenet)
library(pegas)


#' Source filter scripts
source("aux_tstv/funct_and_consts.R")

#' Specify SFS names
set_fasta <- "../example_fasta/Aptenodytes_concat_alignment_clean.fas"

#' Read in data
data1 <- read.dna(file=set_fasta,format="fasta")
gc()

#' Extract outgroup 
out_rows <- c(7,8) #which rows/lines in the fasta file are from outgroup?
rownames(data1)[out_rows]
rownames(data1)[-out_rows]

#' first filter step
keep_sit1 <- apply(as.character(data1),2,ccheck_awful)
data2 <- data1[,keep_sit1]
rm(data1)
gc()

#' Count bases - which ones can be oriented etc
res1 <- apply(as.character(data2),2,ccount_bases)
gc()


#' Extract SFS from transitions and SFS from transversions  
SFS <- sfs_as_sfs(table(res1["nmut",res1["iclean",]=="1"]))


gc()




