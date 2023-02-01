library(adegenet)
library(pegas)


#Read in fasta as DNAbin
  path1 <- "../example_fasta/Aptenodytes_concat_alignment_clean.fas"
trafo_fas <- function(str1){read.dna(file=str1,format="fasta")}

#' Read in .fas files
fas_rdata <- trafo_fas(path1)

#remove outgroup

outgroup_pos <- 7:8
fas_rdata <- fas_rdata[-outgroup_pos,]


#Compute pi per site for all list entries in non-overlapping 50kb windows
wsize <- 15000

pi_wind <- sw(fas_rdata,width = wsize,step = wsize,FUN=nuc.div,
                   rowAverage = TRUE)


save(pi_wind,file = "aptenodytes_res/pi_windowed_aptenodytes.RData")
