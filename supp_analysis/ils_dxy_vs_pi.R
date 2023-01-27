library(PopGenome)
library(sdprisk)

#' Read in fasta alignment of ingroup and outgroup
fasta_file <- "../example_fasta/Aptenodytes_concat_alignment_clean.fas"
spec_name <- "Aptenodytes"
#' specify which rows contain the outgroup sequences
outgroups <- c(7,8)

data1 <- read.big.fasta(filename = fasta_file,big.data = TRUE,FAST = TRUE,window = 50000)


ind_names <- data1@region.data@populations2[[c(1,1)]]
populations <- list(ing=ind_names[-outgroups],
                    outg=ind_names[outgroups])

data1 <- set.populations(data1,new.populations = populations)


#' Compute
#' nucleotide diversity in sample and in outgroup
#' diversity between (dxy) sample and outgroup
#' the estimate for the divergence time between ingroup and outgroup from Section A.5
#' the probability that the MRCA from the sample is bigger than that estimate 
res1 <- data.frame(pi_insample=-1,pi_ogroup=-1,
                   dxy_in_out=-1,est_tdiv=-1,
                   p_morethantdiv=-1)  

gc()

data1 <- diversity.stats.between(data1)
data1 <- neutrality.stats(data1)
res1$pi_insample <- round(sum(data1@theta_Tajima[,1],na.rm = TRUE)/sum(data1@n.sites),8)
res1$pi_ogroup <- round(sum(data1@theta_Tajima[,2],na.rm = TRUE)/sum(data1@n.sites),8)
res1$dxy_in_out <- sum(data1@nuc.diversity.between[,1],na.rm = TRUE)/sum(data1@n.sites)

#' Compute test
tdiv <-  sum(data1@nuc.diversity.between[,1],na.rm = TRUE)/sum(data1@theta_Tajima[,1],na.rm = TRUE) -1
n <- max(sapply(populations,length))
rates <- sapply(n:2,function(i){choose(i,2)})
prob_less_tdiv <- phypoexp(q = tdiv,rate = rates)
rm(data1)
res1$est_tdiv <- round(tdiv,3)
res1$p_morethantdiv <-round(1-prob_less_tdiv,4)




#save(res1,file = "div_inout_1.RData")
write.table(res1,row.names = FALSE,
            file = paste0("div_inout_",spec_name,".txt"))
