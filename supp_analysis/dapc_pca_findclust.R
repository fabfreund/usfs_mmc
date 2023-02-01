library(adegenet)
library(pegas)



#Read in fasta as genind
path1 <- "../example_fasta/"
trafo_fas <- function(str1){out1 <- read.dna(file=str1,format="fasta")
                            return(DNAbin2genind(out1))}
file1 <- "Aptenodytes_concat_alignment_clean.fas"


#' Read in .fas file - takes a while...
  cat(file1,"\n")
  fas_rdata <- trafo_fas(paste0(path1,file1))


#remove outgroup
outgroup_pos <- c(7:8)
fas_rdata <- fas_rdata[-outgroup_pos,drop=TRUE]




#' Double-centered PCA
  cat("DC-PCA \n")
  temp1 <- scaleGen(fas_rdata,NA.method= "mean",center=FALSE,scale=FALSE)
  #Do DC-PCA
  temp1 <- temp1 - matrix(rowMeans(temp1),byrow=FALSE,nrow=nrow(temp1),
                          ncol=ncol(temp1)) - matrix(colMeans(temp1),byrow = TRUE,
                                                     nrow=nrow(temp1),ncol=ncol(temp1)) +
          mean(temp1)
  pca_res <- dudi.pca(temp1,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)

#' Extract number of eigenvectors, we retain all for dapc
eig1 <- length(pca_res$eig)

dapc_res <- find.clusters(fas_rdata,n.pca = eig1,
                                   #choose.n.clust = FALSE,
                                   #criterion = "goodfit",#"diffNgroup",
                                   max.n.clust = as.integer(min(eig1,17,nrow(fas_rdata@tab))))
cat("chosen cluster number:", names(dapc_res$stat),"\n")


pdf("aptenodytes_res/pca_dapc_aptenodytes.pdf")
  col1 <- dapc_res$grp
  col0 <- c("red","black",rainbow(18))
  levels(col1) <- col0[1:length(levels(col1))]
  col1 <- as.character(col1)
  plot(pca_res$li,col=col1,cex=2,pch=3,main="Aptenodytes patagonicus")
dev.off()

pdf("aptenodytes_res/pca_dapc_aptenodytes_2axis.pdf")
  col1 <- dapc_res$grp
  col0 <- c("red","black",rainbow(18))
  levels(col1) <- col0[1:length(levels(col1))]
  col1 <- as.character(col1)
  range1 <- range(pca_res$li[,c(1,2)])
  eig_v <- pca_res$eig
  eig_v <- eig_v/sum(eig_v)
  plot(pca_res$li[,c(1,2)],col=col1,cex=2,pch=3,main="Aptenodytes patagonicus",
       xlab=paste0("PCA 1, frac Var = ",round(eig_v[1],2)),
       ylab=paste0("PCA 2, frac Var = ",round(eig_v[2],2)),
       xlim=range1,ylim=range1)
dev.off()


#' Record PCA results needed for the plots

pca_coord <- pca_res$li
eig_val <- pca_res$eig


save(dapc_res,pca_coord,eig_val,file="aptenodytes_res/res_dapcpca_aptenodytes.RData",
     compress = "xz")

