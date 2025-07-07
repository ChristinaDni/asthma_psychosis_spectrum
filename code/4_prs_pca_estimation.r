#PRS PCA
#First import all the estimated PRS

library(dplyr)

setwd("/dir/revision_03_07_25/prs_pca/profiles")

s1<- read.table("asthma_prs_child.S1.txt", sep="\t",header=T)
s1<- s1%>%select(aln,qlet,SCORE)%>%rename(SCORE_1=SCORE)

s2<- read.table("asthma_prs_child.S2.txt", sep="\t",header=T)
s2<- s2%>%select(aln,qlet,SCORE)%>%rename(SCORE_2=SCORE)

s3<- read.table("asthma_prs_child.S3.txt", sep="\t",header=T)
s3<- s3%>%select(aln,qlet,SCORE)%>%rename(SCORE_3=SCORE)

s4<- read.table("asthma_prs_child.S4.txt", sep="\t",header=T)
s4<- s4%>%select(aln,qlet,SCORE)%>%rename(SCORE_4=SCORE)

s5<- read.table("asthma_prs_child.S5.txt", sep="\t",header=T)
s5<- s5%>%select(aln,qlet,SCORE)%>%rename(SCORE_5=SCORE)

s6<- read.table("asthma_prs_child.S6.txt", sep="\t",header=T)
s6<- s6%>%select(aln,qlet,SCORE)%>%rename(SCORE_6=SCORE)

s7<- read.table("asthma_prs_child.S7.txt", sep="\t",header=T)
s7<- s7%>%select(aln,qlet,SCORE)%>%rename(SCORE_7=SCORE)

s8<- read.table("asthma_prs_child.S8.txt", sep="\t",header=T)
s8<- s8%>%select(aln,qlet,SCORE)%>%rename(SCORE_8=SCORE)

s9<- read.table("asthma_prs_child.S9.txt", sep="\t",header=T)
s9<- s9%>%select(aln,qlet,SCORE)%>%rename(SCORE_9=SCORE)

s10<- read.table("asthma_prs_child.S10.txt", sep="\t",header=T)
s10<- s10%>%select(aln,qlet,SCORE)%>%rename(SCORE_10=SCORE)

s11<- read.table("asthma_prs_child.S11.txt", sep="\t",header=T)
s11<- s11%>%select(aln,qlet,SCORE)%>%rename(SCORE_11=SCORE)

s12<- read.table("asthma_prs_child.S12.txt", sep="\t",header=T)
s12<- s12%>%select(aln,qlet,SCORE)%>%rename(SCORE_12=SCORE)

s13<- read.table("asthma_prs_child.S13.txt", sep="\t",header=T)
s13<- s13%>%select(aln,qlet,SCORE)%>%rename(SCORE_13=SCORE)

#Merge them all 

l <- list(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)
all<- purrr::reduce(.x = l, merge, by = c("aln","qlet"), all = T)

prs.pc <- function(dat){
  xo <- scale(as.matrix(dat[,c(-1,-2)]))  ## scale cols of matrix of only PRSs (remove FID, IID, assuming they are first two cols)
  message("\nGenerating PRS-PCA from the following Pt: ", paste0(colnames(xo), sep = ", "),"\n")
  g <- prcomp(xo)   ## perform principal components
  pca.r2 <- g$sdev^2/sum(g$sdev^2)    ## calculate variance explained by each PC
  pc1.loadings <- g$rotation[,1];     ## loadings for PC1
  pc2.loadings <- g$rotation[,2]      ## loadings for PC2
  ## flip direction of PCs to keep direction of association
  ## (sign of loadings for PC1 is arbitrary so we want to keep same direction)
  if (mean(pc1.loadings>0)==0){     
    pc1.loadings <- pc1.loadings*(-1) 
    pc2.loadings <- pc2.loadings*(-1)
  }
  ## calculate PRS-PCA (outputs PC1 and PC2 even though PC1 sufficient)
  pc1 <- xo %*% pc1.loadings
  pc2 <- xo %*% pc2.loadings
  out <- dat[,c(1,2)]
  out[,"prs.pc1"] <- scale(pc1)   ## rescales PRS-PCA1 
  out[,"prs.pc2"] <- scale(pc2)  ## rescales PRS-PCA2
  return(list(data=out,r2=pca.r2,loadings=pc1.loadings))
}

prs_pc<- prs.pc(all)
prs_pc1<- prs_pc$data%>%select(aln,qlet,prs.pc1)

write.table(prs_pc1, "/dir/revision_03_07_25/prs_pca/prs_pc1/prs_pc1.txt",sep="\t",row.names=F,col.names=T, quote=F)
