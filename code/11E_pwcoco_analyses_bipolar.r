library(dplyr)

#Using the TopSNP file to loop over the datasets

setwd("/dir/revision_03_07_25/pwcoco/trait/data")
data1<- read.table("topregions_asthma_bipolar.txt", sep="\t", header=T)
data2<- read.table("topregions_bip_asthma.txt", sep="\t", header=T)

data<- rbind(data1,data2)

ID<- data %>% split(x=data, f=data$SNP)

lapply(list(ID), function(x){
  
  for(x in ID) {
  
# Chromosome information to build the command 
chr= x$chr.exposure
print(chr)

# MAF to build the command
maf= 0.01

#Using SNP name as pattern to aid extraction from the two directories
pattern=x$SNP
pattern=paste0("*",pattern,"*")
print(pattern)

#List all files for a given snp
A_trait <- list.files(path = "/dir/revision_03_07_25/pwcoco/trait/data/bipolar/trait1", pattern=pattern, recursive = TRUE)
print(A_trait)

B_trait<- list.files(path = "/dir/revision_03_07_25/pwcoco/trait/data/bipolar/trait2", pattern=pattern, recursive = TRUE)
print(B_trait)

#Assign the full directory so that pwcoco can work 
A_trait_dir<- paste0("/dir/revision_03_07_25/pwcoco/trait/data/bipolar/trait1/", A_trait)
B_trait_dir<- paste0("/dir/revision_03_07_25/pwcoco/trait/data/bipolar/trait2/", B_trait)

setwd("/dir/pwcoco_update/pwcoco/build")

#Build the commands
bim_cmd<- paste("--bfile", paste0("/dir/REF_PANEL/refpanel_alspac/CHR", chr))
chr_cmd<- paste("--chr", chr)
maf_cmd<- paste("--maf", maf)

#In order for each file to be used sequentially and matched, building a loop

lapply(list(A_trait_dir), function(r){
  
for(r in A_trait_dir) {

sumstats1<- paste("--sum_stats1", r)
print(sumstats1)

lapply(list(B_trait_dir), function(y){
  
  for(y in B_trait_dir) {

sumstats2<- paste("--sum_stats2", y)
print(sumstats2)

out_cmd<- paste("--out", paste0("/dir/revision_03_07_25/pwcoco/trait/results/asthma_bipolar"))
print(out_cmd)

cmd <- paste("./pwcoco", bim_cmd,
				sumstats1,  
				sumstats2,
				out_cmd,
				chr_cmd, maf_cmd)
print(cmd)
system(cmd)

}})
}})
}})
