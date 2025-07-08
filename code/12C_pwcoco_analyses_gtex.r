setwd("/dir/revision_03_07_25/pwcoco/trait/results")
data<- read.table("passing_h4.txt", sep="\t", header=T)

#Create the path
data$path1<- paste0("/dir/revision_03_07_25/pwcoco/trait/data/trait1/", data$Dataset1)
data$path2<- paste0("/dir/revision_03_07_25/pwcoco/trait/data/trait2/", data$Dataset2)

#Extract the columns with the filename and merge so that everything can be looped

library(dplyr)
dataset1<- data%>%select(path1, SNP, CHR)%>%
				rename(path=path1)
dataset2<- data%>%select(path2, SNP, CHR)%>%
				rename(path=path2)

data<- rbind(dataset1,dataset2)

ID<- data %>% split(x=data, f=data$path)

lapply(list(ID), function(x){
  
  for(x in ID) {
  
# Chromosome information to build the command 
chr= x$CHR
print(chr)

# MAF to build the command
maf= 0.01

#Using SNP name as pattern to aid extraction from the two directories
pattern=x$SNP
pattern=paste0("*",pattern,"*")
print(pattern)

#List all files for a given snp
mol_trait <- list.files(path = "/dir/revision_03_07_25/pwcoco/tissue/data/gtex/clean", pattern=pattern, recursive = TRUE)
print(mol_trait)

#Assign the full directory so that pwcoco can work 
mol_trait_dir<- paste0("/dir/revision_03_07_25/pwcoco/tissue/data/gtex/clean/", mol_trait)

setwd("/dir/pwcoco_update/pwcoco/build")

#Build the commands
bim_cmd<- paste("--bfile", paste0("/dir/REF_PANEL/refpanel_alspac/CHR", chr))
chr_cmd<- paste("--chr", chr)
maf_cmd<- paste("--maf", maf)

#In order for each file to be used sequentially and matched, building a loop

lapply(list(mol_trait_dir), function(r){
  
for(r in mol_trait_dir) {

sumstats1<- paste("--sum_stats1", r)
print(sumstats1)

lapply(list(x$path), function(y){
  
  for(y in x$path) {

sumstats2<- paste("--sum_stats2", y)
print(sumstats2)

out_cmd<- paste("--out", paste0("/dir/revision_03_07_25/pwcoco/tissue/results/gtex"))
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
