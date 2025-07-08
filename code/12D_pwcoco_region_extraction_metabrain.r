library(dplyr)
library(readr)

setwd("/dir/revision_03_07_25/pwcoco/trait/results/")
data<- read.table("passing_h4.txt", sep="\t", header=T)

ID<- data %>% split(x=data, f=data$SNP)

#Start with cortex

lapply(list(ID), function(x){
  
  for(x in ID) {

#Define chromosome, position and window
#These are necessary for extraction and file naming

snp<- x$SNP
  
chr<- x$CHR
print(chr)

pos<- x$BP
print(pos)

window<- 500000

#Read in the GWAS data to extract the regions
gwas_data<- read_table(paste0("/dir/data/metabrain/clean_cd/b37/B37:2021-07-23-cortex-EUR-80PCs-chr", chr, ".txt"))

gwas_data<- as.data.frame(gwas_data)

print(nrow(gwas_data))

#I want to extract for all available genes in the GWAS

region<- gwas_data[which(gwas_data$BP>=max(pos-window,0) & gwas_data$BP<=pos+window),] 

#Create a name for the file

setwd("/dir/revision_03_07_25/pwcoco/tissue/data/metabrain/raw")

file_out_snp<- paste0("cortex_snp","_", snp, ".txt")
print(file_out_snp)

write.table(region, file_out_snp, sep="\t", quote=F, row.names=F, col.names=T)

}})


#Move to cerebellum

lapply(list(ID), function(x){
  
  for(x in ID) {

#Define chromosome, position and window
#These are necessary for extraction and file naming

snp<- x$SNP
  
chr<- x$CHR
print(chr)

pos<- x$BP
print(pos)

window<- 500000

#Read in the GWAS data to extract the regions
gwas_data<- read_table(paste0("/dir/data/metabrain/clean_cd/b37/B37:2021-07-23-cerebellum-EUR-60PCs-chr", chr, ".txt"))

gwas_data<- as.data.frame(gwas_data)

print(nrow(gwas_data))

#I want to extract for all available genes in the GWAS

region<- gwas_data[which(gwas_data$BP>=max(pos-window,0) & gwas_data$BP<=pos+window),] 

#Create a name for the file

setwd("/dir/revision_03_07_25/pwcoco/tissue/data/metabrain/raw")

file_out_snp<- paste0("cerebellum_snp","_", snp, ".txt")
print(file_out_snp)

write.table(region, file_out_snp, sep="\t", quote=F, row.names=F, col.names=T)

}})

#Hippocampus

lapply(list(ID), function(x){
  
  for(x in ID) {

#Define chromosome, position and window
#These are necessary for extraction and file naming

snp<- x$SNP
  
chr<- x$CHR
print(chr)

pos<- x$BP
print(pos)

window<- 500000

#Read in the GWAS data to extract the regions
gwas_data<- read_table(paste0("/dir/data/metabrain/clean_cd/b37/B37:2021-07-23-hippocampus-EUR-30PCs-chr", chr, ".txt"))

gwas_data<- as.data.frame(gwas_data)

print(nrow(gwas_data))

#I want to extract for all available genes in the GWAS

region<- gwas_data[which(gwas_data$BP>=max(pos-window,0) & gwas_data$BP<=pos+window),] 

#Create a name for the file

setwd("/dir/revision_03_07_25/pwcoco/tissue/data/metabrain/raw")

file_out_snp<- paste0("hippocampus_snp","_", snp, ".txt")
print(file_out_snp)

write.table(region, file_out_snp, sep="\t", quote=F, row.names=F, col.names=T)

}})


#Basal ganglia

lapply(list(ID), function(x){
  
  for(x in ID) {

#Define chromosome, position and window
#These are necessary for extraction and file naming

snp<- x$SNP
  
chr<- x$CHR
print(chr)

pos<- x$BP
print(pos)

window<- 500000

#Read in the GWAS data to extract the regions
gwas_data<- read_table(paste0("/dir/data/metabrain/clean_cd/b37/B37:2021-07-23-basalganglia-EUR-30PCs-chr", chr, ".txt"))

gwas_data<- as.data.frame(gwas_data)

print(nrow(gwas_data))

#I want to extract for all available genes in the GWAS

region<- gwas_data[which(gwas_data$BP>=max(pos-window,0) & gwas_data$BP<=pos+window),] 

#Create a name for the file

setwd("/dir/revision_03_07_25/pwcoco/tissue/data/metabrain/raw")

file_out_snp<- paste0("basalganglia_snp","_", snp, ".txt")
print(file_out_snp)

write.table(region, file_out_snp, sep="\t", quote=F, row.names=F, col.names=T)

}})
