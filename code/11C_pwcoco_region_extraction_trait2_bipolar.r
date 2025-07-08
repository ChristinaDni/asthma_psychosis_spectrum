library(dplyr)
library(readr)

setwd("/dir/revision_03_07_25/pwcoco/trait/data")
data<- read.table("topregions_asthma_bipolar.txt", sep="\t", header=T)

ID<- data %>% split(x=data, f=data$SNP)

#Read in the GWAS data to extract the regions
gwas_data<- read_table("/dir/data_repo/clean/neuropsychiatric/B37:BIPOLAR:2025:NOUKB:EUR.txt")

gwas_data<- as.data.frame(gwas_data)

#Here bipolar

lapply(list(ID), function(x){
  
  for(x in ID) {

#Define chromosome, position and window
#These are necessary for extraction and file naming

snp<- x$SNP
  
chr<- x$chr.exposure
print(chr)

pos<- x$pos.exposure
print(pos)

window<- 500000

region_chr<-gwas_data[which(gwas_data$CHR==chr),] 
region<- region_chr[which(region_chr$BP>=max(pos-window,0) & region_chr$BP<=pos+window),] 

#Select the columns I need for pwcoco
region<- region%>%dplyr::select("SNP","A1","A2", "FRQ", "logOR", "SE", "P", "NTOT", "NCAS")


#Create a name for the file

setwd("/dir/revision_03_07_25/pwcoco/trait/data/bipolar/trait2")

file_out_snp<- paste0("bip","_", snp, ".txt")
print(file_out_snp)

write.table(region, file_out_snp, sep="\t", quote=F, row.names=F, col.names=T)

}})

setwd("/dir/revision_03_07_25/pwcoco/trait/data")
data<- read.table("topregions_bip_asthma.txt", sep="\t", header=T)

ID<- data %>% split(x=data, f=data$SNP)

#Read in the GWAS data to extract the regions
gwas_data<- read_table("/dir/data_repo/clean/respiratory/B37:ASTHMA:2022:EUR.txt")

gwas_data<- as.data.frame(gwas_data)

#Here bipolar

lapply(list(ID), function(x){
  
  for(x in ID) {

#Define chromosome, position and window
#These are necessary for extraction and file naming

snp<- x$SNP
  
chr<- x$chr.exposure
print(chr)

pos<- x$pos.exposure
print(pos)

window<- 500000

region_chr<-gwas_data[which(gwas_data$CHR==chr),] 
region<- region_chr[which(region_chr$BP>=max(pos-window,0) & region_chr$BP<=pos+window),] 

#Select the columns I need for pwcoco
region<- region%>%dplyr::select("SNP","A1","A2", "FRQ", "logOR", "SE", "P", "NTOT", "NCAS")


#Create a name for the file

setwd("/dir/revision_03_07_25/pwcoco/trait/data/bipolar/trait2")

file_out_snp<- paste0("asthma","_", snp, ".txt")
print(file_out_snp)

write.table(region, file_out_snp, sep="\t", quote=F, row.names=F, col.names=T)

}})
