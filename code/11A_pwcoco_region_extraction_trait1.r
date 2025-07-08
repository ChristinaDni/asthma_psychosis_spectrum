#Import the packages
library(dplyr)
library(TwoSampleMR)


#Import the asthma GWAS
asthma_gwas<- vroom::vroom("/dir/data_repo/clean/respiratory/B37:ASTHMA:2022:EUR.txt")
asthma_gwas<- as.data.frame(asthma_gwas)

#Remove the wider MHC region

asthma_gwas<- asthma_gwas[!(asthma_gwas$CHR=6&(asthma_gwas$BP>= 28477797 & asthma_gwas$BP<= 33448354)),]

#Extract the asthma instruments

i<- asthma_gwas[which(asthma_gwas$P<=5e-08),] 

i<- format_data(i, type = "exposure", snps = NULL, header = TRUE, 
                snp_col = "SNP", beta_col = "logOR", se_col = "SE", effect_allele_col = "A1", 
                other_allele_col = "A2", pval_col = "P", id_col="pheno", samplesize_col="NTOT", 
				ncase_col = "NCAS",
				ncontrol_col = "NCON",
				phenotype_col="pheno", eaf_col="FRQ",
				chr_col = "CHR",
				pos_col = "BP")
				
i$rsid<- i$SNP
i$pval<- i$pval.exposure

asthma_final<- ieugwasr::ld_clump(
  dat = i,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/dir/REF_PANEL/EUR",
  plink_bin ="/dir/plink_linux_x86_64_20240818/plink") #154
    
write.csv(asthma_final, "/dir/revision_03_07_25/pwcoco/trait/data/asthma_instruments.csv")
  
###

#Import the bipolar GWAS
bip_gwas<- vroom::vroom("/dir/data_repo/clean/neuropsychiatric/B37:BIPOLAR:2025:NOUKB:EUR.txt")
bip_gwas<- as.data.frame(bip_gwas)

#Remove the wider MHC region

bip_gwas<- bip_gwas[!(bip_gwas$CHR=6&(bip_gwas$BP>= 28477797 & bip_gwas$BP<= 33448354)),]

#Extract the bip instruments

i<- bip_gwas[which(bip_gwas$P<=5e-08),] 

i<- format_data(i, type = "exposure", snps = NULL, header = TRUE, 
                snp_col = "SNP", beta_col = "logOR", se_col = "SE", effect_allele_col = "A1", 
                other_allele_col = "A2", pval_col = "P", id_col="pheno", samplesize_col="NTOT", 
				ncase_col = "NCAS",
				ncontrol_col = "NCON",
				phenotype_col="pheno", eaf_col="FRQ",
				chr_col = "CHR",
				pos_col = "BP")
				
i$rsid<- i$SNP
i$pval<- i$pval.exposure

bip_final<- ieugwasr::ld_clump(
  dat = i,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/dir/REF_PANEL/EUR",
  plink_bin ="/dir/plink_linux_x86_64_20240818/plink") #66
    
write.csv(bip_final, "/dir/revision_03_07_25/pwcoco/trait/data/bip_instruments.csv") 

###

#Import the schizophrenia GWAS
scz_gwas<- vroom::vroom("/dir/data_repo/clean/neuropsychiatric/B37:SCHIZOPHRENIA:2022:EUR.txt")
scz_gwas<- as.data.frame(scz_gwas)

#Remove the wider MHC region

scz_gwas<- scz_gwas[!(scz_gwas$CHR=6&(scz_gwas$BP>= 28477797 & scz_gwas$BP<= 33448354)),]

#Extract the scz instruments

i<- scz_gwas[which(scz_gwas$P<=5e-08),] 

i<- format_data(i, type = "exposure", snps = NULL, header = TRUE, 
                snp_col = "SNP", beta_col = "logOR", se_col = "SE", effect_allele_col = "A1", 
                other_allele_col = "A2", pval_col = "P", id_col="pheno", samplesize_col="NTOT", 
				ncase_col = "NCAS",
				ncontrol_col = "NCON",
				phenotype_col="pheno", eaf_col="FRQ",
				chr_col = "CHR",
				pos_col = "BP")
				
i$rsid<- i$SNP
i$pval<- i$pval.exposure

scz_final<- ieugwasr::ld_clump(
  dat = i,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/dir/REF_PANEL/EUR",
  plink_bin ="/dir/plink_linux_x86_64_20240818/plink") #175
    
write.csv(scz_final, "/dir/revision_03_07_25/pwcoco/trait/data/scz_instruments.csv")
  
###

#Merge them to clump again

all<- rbind(asthma_final,scz_final)

#Change the id so that they can be clumped together

all$id.exposure<- 1
all$id<- 1

final<- ieugwasr::ld_clump(
  dat = all,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/dir/REF_PANEL/EUR",
  plink_bin ="/dir/plink_linux_x86_64_20240818/plink") #395 variants
  
 #Split again so that they  can be extracted
 
asthma<- final[which(final$exposure=="asthma_GBMI_2022"),]
  
scz<- final[which(final$exposure=="schizophrenia_2022"),]

###

#Extract regions from each gwas

#asthma

map<- asthma %>% split(x=asthma, f=asthma$SNP)


lapply(list(map), function(x){
  
  for(x in map) {
  
  tryCatch({

#Here creating a document for each SNP for which I extract- ie they are the instruments

region<- x[c("SNP","chr.exposure","pos.exposure")]
region$TopSNP<- paste0((region$SNP),
                     ":",
                     (region$chr.exposure),
                     ":",
                     (region$pos.exposure))

setwd("/dir/revision_03_07_25/pwcoco/trait/data")
if(file.exists("topregions_asthma_schizophrenia.txt") == TRUE) {
  write.table(region,"topregions_asthma_schizophrenia.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_asthma_schizophrenia.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#This indicator will be useful to name each region dataframe
snp_list<- paste0(region$SNP,":",region$chr.exposure, ":", region$pos.exposure)
snp_list<- as.data.frame(snp_list)

#Define the extraction window
window<- 500000

#Extract on the chromosome
all_chr<- asthma_gwas[which(asthma_gwas$CHR==region$chr.exposure),]

#Extract on the BP window
all_pos<- all_chr[which(all_chr$BP>=max(region$pos.exposure-window, 0) & all_chr$BP<=region$pos.exposure+window),]

#Keep all the necessary information 
variants_region <- all_pos %>% 
    dplyr::select("SNP", "A1", "A2", "FRQ", "logOR", "SE", "P", "NTOT", "NCAS") 
    	
variants_region<- as.data.frame(variants_region)

print(nrow(variants_region))

setwd ("/dir/revision_03_07_25/pwcoco/trait/data/schizophrenia/trait1")

#Save the output

snp<- snp_list$snp_list

file_out_snp<- paste0("asthma_", snp, ".txt")
for(o in 1:length(snp)) {

write.table(variants_region, file_out_snp[o], sep="\t", quote=F, row.names=F, col.names=T)}

}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}})


#schizophrenia
map<- scz %>% split(x=scz, f=scz$SNP)

lapply(list(map), function(x){
  
  for(x in map) {
  
  tryCatch({

#Here creating a document for each SNP for which I extract- ie they are the instruments

region<- x[c("SNP","chr.exposure","pos.exposure")]
region$TopSNP<- paste0((region$SNP),
                     ":",
                     (region$chr.exposure),
                     ":",
                     (region$pos.exposure))

setwd("/dir/revision_03_07_25/pwcoco/trait/data")
if(file.exists("topregions_scz_asthma.txt") == TRUE) {
  write.table(region,"topregions_scz_asthma.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_scz_asthma.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#This indicator will be useful to name each region dataframe
snp_list<- paste0(region$SNP,":",region$chr.exposure, ":", region$pos.exposure)
snp_list<- as.data.frame(snp_list)

#Define the extraction window
window<- 500000

#Extract on the chromosome
all_chr<- scz_gwas[which(scz_gwas$CHR==region$chr.exposure),]

#Extract on the BP window
all_pos<- all_chr[which(all_chr$BP>=max(region$pos.exposure-window, 0) & all_chr$BP<=region$pos.exposure+window),]

#Keep all the necessary information 
variants_region <- all_pos %>% 
    dplyr::select("SNP", "A1", "A2", "FRQ", "logOR", "SE", "P", "NTOT", "NCAS") 
    	
variants_region<- as.data.frame(variants_region)

print(nrow(variants_region))

setwd ("/dir/revision_03_07_25/pwcoco/trait/data/schizophrenia/trait1")

#Save the output

snp<- snp_list$snp_list

file_out_snp<- paste0("scz_", snp, ".txt")
for(o in 1:length(snp)) {

write.table(variants_region, file_out_snp[o], sep="\t", quote=F, row.names=F, col.names=T)}

}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}})

###

#Repeat the same process for bipolar

#Merge them to clump again

all<- rbind(asthma_final,bip_final)

#Change the id so that they can be clumped together

all$id.exposure<- 1
all$id<- 1

final<- ieugwasr::ld_clump(
  dat = all,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/dir/REF_PANEL/EUR",
  plink_bin ="/dir/plink_linux_x86_64_20240818/plink") #395 variants
  
 #Split again so that they  can be extracted
 
asthma<- final[which(final$exposure=="asthma_GBMI_2022"),]
  
bip<- final[which(final$exposure=="bipolar_2025_eur_noukb"),]

###

#Extract regions from each gwas

#asthma

map<- asthma %>% split(x=asthma, f=asthma$SNP)


lapply(list(map), function(x){
  
  for(x in map) {
  
  tryCatch({

#Here creating a document for each SNP for which I extract- ie they are the instruments

region<- x[c("SNP","chr.exposure","pos.exposure")]
region$TopSNP<- paste0((region$SNP),
                     ":",
                     (region$chr.exposure),
                     ":",
                     (region$pos.exposure))

setwd("/dir/revision_03_07_25/pwcoco/trait/data")
if(file.exists("topregions_asthma_bipolar.txt") == TRUE) {
  write.table(region,"topregions_asthma_bipolar.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_asthma_bipolar.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#This indicator will be useful to name each region dataframe
snp_list<- paste0(region$SNP,":",region$chr.exposure, ":", region$pos.exposure)
snp_list<- as.data.frame(snp_list)

#Define the extraction window
window<- 500000

#Extract on the chromosome
all_chr<- asthma_gwas[which(asthma_gwas$CHR==region$chr.exposure),]

#Extract on the BP window
all_pos<- all_chr[which(all_chr$BP>=max(region$pos.exposure-window, 0) & all_chr$BP<=region$pos.exposure+window),]

#Keep all the necessary information 
variants_region <- all_pos %>% 
    dplyr::select("SNP", "A1", "A2", "FRQ", "logOR", "SE", "P", "NTOT", "NCAS") 
    	
variants_region<- as.data.frame(variants_region)

print(nrow(variants_region))

setwd ("/dir/revision_03_07_25/pwcoco/trait/data/bipolar/trait1")

#Save the output

snp<- snp_list$snp_list

file_out_snp<- paste0("asthma_", snp, ".txt")
for(o in 1:length(snp)) {

write.table(variants_region, file_out_snp[o], sep="\t", quote=F, row.names=F, col.names=T)}

}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}})


#Bipolar
map<- bip %>% split(x=bip, f=bip$SNP)

lapply(list(map), function(x){
  
  for(x in map) {
  
  tryCatch({

#Here creating a document for each SNP for which I extract- ie they are the instruments

region<- x[c("SNP","chr.exposure","pos.exposure")]
region$TopSNP<- paste0((region$SNP),
                     ":",
                     (region$chr.exposure),
                     ":",
                     (region$pos.exposure))

setwd("/dir/revision_03_07_25/pwcoco/trait/data")
if(file.exists("topregions_bip_asthma.txt") == TRUE) {
  write.table(region,"topregions_bip_asthma.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_bip_asthma.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#This indicator will be useful to name each region dataframe
snp_list<- paste0(region$SNP,":",region$chr.exposure, ":", region$pos.exposure)
snp_list<- as.data.frame(snp_list)

#Define the extraction window
window<- 500000

#Extract on the chromosome
all_chr<- bip_gwas[which(bip_gwas$CHR==region$chr.exposure),]

#Extract on the BP window
all_pos<- all_chr[which(all_chr$BP>=max(region$pos.exposure-window, 0) & all_chr$BP<=region$pos.exposure+window),]

#Keep all the necessary information 
variants_region <- all_pos %>% 
    dplyr::select("SNP", "A1", "A2", "FRQ", "logOR", "SE", "P", "NTOT", "NCAS") 
    	
variants_region<- as.data.frame(variants_region)

print(nrow(variants_region))

setwd ("/dir/revision_03_07_25/pwcoco/trait/data/bipolar/trait1")

#Save the output

snp<- snp_list$snp_list

file_out_snp<- paste0("bip_", snp, ".txt")
for(o in 1:length(snp)) {

write.table(variants_region, file_out_snp[o], sep="\t", quote=F, row.names=F, col.names=T)}

}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}})

