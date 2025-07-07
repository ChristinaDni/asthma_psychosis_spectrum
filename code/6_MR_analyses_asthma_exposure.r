#MR analysis between asthma, schizophrenia, bipolar
#Forward direction
#This is an updated script: 26/06/2024
#This is a script for the revision: 03/07/2025

###

exposure_file <- "/dir/data_repo/clean/respiratory/B37:ASTHMA:2022:EUR.txt"

outcome_files <- list.files(path="/dir/data_repo/clean/neuropsychiatric", 
pattern="*.txt", full.names=TRUE, recursive=FALSE) 

outcome_files<- outcome_files[c(6,19)]


library(TwoSampleMR)
library(readr)

lapply(list(exposure_file), function(i){
  
  for(i in exposure_file) {
  
  exposure_file<- read_table(i)
  head(exposure_file)
  
iterations= 1:2

for (var in iterations){

  outcome_file<- read_table(outcome_files[[var]])
  head(outcome_file)
  
  #Extract common set of SNPs
  merged<- merge(exposure_file, outcome_file, by.x="SNP", by.y="SNP") 
  
  filename<- merged$pheno.y
  
  exposure<- exposure_file[which(exposure_file$SNP%in%merged$SNP),]
  
  exposure$extracted_from<- filename

	i<- exposure[which(exposure$P<=5e-08),] 

	i<- format_data(i, type = "exposure", snps = NULL, header = TRUE, 
                snp_col = "SNP", beta_col = "logOR", se_col = "SE", effect_allele_col = "A1", 
                other_allele_col = "A2", pval_col = "P", id_col="pheno", samplesize_col="NTOT", phenotype_col="extracted_from", eaf_col="FRQ",
		chr_col= "CHR", pos_col="BP")
				
i$rsid<- i$SNP
i$pval<- i$pval.exposure


exposure_clumped<- ieugwasr::ld_clump(
  dat = i,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/user/work/cd17108/REF_PANEL/EUR",
  plink_bin ="/user/home/cd17108/plink_linux_x86_64_20240818/plink")
  
  
setwd("/dir/revision_03_07_25/mr/results")
if(file.exists("instruments_asthma.txt") == TRUE) {
  write.table(exposure_clumped,"instruments_asthma.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(exposure_clumped, "instruments_asthma.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}
	
out<- read_outcome_data(snps = exposure_clumped$SNP, filename = outcome_files[[var]], 
                        sep = "\t", snp_col = "SNP", beta_col = "logOR", 
                        se_col = "SE", eaf_col = "FRQ", effect_allele_col = "A1", 
                        other_allele_col = "A2", pval_col = "P", id_col="pheno", samplesize_col="NTOT", phenotype_col="pheno", ncase_col="NCAS", ncontrol_col="NCON",
			chr_col="CHR", pos_col="BP")
						
dat<- harmonise_data(exposure_clumped, out) 

setwd("/dir/revision_03_07_25/mr/results")
if(file.exists("harmonised_forward.txt") == TRUE) {
  write.table(dat,"harmonised_forward.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(dat, "harmonised_forward.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

dat_final<- dat[which(dat$mr_keep=="TRUE"),]

res<- mr(dat_final)
or<- generate_odds_ratios(res)

setwd("/dir/revision_03_07_25/mr/results")
if(file.exists("forward_results_asthma_psychoses.txt") == TRUE) {
  write.table(or,"forward_results_asthma_psychoses.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(or, "forward_results_asthma_psychoses.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

} 


  }
  })
