#Sensitivity analysis LHC MR

#Import the data
asthma<- vroom::vroom("/dir/data_repo/clean/respiratory/B37:ASTHMA:2022:EUR.txt")
schizophrenia<- vroom::vroom("/dir/data_repo/clean/neuropsychiatric/B37:SCHIZOPHRENIA:2022:EUR.txt")

asthma<- as.data.frame(asthma)
schizophrenia<- as.data.frame(schizophrenia)

#Transfomr the estimates on the liability scale

liability_scale<- function(dat){
get_k <- function(pheno) {
  if (pheno == "schizophrenia_2022") {
    return(0.01)
  } else {
    return(0.13)
  }
}

K<- get_k(dat$pheno[1])
print(K)
#Population prevalence- for schizophrenia was extracted from the Ripke paper of the GWAS
#For asthma it is tricky due to highly varying prevalence
#Used the approach in the manuscript and used the UKB asthma prevalence considering that the data here are EUR ancestry

get_p <- function(pheno) {
  if (pheno == "schizophrenia_2022") {
    return(0.41)
  } else {
    return(0.09)
  }
}

P<- get_p(dat$pheno[1])
print(P)
# P: sample prevalence
#Estimated from a 53386 cases in the GWAS and 77258 controls in schizophrenia
#Estimated from 121940 cases in the GWAS and 1254131 controls in asthma

t <- qnorm(1 - K)                 # Threshold on the standard normal scale
z <- dnorm(t) / K               # Standard normal PDF at the threshold divided by K
f <- z / sqrt(K * (1 - K)) * sqrt(P * (1 - P)) #Formula to convert on the liability scale

dat$liability_beta <- dat$logOR * f
dat$liability_se<- dat$SE *f
return(dat)
}

asthma<- liability_scale(asthma)
schizophrenia<- liability_scale(schizophrenia)

#Rename the columns
library(dplyr)

schizophrenia_lhc<- schizophrenia%>%select(SNP,liability_beta,liability_se,P,A1,A2,NTOT,CHR,BP)%>%rename(
							LOG_ODDS=	liability_beta,
							STANDARD_ERROR=liability_se,
							N=NTOT
							)
asthma_lhc<- asthma%>%select(SNP,logOR,SE,P,A1,A2,NTOT,CHR,BP)%>%rename(
							LOG_ODDS=	logOR,
							STANDARD_ERROR=SE,
							N=NTOT
							)

#File paths needed for the analysis
LD.filepath = "/dir/files/LDscores_filtered.csv" # LD scores
rho.filepath = "/dir/files/LD_GM2_2prm.csv" # local/SNP-specfic LD scores

ld = "/dir/ldsc-master/eur_w_ld_chr/"
hm3 = "/dir/ldsc-master/eur_w_ld_chr/w_hm3.snplist"

library(lhcMR)

setwd("/dir/revision_03_07_25/lhc/output1")

## Step 1
trait.names=c("ASTHMA","SCHIZOPHRENIA")
input.files = list(asthma_lhc,schizophrenia_lhc)
df = merge_sumstats(input.files,trait.names,LD.filepath,rho.filepath)

## Step 2

SP_list = calculate_SP(df,trait.names,run_ldsc=TRUE,run_MR=TRUE,hm3=hm3,ld=ld,nStep = 2,
                       SP_single=3,SP_pair=50,SNP_filter=10)

## Step 3					   
res = lhc_mr(SP_list, trait.names, paral_method="lapply", nBlock=200)	

####