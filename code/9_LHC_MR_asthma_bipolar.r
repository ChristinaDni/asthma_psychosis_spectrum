#Sensitivity analysis LHC MR

#Import the data
asthma<- vroom::vroom("/dir/data_repo/clean/respiratory/B37:ASTHMA:2022:EUR.txt")
bipolar<- vroom::vroom("/dir/data_repo/clean/neuropsychiatric/B37:BIPOLAR:2025:NOUKB:EUR.txt")

asthma<- as.data.frame(asthma)
bipolar<- as.data.frame(bipolar)

#Transfomr the estimates on the liability scale

liability_scale<- function(dat){
get_k <- function(pheno) {
  if (pheno == "bipolar_2025_eur_noukb") {
    return(0.02)
  } else {
    return(0.13)
  }
}

K<- get_k(dat$pheno[1])
print(K)
#Popilation prevalence- for bipolar was extracted from the paper of the GWAS
#For asthma it is tricky due to highly varying prevalence
#Used the approach in the manuscript and used the UKB asthma prevalence considering that the data here are EUR ancestry

get_p <- function(pheno) {
  if (pheno == "bipolar_2025_eur_noukb") {
    return(0.07)
  } else {
    return(0.09)
  }
}

P<- get_p(dat$pheno[1])
print(P)
# P: sample prevalence
#Estimated from a 57833 cases in the GWAS and 722909 controls in bipolar
#Estimated from 121940 cases in the GWAS and 1254131 controls in asthma

t <- qnorm(1 - K)                 # Threshold on the standard normal scale
z <- dnorm(t) / K               # Standard normal PDF at the threshold divided by K
f <- z / sqrt(K * (1 - K)) * sqrt(P * (1 - P)) #Formula to convert on the liability scale

dat$liability_beta <- dat$logOR * f
dat$liability_se<- dat$SE *f
return(dat)
}

asthma<- liability_scale(asthma)
bipolar<- liability_scale(bipolar)

#Rename the columns
library(dplyr)

bipolar_lhc<- bipolar%>%select(SNP,liability_beta,liability_se,P,A1,A2,NTOT,CHR,BP)%>%rename(
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

setwd("/dir/revision_03_07_25/lhc/output2")

## Step 1
trait.names=c("ASTHMA","BIPOLAR")
input.files = list(asthma_lhc,bipolar_lhc)
df = merge_sumstats(input.files,trait.names,LD.filepath,rho.filepath)

## Step 2

SP_list = calculate_SP(df,trait.names,run_ldsc=TRUE,run_MR=TRUE,hm3=hm3,ld=ld,nStep = 2,
                       SP_single=3,SP_pair=50,SNP_filter=10)

## Step 3					   
res = lhc_mr(SP_list, trait.names, paral_method="lapply", nBlock=200)	

####