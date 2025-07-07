#Sensitivity analysis LHC MR

#Import the data
asthma<- vroom::vroom("/dir/data_repo/clean/respiratory/B37:ASTHMA:2022:EUR.txt")
bipolar<- vroom::vroom("/dir/data_repo/clean/neuropsychiatric/B37:BIPOLAR:2025:NOUKB:EUR.txt")

asthma<- as.data.frame(asthma)
bipolar<- as.data.frame(bipolar)

#Rename the columns
library(dplyr)

bipolar_lhc<- bipolar%>%select(SNP,logOR,SE,P,A1,A2,NTOT,CHR,BP)%>%rename(
							LOG_ODDS=	logOR,
							STANDARD_ERROR=SE,
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
