!/bin/bash

#Extract the regions with evidence of colocalisation from gtex
DIR="/dir/revision_03_07_25/pwcoco/trait/results/passing_h4.txt" 

cd /dir/smr-1.3.1-linux-x86_64/

for j in $(cat $DIR | tr -s ' ' | cut -f15)
do
./smr-1.3.1 --beqtl-summary /dir/data/gtexv7/GTEx_V7_cis_eqtl_summary/Lung \
--query 1 --snp $j --snp-wind 500 --out /dir/revision_03_07_25/pwcoco/tissue/data/gtex/raw/lung_snp_$j
done 

