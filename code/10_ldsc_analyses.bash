module load ldsc/2.0.1-openblas
source activate ldsc

cd /dir/ldsc-master

munge_sumstats.py \
--sumstats /dir/data_repo/clean/respiratory/B37:ASTHMA:2022:EUR.txt \
--N-col NTOT \
--signed-sumstats logOR,0 \
--out /dir/revision_03_07_25/ldsc/munged/B37:ASTHMA:2022:EUR.txt \
--merge-alleles w_hm3.snplist

munge_sumstats.py \
--sumstats /dir/data_repo/clean/neuropsychiatric/B37:BIPOLAR:2025:NOUKB:EUR.txt \
--N-col NTOT \
--signed-sumstats logOR,0 \
--out /dir/revision_03_07_25/ldsc/munged/B37:BIPOLAR:2025:NOUKB:EUR \
--merge-alleles w_hm3.snplist

munge_sumstats.py \
--sumstats /dir/data_repo/clean/neuropsychiatric/B37:SCHIZOPHRENIA:2022:EUR.txt \
--N-col NTOT \
--signed-sumstats logOR,0 \
--out /dir/revision_03_07_25/ldsc/munged/B37:SCHIZOPHRENIA:2022:EUR \
--merge-alleles w_hm3.snplist


#Analysis

ldsc.py \
--rg /dir/revision_03_07_25/ldsc/munged/B37:ASTHMA:2022:EUR.txt.sumstats.gz,/dir/revision_03_07_25/ldsc/munged/B37:BIPOLAR:2025:NOUKB:EUR.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out /dir/revision_03_07_25/ldsc/results/results_asthma_bipolar

ldsc.py \
--rg /dir/revision_03_07_25/ldsc/munged/B37:ASTHMA:2022:EUR.txt.sumstats.gz,/dir/revision_03_07_25/ldsc/munged/B37:SCHIZOPHRENIA:2022:EUR.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out /dir/revision_03_07_25/ldsc/results/results_asthma_schizophrenia

