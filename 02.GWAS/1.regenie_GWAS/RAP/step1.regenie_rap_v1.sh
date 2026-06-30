#!/bin/bash

## upload files: /home/unix/muddin/.local/bin/dx upload /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/step1.regenie_rap.sh --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/ch_pheno_8Mar2023/
################
cohort=${1} #"ukb200k"
PLINK_Prefix=${2} #"ukb_v2_GT_chr1_22"
KEEP_SNPs=${3} #"qc_pass_snps.ukb_v2.list"
KEEP_SAMPLE=${4} # "qc_pass_ukb_all_samples.id"
PHENO_FILE=${5} #"CH_phenoCovar.ukb200k_ukb250k_N421316.27cols.8Mar2023.tsv.gz"
PHENO_COL_LIST=${6} # "hasCH,hasCHvaf05,hasCHvaf10,hasDNMT3A,hasTET2,hasDTA"
COVAR_COL_LIST=${7} #"Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10}"
CAT_COVAR_COL_LIST=${8} # "Genetic_Sex,Ethnic_Background,GenoBatch"
BIN_SIZE=${9} # 1000
docker_image=${10} # mmuddin2020/regenie:v3.2.5.2.gz
minCase=${11} # minimum case count
#OUTPUT_PATH="/ukbb_CH/ch_gwas_2023/step1"

# cd /home/dnanexus/
# ls -lhrt

# pwd

# ls -lhrt /tmp

# ls -lhrt /mnt/project

# find . | grep -E 'ukb21007_c5_b0_v1|CH_phenoCovar'

# docker pull ${docker_image}
##
docker run -w /tmp -v /home/dnanexus/out/out:/tmp -v /mnt/project:/mnt/project ${docker_image} regenie \
	--step 1 \
	--loocv \
	--bed ${PLINK_Prefix} \
	--extract ${KEEP_SNPs} \
	--keep ${KEEP_SAMPLE} \
	--phenoFile ${PHENO_FILE} \
	--phenoColList ${PHENO_COL_LIST} \
	--covarFile ${PHENO_FILE} \
	--covarColList ${COVAR_COL_LIST} \
	--catCovarList ${CAT_COVAR_COL_LIST} \
	--bt \
	--minCaseCount ${minCase} \
	--bsize ${BIN_SIZE} \
	--lowmem \
	--threads $(nproc --all) \
	--out NULL_MODEL.${cohort}


