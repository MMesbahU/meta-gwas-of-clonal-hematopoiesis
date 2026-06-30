#!/bin/bash

##################################################################
############ Author: Md Mesbah Uddin <mdmesbah@gmail.com>  #######
############ Aug 26, 2023				   #######
#################################################################

########## 250k ############
## ## AoU 250k CHIP
# Multi-ancestry:
## while read lines; do bash /home/jupyter/ch_gwas/aou250k/run_regenie_step2_LOOCV.AoU_2023.preemptible.sh CH_DTA gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/gwas gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/scripts_v7/regenie_step2_LOOCV.AoU_2023.sh AoU_multiAnc_2023 "hasCH,hasCHvaf10,hasDTA,hasDDR,hasSF,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1" "Age_biosample_collection,sqrAge_biosample_collection,PC{1:10}" "dragen_sex_ploidy,ancestry_pred_other,site_id,Batch_CH" 400 0.01 100 $(echo ${lines} | awk '{print $1}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/gwas/AoU_multiAnc_2023/step1/NULL_MODEL_pred_modified.list gcr.io/bick-aps2/regenie:v3.2.8 gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/my_CH/phenoCH_AoU250k.fid0_iid.qcd_myeloidCA_rel_NA.26Aug2023.tsv 0.3 $(echo ${lines} | awk '{print $1":"$2}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/gwas/AoU_multiAnc_2023/step1/NULL_MODEL gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold/bgen/acaf_threshold.chr$(echo ${lines} | awk '{print $1}').bgen gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold/bgen/acaf_threshold.chr$(echo ${lines} | awk '{print $1}').sample e2-highcpu-32; done </home/jupyter/ch_gwas/hg38_chrom_interval.tsv

###########################


### Notes:
## modify begn sample files by replacing fid by "0" to batch with array genotypes, null step, and phenotype file.
## awk '(NR<3){print $0}(NR>2){print 0,$2,$3}' aou-alpha3-filtered-chr22.sample > aou-alpha3-filtered-chr22_fid0_iid.sample
# test run worked: regenie --step 2 --bgen aou-alpha3-filtered-chr22.bgen --sample aou-alpha3-filtered-chr22_fid0_iid.sample --bt --htp AFR --phenoFile aou_ch_phenotype.AFR.fid0_iid.14Mar2023.tsv.gz --phenoColList hasCH,hasDNMT3A,hasTET2,hasDTA --covarFile aou_ch_phenotype.AFR.fid0_iid.14Mar2023.tsv.gz --covarColList Age_at_survey,Age_sqr,PC{1:10} --catCovarList Sex_at_birth --bt --bsize 400 --firth --ref-first --approx --gz --write-samples --print-pheno --pThresh 0.01 --pred NULL_MODEL_pred_nopath.list --minMAC 40 --minINFO 0.3 --threads $(nproc --all) --range chr22:28619558-29203497 --out test.chr22.afr
#########

## AoU 100k CHIP
# Multi-ancestry: 
## while read lines; do bash run_regenie_step2_LOOCV.AoU_2023.preemptible.sh gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/scripts/regenie_step2_LOOCV.AoU_2023.sh AoU_multiAnc_2023 "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_survey,Age_sqr,PC{1:10}" "Sex_at_birth,Ancestry" 400 0.01 40 $(echo ${lines} | awk '{print $1}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/AoU_multiAnc_2023/step1/NULL_MODEL_pred_modified.list gcr.io/bick-aps2/ghcr.io/rgcgithub/regenie/regenie:v3.2.4.gz gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/aou_ch_phenotype.fid0_iid.14Mar2023.tsv.gz 0.30 $(echo ${lines} | awk '{print $1":"$2}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/AoU_multiAnc_2023/step1/NULL_MODEL gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr$(echo ${lines} | awk '{print $1}').bgen gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr$(echo ${lines} | awk '{print $1}').sample e2-highcpu-32; done <hg38_chrom_interval.tsv

## EUR: aou_ch_phenotype.EUR.fid0_iid.14Mar2023.tsv.gz	
# while read lines; do bash run_regenie_step2_LOOCV.AoU_2023.preemptible.sh gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/scripts/regenie_step2_LOOCV.AoU_2023.sh AoU_EUR_2023 "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_survey,Age_sqr,PC{1:10}" Sex_at_birth 400 0.01 40 $(echo ${lines} | awk '{print $1}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/AoU_EUR_2023/step1/NULL_MODEL_pred_modified.list gcr.io/bick-aps2/ghcr.io/rgcgithub/regenie/regenie:v3.2.4.gz gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/aou_ch_phenotype.EUR.fid0_iid.14Mar2023.tsv.gz 0.30 $(echo ${lines} | awk '{print $1":"$2}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/AoU_EUR_2023/step1/NULL_MODEL gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr$(echo ${lines} | awk '{print $1}').bgen gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr$(echo ${lines} | awk '{print $1}').sample n2-standard-8; done <hg38_chrom_interval.tsv

## AFR:
## while read lines; do bash run_regenie_step2_LOOCV.AoU_2023.preemptible.sh gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/scripts/regenie_step2_LOOCV.AoU_2023.sh AoU_AFR_2023 "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_survey,Age_sqr,PC{1:10}" Sex_at_birth 400 0.01 40 $(echo ${lines} | awk '{print $1}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/AoU_AFR_2023/step1/NULL_MODEL_pred_modified.list gcr.io/bick-aps2/ghcr.io/rgcgithub/regenie/regenie:v3.2.4.gz gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/aou_ch_phenotype.AFR.fid0_iid.14Mar2023.tsv.gz 0.30 $(echo${lines} | awk '{print $1":"$2}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/AoU_AFR_2023/step1/NULL_MODEL gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr$(echo ${lines} | awk '{print $1}').bgen gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr$(echo ${lines} | awk '{print $1}').sample n2-standard-8; done <hg38_chrom_interval.tsv 1>>job_log.txt 2>>job_err.txt &

## AMR
# while read lines; do bash run_regenie_step2_LOOCV.AoU_2023.preemptible.sh gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/scripts/regenie_step2_LOOCV.AoU_2023.sh AoU_AMR_2023 "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_survey,Age_sqr,PC{1:10}" Sex_at_birth 400 0.01 40 $(echo ${lines} | awk '{print $1}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/AoU_AMR_2023/step1/NULL_MODEL_pred_modified.list gcr.io/bick-aps2/ghcr.io/rgcgithub/regenie/regenie:v3.2.4.gz gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/aou_ch_phenotype.AMR.fid0_iid.14Mar2023.tsv.gz 0.30 $(echo ${lines} | awk '{print $1":"$2}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/AoU_AMR_2023/step1/NULL_MODEL gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr$(echo ${lines} | awk '{print $1}').bgen gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr$(echo ${lines} | awk '{print $1}').sample n2-standard-8; done <hg38_chrom_interval.tsv 1>>job_log.txt 2>>job_err.txt &
#########################################

## rerun failed jobs
# while read lines; do bash run_regenie_step2_LOOCV.AoU_2023.preemptible.sh gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/scripts/regenie_step2_LOOCV.AoU_2023.sh $(echo ${lines} | awk '{print $1}') "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_survey,Age_sqr,PC{1:10}" Sex_at_birth 400 0.01 40 $(echo ${lines} | awk '{print $2}' | awk -F':''{print $1}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/$(echo ${lines} | awk '{print $1}')/step1/NULL_MODEL_pred_modified.list gcr.io/bick-aps2/ghcr.io/rgcgithub/regenie/regenie:v3.2.4.gz gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/aou_ch_phenotype.$(echo ${lines} | awk '{print $1}' | sed -e 's:AoU_::g' -e 's:_2023::g').fid0_iid.14Mar2023.tsv.gz 0.30 $(echo ${lines} | awk '{print $2}') gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/$(echo ${lines} | awk '{print $1}')/step1/NULL_MODEL gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr$(echo ${lines} | awk '{print $2}' | awk -F':''{print $1}').bgen gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr$(echo ${lines} | awk '{print $2}' | awk -F':''{print $1}').sample n2-standard-8; done < <(grep -v '22:1-10163693'  missed_run.sh) 1>>rerun_failed_job_log.txt 2>>rerun_failed_job_err.txt &

## 
pheno_name=${1} # "CH_DTA"

### AoU specific
#
PROJECT_PATH=${2} # "gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas"

##
AOU_NETWORK=network

# local
AOU_SUBNETWORK=subnetwork
DSUB_USER_NAME="$(echo "${OWNER_EMAIL}" | cut -d@ -f1)"

## e2-highcpu-32: spot $0.23744 vs $0.791488
## https://cloud.google.com/compute/vm-instance-pricing
##########
# User Input
REGENIE_SCRIPT=${3} # gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/scripts/regenie_step2_LOOCV.AoU_2023.sh

job_name=${4} # AoU_multiAnc_2023

phenoColList=${5} # hasCH,hasDNMT3A,hasTET2,hasDTA

covar_col_list=${6} # Age_at_survey,Age_sqr,PC{1:10}

# categorical Covariates
catCovarList=${7} # Sex_at_birth,Ancestry

geno_block_size=${8} # 400

p_threshold=${9} # 0.01

min_mac=${10} # 40

chr=${11} # 4

step1_list_file=${12} # 

# Regenie Docker
REGENIE_DOCKER=${13} # gcr.io/bick-aps2/ghcr.io/rgcgithub/regenie/regenie:v3.2.4.gz

    # Pheno and cover combined file
pheno_file=${14} # pheno_file="gs://ukbb_v2/projects/muddin/UKB200K/pheno_covar/UKB200k_CHIP_PhenoCovar.20Aug2021.WB.tsv.gz"

# imputation rsq>0.3
INFO=${15}

chr_range=${16}
##
   # Step 1 Loco files
step1_loco_files="${17}_*.loco" # gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/AoU_AFR_2023/step1/NULL_MODEL 
	# WGS bgen from Bick lab
bgen_file=${18} # gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr${chr}.bgen

bgen_sample=${19} # gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/data/aou/geno/20221116/aou-alpha3-filtered-chr${chr}.sample

machine_type=${20} # n2-standard-8

# Enable exit on error
set -o errexit

# Launch step 2 and block until completion
echo -e "\n"
echo -e "Launching REGENIE Step 2: ${chr_range}\n"

#--wait --skip 	--input COVAR_FILE=${covar_file} \
# --regions us-central1 \
# MAX_PREEMPTION=6
# MAX_Retries=7
####
dsub \
	--user-project ${GOOGLE_PROJECT} \
	--project ${GOOGLE_PROJECT} \
	--network ${AOU_NETWORK} \
	--subnetwork ${AOU_SUBNETWORK} \
	--service-account $(gcloud config get-value account) \
	--user ${DSUB_USER_NAME} \
	--provider google-cls-v2 \
	--regions us-central1 \
	--disk-type pd-standard \
	--disk-size 200 \
	--machine-type=${machine_type} \
	--image ${REGENIE_DOCKER} \
	--logging ${PROJECT_PATH}/${job_name}/dsub-logs/step2/regenie-step2-chr${chr}_$(echo ${chr_range} | awk -F ':' '{print $2}') \
	--input PHENO_FILE=${pheno_file} \
	--input BGEN_SAMPLE=${bgen_sample} \
	--input BGEN_FILE=${bgen_file} \
	--input STEP1_LOCO_FILE=${step1_loco_files} \
	--input STEP1_LIST_FILE=${step1_list_file} \
	--output-recursive GWAS_OUT_PATH=${PROJECT_PATH}/${job_name}/step2 \
	--env PHENO_COL_LIST=${phenoColList} \
	--env COVAR_COL_LIST=${covar_col_list} \
	--env CAT_COVAR_COL_LIST=${catCovarList} \
	--env GENO_BLOCK_SIZE=${geno_block_size} \
	--env P_THRESH=${p_threshold} \
	--env MIN_MAC=${min_mac} \
	--env CHR=${chr} \
	--env COHORT=${job_name} \
	--env INFO=${INFO} \
	--env CHR_RANGE=${chr_range} \
	--timeout '1w' \
	--name ${job_name}-step2-chr${chr}-$(echo ${chr_range} | awk -F ':' '{print$2}') \
	--label "pheno=$(echo ${pheno_name} | awk '{print tolower($0)}')" \
	--label 'step=2' \
	--label 'script=chip_gwas' \
	--label "chromosome=chr${chr}" \
	--script ${REGENIE_SCRIPT}

#####
echo -e "\n"

