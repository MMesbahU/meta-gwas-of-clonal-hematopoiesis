#!/bin/bash

# Author: Md Mesbah Uddin <mdmesbah@gmail.com>
# Agu 11, 2021
# Run: for chr in {1..22}; do bash run_regenie_step2_LOOCV_Jun23_2021.sh regenie_step2_LOOCV_Jun23_2021.sh UKB200k_11Aug2021 CHIP_ExpandedCHIPDnmtTetAsxlJak "hasCHIP,hasExpandedCHIP,hasDNMT3A,hasTET2,hasASXL1,hasJAK2,hasexpDNMT3A,hasexpTET2,hasexpASXL1" "Age,Age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" "SEX,Ethnic_Background,GenoBatchUKBiLEVE" 400 0.01 20 ${chr} gs://ukbb_v2/projects/muddin/chip_gwas/UKB200k_11Aug2021/step1/NULL_MODEL_pred_modified.list additive gcr.io/ukbb-analyses/regenie:v2.0.2.gz gs://ukbb_v2/projects/muddin/UKB200K/pheno_covar/UKB200k_CHIP_PhenoCovar.11Aug2021.tsv.gz; done

#   
# Jun 23, 2021
# sed 's:output:input:g' NULL_MODEL_pred.list > NULL_MODEL_pred_modified.list

# for chr in {1..22}; do bash run_regenie_step2_LOOCV_Jun23_2021.sh regenie_step2_LOOCV_Jun23_2021.sh UKB200k_23Jun2021 CHIP_ExpandedCHIP "hasCHIP,hasExpandedCHIP,hasDNMT3A,hasTET2,hasASXL1" "AGE_assessment,Age2,Sex_Genetic,GenoBatchUKBiLEVE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" "ever_nerver_smoked,Ethnic_Background" 400 0.01 20 ${chr} gs://ukbb_v2/projects/muddin/chip_gwas/UKB200k_23Jun2021/step1/NULL_MODEL_pred_modified.list additive gcr.io/ukbb-analyses/regenie:v2.0.2.gz gs://ukbb_v2/projects/muddin/UKB200K/pheno_covar/ukb_CHIP_CoVar_v23Jun2021.tsv; done

# User Input
regenie_script=${1} # regenie_step1_LOOCV_Feb27_2021.sh

job_name=${2} # UKB200k_23Jun2021

pheno_name=${3} # lCHIP_VAF10_mCHIP_VAF10

phenoColList=${4} # "hasCHIP,hasExpandedCHIP,hasDNMT3A,hasTET2,hasASXL1"

covar_col_list=${5} # "AGE_assessment,Age2,Sex_Genetic,GenoBatchUKBiLEVE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"

# categorical Covariates
catCovarList=${6} # "ever_nerver_smoked,Ethnic_Background"

geno_block_size=${7} # 400

p_threshold=${8} # 0.01

min_mac=${9} # 20

chr=${10} # 4

step1_list_file=${11} # gs://ukbb_v2/projects/muddin/chip_gwas/UKB200k_23Jun2021/step1/NULL_MODEL_pred_modified.list

genetic_test_name=${12} # additive | dominant | recessive

# Regenie Docker
regenie_docker=${13} # gcr.io/ukbb-analyses/regenie:v2.0.2.gz

    # Pheno and cover combined file
pheno_file=${14} # v

 ## inluded variables
BILLING_PROJECT="ukbb-analyses"

PROJECT_PATH="gs://ukbb_v2/projects/muddin/chip_gwas"
	
	# Imputed GWAS data
geno_path="gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183"

bgen_file="ukb_imp_chr${chr}_v3.bgen"

bgi_file="ukb_imp_chr${chr}_v3.bgen.bgi"

bgen_sample="gs://ukbb_v2/projects/muddin/chip_gwas/ukb7089_imp_chr3_v3_s487395.sample"

    # Step 1 Loco files
# step1_loco_files="gs://ukbb_v2/projects/muddin/chip_gwas/UKB200k_23Jun2021/step1/NULL_MODEL_*.loco"
step1_loco_files="gs://ukbb_v2/projects/muddin/chip_gwas/UKB200k_11Aug2021/step1/NULL_MODEL_*.loco"
# Enable exit on error
set -o errexit

# Launch step 2 and block until completion
echo "Launching REGENIE Step 2"
#--wait --skip 	--input COVAR_FILE=${covar_file} \
# --regions us-central1 \
dsub \
	--project ${BILLING_PROJECT} \
	--provider google-cls-v2 \
	--use-private-address \
	--zone=us-central1-c \
	--disk-type pd-standard \
	--machine-type=e2-highcpu-32 \
	--image ${regenie_docker} \
	--logging ${PROJECT_PATH}/${job_name}/dsub-logs/regenie-step2-chr${chr}-${genetic_test_name} \
	--mount GENO_PATH=${geno_path} \
	--input PHENO_FILE=${pheno_file} \
	--input BGEN_SAMPLE=${bgen_sample} \
	--input STEP1_LOCO_FILE=${step1_loco_files} \
	--input STEP1_LIST_FILE=${step1_list_file} \
	--output-recursive GWAS_OUT_PATH=${PROJECT_PATH}/${job_name}/${genetic_test_name} \
	--env GENETIC_TEST=${genetic_test_name} \
	--env BGEN_FILE=${bgen_file} \
	--env BGI_FILE=${bgi_file} \
	--env PHENO_COL_LIST=${phenoColList} \
	--env COVAR_COL_LIST=${covar_col_list} \
	--env CAT_COVAR_COL_LIST=${catCovarList} \
	--env GENO_BLOCK_SIZE=${geno_block_size} \
	--env P_THRESH=${p_threshold} \
	--env MIN_MAC=${min_mac} \
	--env CHR=${chr} \
	--timeout '5d' \
	--name ${job_name}-step2-chr${chr}-${genetic_test_name} \
	--label "pheno=$(echo ${pheno_name} | awk '{print tolower($0)}')" \
	--label 'step=2' \
	--label 'script=chip_gwas' \
	--label "chromosome=chr${chr}" \
	--script ${regenie_script}


