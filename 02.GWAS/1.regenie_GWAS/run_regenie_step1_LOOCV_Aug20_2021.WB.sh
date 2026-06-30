#!/bin/bash

# Author: Md Mesbah Uddin <mdmesbah@gmail.com>
# Aug 20, 2021

## Run Aug 20, 2021: bash /home/muddin/ukb200k/run_regenie_step1_LOOCV_Aug20_2021.WB.sh /home/muddin/ukb200k/regenie_step1_LOOCV_Aug20_2021.WB.sh eurUKB200k_20Aug2021 "hasCHIP,hasDNMT3A,hasTET2" "Age,Age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" "SEX,GenoBatchUKBiLEVE" gcr.io/ukbb-analyses/regenie:v2.0.2.gz

## Run August 11, 2021: 
## bash /home/muddin/ukb200k/run_regenie_step1_LOOCV_Jun23_2021.sh /home/muddin/ukb200k/regenie_step1_LOOCV_Jun23_2021.sh UKB200k_11Aug2021 "hasCHIP,hasExpandedCHIP,hasDNMT3A,hasTET2,hasASXL1,hasJAK2,hasexpDNMT3A,hasexpTET2,hasexpASXL1" "Age,Age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" "SEX,Ethnic_Background,GenoBatchUKBiLEVE" gcr.io/ukbb-analyses/regenie:v2.0.2.gz

#####
## bash /home/muddin/ukb200k/run_regenie_step1_LOOCV_Jun23_2021.sh /home/muddin/ukb200k/regenie_step1_LOOCV_Jun23_2021.sh UKB200k_23Jun2021 "hasCHIP,hasExpandedCHIP,hasDNMT3A,hasTET2,hasASXL1" "AGE_assessment,Age2,Sex_Genetic,GenoBatchUKBiLEVE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" "ever_nerver_smoked,Ethnic_Background" gcr.io/ukbb-analyses/regenie:v2.0.2.gz

# This script runs a UK Biobank GWAS using REGENIE on GCP using dsub. 
# First, please install dsub: https://github.com/DataBiosphere/dsub#install-dsub .

# This script is written with reference to the Broad-wide UK Biobank genetic
# data stored under gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183 
## bash run_regenie_step1_LOOCV.sh regenie_step1_LOOCV.sh JAK2_Sep8_2020 JAK2 
BILLING_PROJECT="ukbb-analyses"
PROJECT_PATH="gs://ukbb_v2/projects/muddin/chip_gwas"

## Input files
REGENIE_SCRIPT=${1}
job_name=${2} # "CHIP_Jun23_2021"
	# Load Genotype Plink files 
plink_files="gs://ukbb_v2/projects/muddin/ukb_v2_data/ukb_v2_GT_chr1_22*"
	
	# SNPs to keep
	## script used: "gs://ukbb_v2/projects/muddin/ukb_v2_data/ukb_qc_step1/Jun24_2021/combine_ukb_BEDs.sh"
snps_2_keep="gs://ukbb_v2/projects/muddin/ukb_v2_data/ukb_qc_step1/Jun24_2021/qc_pass.chr1_22.24jun2021.snplist" #"gs://ukbb_v2/projects/muddin/ukb_v2_data/ukb_qc_step1/qc_pass.chr1_22.snplist"
	# Samples to keep: one sample with >10% missing genotype "4523546"
samples_2_keep="gs://ukbb_v2/projects/muddin/ukb_v2_data/ukb_qc_step1/Jun24_2021/qc_pass.chr1_22.24jun2021.id" # "gs://ukbb_v2/projects/muddin/ukb_v2_data/ukb_qc_step1/qc_pass.chr1_22.id"

# Load Pheno file 
# pheno_file="gs://ukbb_v2/projects/muddin/chip_gwas/pheno_file_ukb_CHIP.tsv" 
# pheno_file="gs://ukbb_v2/projects/muddin/chip_gwas/traits/WhiteBrt1001.ML_CHIP.UKB_pheno.27Feb2021.tsv"
# pheno_file="gs://ukbb_v2/projects/muddin/UKB200K/pheno_covar/ukb_CHIP_CoVar_v23Jun2021.tsv"
# pheno_file="gs://ukbb_v2/projects/muddin/UKB200K/pheno_covar/UKB200k_CHIP_PhenoCovar.11Aug2021.tsv.gz"
pheno_file="gs://ukbb_v2/projects/muddin/UKB200K/pheno_covar/UKB200k_CHIP_PhenoCovar.20Aug2021.WB.tsv.gz"

# space/tab separated file with header: FID IID Y1 Y2 Y3 Y4 Y5 
#phenoCols="CHIP,VAF,DNMT3A,TET2,ASXL1,JAK2"  # comma separated list of pheno to include, e.g. Y1,Y3,Y5  
# CHIP,VAF10,DNMT3A,TET2,ASXL1,JAK2,TP53,PPM1D,SRSF2,GNB1,SF3B1
phenoColList=${3}
# Covar files
#covar_file="gs://ukbb_v2/projects/muddin/chip_gwas/covar_file_ukb_CHIP_sex01.tsv"

# numerical Covariates
covarColList=${4}

# categorical Covariates
catCovarList=${5} #"Ever_Smoked,Race"

# Bin size
bsize=1000
# num_threads= | nproc --all
# out_file_name="CHIP_GWAS_v1"

# Regenie Docker: gcr.io/ukbb-analyses/regenie:v1.0.6.0.gz
regenie_docker=${6} # "gcr.io/ukbb-analyses/regenie:v1.0.5.8"

# Launch step 1 and block until completion
echo "Launching REGENIE Step 1"
dsub \
	--project ${BILLING_PROJECT} \
	--provider google-cls-v2 \
	--use-private-address \
	--regions us-central1 us-east1 us-west1 \
	--disk-type pd-standard \
	--disk-size 200 \
	--machine-type=e2-highcpu-32 \
	--image ${regenie_docker} \
	--skip \
	--logging ${PROJECT_PATH}/${job_name}/dsub-logs/regenie-step1 \
	--input PHENO_FILE=${pheno_file} \
	--input PLINK_FILES=${plink_files} \
	--input KEEP_SAMPLE=${samples_2_keep} \
	--input KEEP_SNPs=${snps_2_keep} \
	--output-recursive OUTPUT_PATH=${PROJECT_PATH}/${job_name}/step1 \
	--env PHENO_COL_LIST=${phenoColList} \
	--env COVAR_COL_LIST=${covarColList} \
	--env CAT_COVAR_COL_LIST=${catCovarList} \
	--env BIN_SIZE=${bsize} \
	--timeout '1w' \
	--name ${job_name}-step1 \
	--label 'step=1' \
	--label 'script=chip_gwas_null_model' \
	--script ${REGENIE_SCRIPT}


