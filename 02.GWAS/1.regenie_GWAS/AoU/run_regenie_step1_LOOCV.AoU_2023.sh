#!/bin/bash

# Author: Md Mesbah Uddin <mdmesbah@gmail.com>
# March 26, 2023
#######################
################################
###### 250k samples
#### Aug 26, 2023
## unziped pheno file needed
#### MultiAnc
## bash /home/jupyter/ch_gwas/aou250k/run_regenie_step1_LOOCV.AoU_2023.sh gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/gwas gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/scripts_v7/regenie_step1_LOOCV.AoU_2023.sh AoU_multiAnc_2023 "gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno/QC_2023Aug11/qc_arrays*" gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno/QC_2023Aug11/aou_qc_pass.snplist gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno/QC_2023Aug11/aou_qc_pass.id gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/my_CH/phenoCH_AoU250k.fid0_iid.qcd_myeloidCA_rel_NA.26Aug2023.tsv "hasCH,hasCHvaf10,hasDTA,hasDDR,hasSF,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1" "Age_biosample_collection,sqrAge_biosample_collection,PC{1:10}" "dragen_sex_ploidy,ancestry_pred_other,site_id,Batch_CH" gcr.io/bick-aps2/regenie:v3.2.8

## AFR, AMR, EUR, EAS, SAS, MID
## while read batch; do bash /home/jupyter/ch_gwas/aou250k/run_regenie_step1_LOOCV.AoU_2023.sh gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/gwas gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/scripts_v7/regenie_step1_LOOCV.AoU_2023.sh AoU_${batch}_2023 "gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno/QC_2023Aug11/qc_arrays*" gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno/QC_2023Aug11/aou_qc_pass.snplist gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno/QC_2023Aug11/aou_qc_pass.id gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/my_CH/phenoCH_AoU250k.fid0_iid.${batch}.qcd_myeloidCA_rel_NA.26Aug2023.tsv "hasCH,hasCHvaf10,hasDTA,hasDDR,hasSF,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1" "Age_biosample_collection,sqrAge_biosample_collection,PC{1:10}" "dragen_sex_ploidy,site_id,Batch_CH" gcr.io/bick-aps2/regenie:v3.2.8; done < <(echo -e "AFR\nAMR\nEUR\nEAS\nSAS\nMID")

#########################

## Notes: 
# 1. have to remove chrx and chry from genotype file
# 2. FID==0;
## keep only autosomes
## plink2 --bfile /home/jupyter/qc_geno/arrays/arrays --extract /home/jupyter/qc_geno/arrays/aou_qc_pass.snplist --make-bed --out /home/jupyter/qc_geno/arrays/qc_arrays --threads 4  1>>log.txt 2>>err.txt  &
## test regenie: 
## regenie --step 1 --loocv --bed ~/qc_geno/arrays/arrays --extract ~/qc_geno/arrays/aou_qc_pass.snplist --keep ~/qc_geno/arrays/aou_qc_pass.id --phenoFile aou_ch_phenotype.14Mar2023.tsv.gz --phenoColList hasCH,hasDNMT3A,hasTET2,hasDTA --covarFile aou_ch_phenotype.14Mar2023.tsv.gz --covarColList Age_at_survey,Age_sqr,PC{1:10} --catCovarList Sex_at_birth,Ancestry --bt --bsize 1000 --lowmem --threads $(nproc --all) --out NULL_MODEL
# worked: regenie --step 1 --loocv --bed ~/qc_geno/arrays/qc_arrays --extract ~/qc_geno/arrays/aou_qc_pass.snplist --keep ~/qc_geno/arrays/aou_qc_pass.id --phenoFile aou_ch_phenotype.EUR.fid0_iid.14Mar2023.tsv.gz --phenoColList hasCH,hasDNMT3A,hasTET2,hasDTA --covarFile aou_ch_phenotype.EUR.fid0_iid.14Mar2023.tsv.gz --covarColList Age_at_survey,Age_sqr,PC{1:10} --catCovarList Sex_at_birth,Ancestry --bt --bsize 1000 --lowmem --threads $(nproc --all) --out NULL_MODEL
#########################


####### Aug 6, 2022
## AoU 100k CHIP
# bash /home/jupyter/ch_gwas/run_regenie_step1_LOOCV.AoU_2023.sh gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/scripts/regenie_step1_LOOCV.AoU_2023.sh AoU_multiAnc_2023 gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/qc_geno/arrays/aou_qc_pass.snplist gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/qc_geno/arrays/aou_qc_pass.id gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/aou_ch_phenotype.fid0_iid.14Mar2023.tsv.gz "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_survey,Age_sqr,PC{1:10}" "Sex_at_birth,Ancestry" gcr.io/bick-aps2/ghcr.io/rgcgithub/regenie/regenie:v3.2.4.gz

## gsutil -m cp aou_ch_phenotype.{EUR,AFR,AMR}.14Mar2023.tsv.gz gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/
## EUR
# filter monoallelic snps
# plink2 --bfile ~/qc_geno/arrays/qc_arrays --mac 100 --keep <(zcat aou_ch_phenotype.EUR.fid0_iid.14Mar2023.tsv.gz | cut -f1,2) --write-snplist --out snps_pass.eur
# gsutil -m cp snps_pass.eur* gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/qc_geno/arrays/
# bash /home/jupyter/ch_gwas/run_regenie_step1_LOOCV.AoU_2023.sh gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/scripts/regenie_step1_LOOCV.AoU_2023.sh AoU_EUR_2023 gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/qc_geno/arrays/snps_pass.eur.snplist gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/qc_geno/arrays/aou_qc_pass.id gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/aou_ch_phenotype.EUR.fid0_iid.14Mar2023.tsv.gz "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_survey,Age_sqr,PC{1:10}" "Sex_at_birth" gcr.io/bick-aps2/ghcr.io/rgcgithub/regenie/regenie:v3.2.4.gz

## AFR
# bash /home/jupyter/ch_gwas/run_regenie_step1_LOOCV.AoU_2023.sh gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/scripts/regenie_step1_LOOCV.AoU_2023.sh AoU_AFR_2023 gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/qc_geno/arrays/aou_qc_pass.snplist gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/qc_geno/arrays/aou_qc_pass.id gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/aou_ch_phenotype.AFR.fid0_iid.14Mar2023.tsv.gz "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_survey,Age_sqr,PC{1:10}" "Sex_at_birth" gcr.io/bick-aps2/ghcr.io/rgcgithub/regenie/regenie:v3.2.4.gz

## AMR
# bash /home/jupyter/ch_gwas/run_regenie_step1_LOOCV.AoU_2023.sh gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/scripts/regenie_step1_LOOCV.AoU_2023.sh AoU_AMR_2023 gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/qc_geno/arrays/aou_qc_pass.snplist gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/qc_geno/arrays/aou_qc_pass.id gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas/aou_ch_phenotype.AMR.fid0_iid.14Mar2023.tsv.gz "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_survey,Age_sqr,PC{1:10}" "Sex_at_birth" gcr.io/bick-aps2/ghcr.io/rgcgithub/regenie/regenie:v3.2.4.gz

#################################

# This script runs a UK Biobank GWAS using REGENIE on GCP using dsub. 
# First, please install dsub: https://github.com/DataBiosphere/dsub#install-dsub .

# This script is written with reference to the Broad-wide UK Biobank genetic
# data stored under gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183 
## bash run_regenie_step1_LOOCV.sh regenie_step1_LOOCV.sh JAK2_Sep8_2020 JAK2 
# BILLING_PROJECT="ukbb-analyses"
PROJECT_PATH=${1} #"gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/ch_gwas"

## Input files
REGENIE_SCRIPT=${2}

job_name=${3} # "CHIP_Jun23_2021"
	# Load Genotype Plink files 
# plink_files="gs://fc-aou-datasets-controlled/v6/microarray/plink/arrays*"
plink_files=${4} # gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno/QC_2023Aug11/qc_arrays.bim; "gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/qc_geno/arrays/qc_arrays*"	
	# SNPs to keep
	## script used: "gs://ukbb_v2/projects/muddin/ukb_v2_data/ukb_qc_step1/Jun24_2021/combine_ukb_BEDs.sh"
snps_2_keep=${5} # "gs://ukbb_v2/projects/muddin/ukb_v2_data/qced_ukb500k/qc_pass_snps.ukb_v2.list"
# snps_2_keep="gs://ukbb_v2/projects/muddin/ukb_v2_data/ukb_qc_step1/Jun24_2021/qc_pass.chr1_22.24jun2021.snplist" #"gs://ukbb_v2/projects/muddin/ukb_v2_data/ukb_qc_step1/qc_pass.chr1_22.snplist"
	# Samples to keep: one sample with >10% missing genotype "4523546"
samples_2_keep=${6} # "gs://ukbb_v2/projects/muddin/ukb_v2_data/qced_ukb500k/qc_pass_ukb_all_samples.id"
# samples_2_keep="gs://ukbb_v2/projects/muddin/ukb_v2_data/ukb_qc_step1/Jun24_2021/qc_pass.chr1_22.24jun2021.id" # "gs://ukbb_v2/projects/muddin/ukb_v2_data/ukb_qc_step1/qc_pass.chr1_22.id"

# Load Pheno file 
# pheno_file="gs://ukbb_v2/projects/muddin/chip_gwas/pheno_file_ukb_CHIP.tsv" 
# pheno_file="gs://ukbb_v2/projects/muddin/chip_gwas/traits/WhiteBrt1001.ML_CHIP.UKB_pheno.27Feb2021.tsv"
# pheno_file="gs://ukbb_v2/projects/muddin/UKB200K/pheno_covar/ukb_CHIP_CoVar_v23Jun2021.tsv"
# pheno_file="gs://ukbb_v2/projects/muddin/UKB200K/pheno_covar/UKB200k_CHIP_PhenoCovar.11Aug2021.tsv.gz"
pheno_file=${7}

# space/tab separated file with header: FID IID Y1 Y2 Y3 Y4 Y5 
#phenoCols="CHIP,VAF,DNMT3A,TET2,ASXL1,JAK2"  # comma separated list of pheno to include, e.g. Y1,Y3,Y5  
# CHIP,VAF10,DNMT3A,TET2,ASXL1,JAK2,TP53,PPM1D,SRSF2,GNB1,SF3B1
phenoColList=${8}
# Covar files
#covar_file="gs://ukbb_v2/projects/muddin/chip_gwas/covar_file_ukb_CHIP_sex01.tsv"

# numerical Covariates
covarColList=${9}

# categorical Covariates
catCovarList=${10} #"Ever_Smoked,Race"

# Bin size
bsize=1000
# num_threads= | nproc --all
# out_file_name="CHIP_GWAS_v1"

# Regenie Docker: gcr.io/ukbb-analyses/regenie:v1.0.6.0.gz
# regenie_docker=${9} # "gcr.io/ukbb-analyses/regenie:v1.0.5.8"
REGENIE_DOCKER=${11} # gcr.io/bick-aps2/ghcr.io/rgcgithub/regenie/regenie:v3.2.4.gz
# Launch step 1 and block until completion
echo "Launching REGENIE Step 1"
###
#local 
AOU_NETWORK=network
# local 
AOU_SUBNETWORK=subnetwork
DSUB_USER_NAME="$(echo "${OWNER_EMAIL}" | cut -d@ -f1)"
## e2-highcpu-32: spot $0.23744 vs $0.791488
## https://cloud.google.com/compute/vm-instance-pricing
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
	--machine-type=e2-highcpu-32 \
	--image ${REGENIE_DOCKER} \
	--skip \
	--logging ${PROJECT_PATH}/${job_name}/dsub-logs/step1 \
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



