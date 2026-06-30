#!/bin/bash

## modify bgen sample files by replacing fid by "0" to batch with array genotypes, null step, and phenotype file.
## awk '(NR<3){print $0}(NR>2){print 0,$2,$3}' aou-alpha3-filtered-chr22.sample > aou-alpha3-filtered-chr22_fid0_iid.sample
# test run worked: regenie --step 2 --bgen aou-alpha3-filtered-chr22.bgen --sample aou-alpha3-filtered-chr22_fid0_iid.sample --bt --htp AFR --phenoFile aou_ch_phenotype.AFR.fid0_iid.14Mar2023.tsv.gz --phenoColList hasCH,hasDNMT3A,hasTET2,hasDTA --covarFile aou_ch_phenotype.AFR.fid0_iid.14Mar2023.tsv.gz --covarColList Age_at_survey,Age_sqr,PC{1:10} --catCovarList Sex_at_birth --bt --bsize 400 --firth --ref-first --approx --gz --write-samples --print-pheno --pThresh 0.01 --pred NULL_MODEL_pred_nopath.list --minMAC 40 --minINFO 0.3 --threads $(nproc --all) --range chr22:28619558-29203497 --out test.chr22.afr
##

## Prepare loco_list file
#loco_path=$(dirname ${STEP1_LOCO_FILE})
#loco_file_name=$(basename ${STEP1_LOCO_FILE})
#Loco_File=${loco_path}/${loco_file_name}
#echo "${PHENO_COL} ${Loco_File}" >> Step1_pred.list
#STEP1_LIST_FILE="Step1_pred.list"
####### 
## modify bgen sample files by replacing fid by "0" to batch with array genotypes, null step, and phenotype file.

awk '(NR<3){print $0}(NR>2){print 0,$2,$3}' ${BGEN_SAMPLE} > acaf_threshold.chr${CHR}_fid0_iid.sample

SAMPLE_FILE="acaf_threshold.chr${CHR}_fid0_iid.sample"

#######
echo -e "Starting REGENIE GWAS Run\n"
##
regenie \
	--step 2 \
	--bgen ${BGEN_FILE} \
	--sample ${SAMPLE_FILE} \
	--bt \
	--htp ${COHORT} \
	--phenoFile ${PHENO_FILE} \
	--phenoColList ${PHENO_COL_LIST} \
	--covarFile ${PHENO_FILE} \
	--covarColList ${COVAR_COL_LIST} \
	--catCovarList ${CAT_COVAR_COL_LIST} \
	--bsize ${GENO_BLOCK_SIZE} \
	--firth \
	--ref-first \
	--approx \
	--gz \
	--write-samples \
	--print-pheno \
	--pThresh ${P_THRESH} \
	--pred ${STEP1_LIST_FILE} \
	--minMAC ${MIN_MAC} \
	--minINFO ${INFO} \
	--range ${CHR_RANGE} \
	--threads $(nproc --all) \
	--out ${GWAS_OUT_PATH}/chr${CHR}_${COHORT}.${CHR}_$(echo ${CHR_RANGE} | awk -F ':' '{print $2}')


