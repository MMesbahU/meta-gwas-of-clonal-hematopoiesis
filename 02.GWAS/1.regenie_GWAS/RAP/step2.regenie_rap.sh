#!/bin/bash

## upload files: /home/unix/muddin/.local/bin/dx upload /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/step2.regenie_rap*.sh --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/ch_pheno_8Mar2023/
################
COHORT=${1}

BGEN=${2}
BGEN_SAMPLE=${3}
PHENO_FILE=${4}

PHENO_COL_LIST=${5}
COVAR_COL_LIST=${6}
CAT_COVAR_COL_LIST=${7}

GENO_BLOCK_SIZE=${8}
P_THRESH=${9}
STEP1_LIST_FILE=${10}
MIN_MAC=${11}

INFO=${12}
CHR=${13} # 1
CHR_RANGE=${14} # chr1:1-10000

DOCKER_IMAGE=${15} # mmuddin2020/regenie:v3.2.5.2.gz

minCase=${16}
###

# ls -lhrt /tmp

# ls -lhrt /mnt/project

# find . | grep -E 'ukb21007_c|CH_phenoCovar'

# docker pull ${docker_image}
docker run -w /tmp -v /home/dnanexus/out/out:/tmp -v /mnt/project:/mnt/project ${DOCKER_IMAGE} regenie \
	--step 2 \
	--bgen ${BGEN} \
	--sample ${BGEN_SAMPLE} \
	--bt \
	--minCaseCount ${minCase} \
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
	--out chr${CHR}_${COHORT}.${CHR}_$(echo ${CHR_RANGE} | awk -F ':' '{print $2}')

