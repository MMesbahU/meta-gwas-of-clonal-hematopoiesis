#!/bin/bash




## use Anaconda3
# reuse Python-3.9
# project-G9Yxy98JffqqBb9Y4j3PfK78:/Bulk/Imputation/Imputation from genotype (TOPmed)/

###################################
#### 2024-09-26
## upload scripts: ~/.local/bin/dx upload /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/step2.regenie_rap.MultiAnc.sh  --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/2024/scripts/
## MultiAnc CHIP types GWAS
# while read ukb_batch; do while read ANC; do echo ${ANC}.${ukb_batch}; while read lines; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.MultiAnc.sh ${ANC} "project-G9Yxy98JffqqBb9Y4j3PfK78:/Bulk/Imputation/Imputation from genotype (TOPmed)/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/pheno/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/gwas/step2/${ANC}/${ukb_batch}/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/gwas/step1/${ANC}/${ukb_batch}/" ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen.bgi ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.sample CH_phenoCovar.${ukb_batch}.28cols.03_05_2024.tsv.gz "hasCH,hasCHvaf05,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_at_recruitment,sqrtAge_at_recruitment" "Genetic_Sex,knn,GenoBatch" 400 0.01 20 0.3 $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print "chr"$1":"$2}') mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.MultiAnc.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/scripts/" 1 ${ukb_batch}; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv; done < <(echo -e "MultiAnc"); done < <(echo -e "ukb200k_N193342\nukb250k_N243350") 1>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/MultiAnc.rapJobs.20240926.log 2>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/MultiAnc.rapJobs.20240926.err &
######################################################
###
COHORT=${1} #"ukb200k"

GENO_PATH=${2} # /Bulk/Imputation/Imputation_from_genotype_TOPmed/
PHENO_DIR=${3} # /ukbb_CH/2024/pheno/
OUT_DIR=${4}  # /ukbb_CH/2024/gwas/step2/
STEP1_DIR=${5} # /ukbb_CH/2024/gwas/step1/

BGEN=${6} # ukb21007_c${chr}_b0_v1.bgen
BGEN_INDEX=${7} # ukb21007_c${chr}_b0_v1.bgen.bgi
BGEN_SAMPLE=${8} # ukb21007_c${chr}_b0_v1.sample 

PHENO_FILE=${9} #"CH_phenoCovar.${ukb_batch}.28cols.03_05_2024.tsv.gz"
PHENO_COL_LIST=${10} # "hasCH,hasCHvaf05,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR"
COVAR_COL_LIST=${11} #"Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10}"
CAT_COVAR_COL_LIST=${12} # "Genetic_Sex,knn,GenoBatch"

GENO_BLOCK_SIZE=${13} # 400
P_THRESH=${14} # 0.01
MIN_MAC=${15} # MAC>=20
INFO=${16} # INFO score >=0.3
CHR=${17} # 1, 2, ,3 
CHR_RANGE=${18} # chr1:1-10000

DOCKER_IMAGE=${19} # mmuddin2020/regenie:v3.2.5.2.gz
BASH_SCRIPT=${20} # step2.regenie_rap.sh
MEMORY_TYPE=${21} # mem1_ssd1_v2_x36 70/640gb/0.8928/hr  mem2_ssd1_v2_x32 125/1116/.9056/h
PRIORITY=${22} # normal/ high/low
# 
#####
DX=${23} # "/home/unix/muddin/.local/bin/dx"

BASH_DIR=${24} # /ukbb_CH/2024/scripts

minCase=${25}

BATCH=${26} # ukb200k; ukb250k
## 
STEP1_LIST_FILE="NULL_MODEL.${COHORT}_pred.list"
##
${DX} run swiss-army-knife \
	-icmd="bash ${BASH_SCRIPT} ${COHORT} ${BGEN} ${BGEN_SAMPLE} ${PHENO_FILE} ${PHENO_COL_LIST} ${COVAR_COL_LIST} ${CAT_COVAR_COL_LIST} ${GENO_BLOCK_SIZE} ${P_THRESH} ${STEP1_LIST_FILE} ${MIN_MAC} ${INFO} ${CHR} ${CHR_RANGE} ${DOCKER_IMAGE} ${minCase} ${BATCH}" \
	-iin="${GENO_PATH}/${BGEN}" \
	-iin="${GENO_PATH}/${BGEN_INDEX}" \
	-iin="${GENO_PATH}/${BGEN_SAMPLE}" \
	-iin="${STEP1_DIR}/${STEP1_LIST_FILE}" \
	-iin="${STEP1_DIR}/NULL_MODEL.${COHORT}_1.loco" \
	-iin="${STEP1_DIR}/NULL_MODEL.${COHORT}_2.loco" \
	-iin="${STEP1_DIR}/NULL_MODEL.${COHORT}_3.loco" \
	-iin="${STEP1_DIR}/NULL_MODEL.${COHORT}_4.loco" \
	-iin="${STEP1_DIR}/NULL_MODEL.${COHORT}_5.loco" \
	-iin="${STEP1_DIR}/NULL_MODEL.${COHORT}_6.loco" \
	-iin="${STEP1_DIR}/NULL_MODEL.${COHORT}_7.loco" \
	-iin="${STEP1_DIR}/NULL_MODEL.${COHORT}_8.loco" \
	-iin="${PHENO_DIR}/${PHENO_FILE}" \
	-iin="${BASH_DIR}/${BASH_SCRIPT}" \
	--name="step2.${COHORT}.${BATCH}" \
	--tag="step2" \
	--instance-type=${MEMORY_TYPE} \
	--priority=${PRIORITY} \
	--destination=${OUT_DIR} \
	--brief --yes

