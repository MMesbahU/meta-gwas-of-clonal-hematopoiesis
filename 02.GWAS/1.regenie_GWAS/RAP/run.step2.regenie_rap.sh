#!/bin/bash




## use Anaconda3
# reuse Python-3.9
# project-G9Yxy98JffqqBb9Y4j3PfK78:/Bulk/Imputation/Imputation from genotype (TOPmed)/
############ ### 2024: Ancestry starified GWAS
# while read ANC; do while read lines; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.sh ${ANC} "project-G9Yxy98JffqqBb9Y4j3PfK78:/Bulk/Imputation/Imputation from genotype (TOPmed)/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/pheno/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/gwas/step2/${ANC}/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/gwas/step1/${ANC}/" ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen.bgi ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.sample CH_phenoCovar.ukb200k_n_250k_${ANC}_N*.28cols.03_05_2024.tsv.gz "hasCH,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_at_recruitment,sqrtAge_at_recruitment" "Genetic_Sex,Batch,GenoBatch" 400 0.01 20 0.3 $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print "chr"$1":"$2}') mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/scripts/" 1; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv; done < <(echo -e "AMR\nSAS\nEAS\nAFR") 1>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/log_rapJobs.txt 2>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/err_rapJobs.err &

### EUR:
# while read ANC; do while read lines; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.sh ${ANC} "project-G9Yxy98JffqqBb9Y4j3PfK78:/Bulk/Imputation/Imputation from genotype (TOPmed)/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/pheno/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/gwas/step2/${ANC}/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/gwas/step1/${ANC}/" ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen.bgi ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.sample CH_phenoCovar.ukb200k_n_250k_${ANC}_N*.28cols.03_05_2024.tsv.gz "hasCH,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_at_recruitment,sqrtAge_at_recruitment" "Genetic_Sex,Batch,GenoBatch" 400 0.01 20 0.3 $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print "chr"$1":"$2}') mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/scripts/" 1; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv; done < <(echo -e "EUR") 1>>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/log_rapJobs.txt 2>>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/err_rapJobs.err &

### Male/Female: adjusted for knn
## rerun: while read SEX; do while read lines; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.sh ${SEX} "project-G9Yxy98JffqqBb9Y4j3PfK78:/Bulk/Imputation/Imputation from genotype (TOPmed)/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/pheno/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/gwas/step2/${SEX}/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/gwas/step1/${SEX}/" ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen.bgi ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.sample CH_phenoCovar.ukb200k_n_250k_${SEX}_N*.28cols.03_05_2024.tsv.gz "hasCH,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_at_recruitment,sqrtAge_at_recruitment" "knn,Batch,GenoBatch" 400 0.01 40 0.3 $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print "chr"$1":"$2}') mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/scripts/" 1; done < <(grep -w '19829556-39659111' /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv); done < <(echo -e "Male") 1>>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/log_rapJobs.txt 2>>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/err_rapJobs.err

# while read SEX; do while read lines; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.sh ${SEX} "project-G9Yxy98JffqqBb9Y4j3PfK78:/Bulk/Imputation/Imputation from genotype (TOPmed)/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/pheno/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/gwas/step2/${SEX}/" "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/gwas/step1/${SEX}/" ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen.bgi ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.sample CH_phenoCovar.ukb200k_n_250k_${SEX}_N*.28cols.03_05_2024.tsv.gz "hasCH,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_at_recruitment,sqrtAge_at_recruitment" "knn,Batch,GenoBatch" 400 0.01 40 0.3 $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print "chr"$1":"$2}') mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx "project-G9Yxy98JffqqBb9Y4j3PfK78:/ukbb_CH/2024/scripts/" 1; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv; done < <(echo -e "Male\nFemale") 1>>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/log_rapJobs.txt 2>>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/err_rapJobs.err
#### 

##########################################



## upload files: /home/unix/muddin/.local/bin/dx upload /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/step2.regenie_rap.sh --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/ch_pheno_8Mar2023/
### run
	## Test run: 
	## TERT region: bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.sh ukb200k /Bulk/Imputation/Imputation_from_genotype_TOPmed/ /ukbb_CH/ch_pheno_8Mar2023/ /ukbb_CH/ch_gwas_2023/step2/ /ukbb_CH/ch_gwas_2023/step1/ ukb21007_c5_b0_v1.bgen ukb21007_c5_b0_v1.bgen.bgi ukb21007_c5_b0_v1.sample CH_phenoCovar.ukb200k_N193342.27cols.8Mar2023.tsv.gz hasCH,hasDNMT3A,hasTET2,hasDTA Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} Genetic_Sex,Ethnic_Background,GenoBatch 400 0.01 40 0.3 5 chr5:1097331-1456860 mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.tert_test.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx

################ Multi-Ancestry ####################
    ## ukb200k
      # while read lines; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.sh ukb200k /Bulk/Imputation/Imputation_from_genotype_TOPmed/ /ukbb_CH/ch_pheno_8Mar2023/ /ukbb_CH/ch_gwas_2023/step2/ /ukbb_CH/ch_gwas_2023/step1/ ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen.bgi ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.sample CH_phenoCovar.ukb200k_N193342.27cols.8Mar2023.tsv.gz hasCH,hasDNMT3A,hasTET2,hasDTA Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} Genetic_Sex,Ethnic_Background,GenoBatch 400 0.01 40 0.3 $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print "chr"$1":"$2}') mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv

	## ukb250k
	# while read lines; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.sh ukb250k /Bulk/Imputation/Imputation_from_genotype_TOPmed/ /ukbb_CH/ch_pheno_8Mar2023/ /ukbb_CH/ch_gwas_2023/step2/ /ukbb_CH/ch_gwas_2023/step1/ ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen.bgi ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.sample CH_phenoCovar.ukb250k_N243350.27cols.8Mar2023.tsv.gz hasCH,hasDNMT3A,hasTET2,hasDTA Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} Genetic_Sex,Ethnic_Background,GenoBatch 400 0.01 40 0.3 $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print "chr"$1":"$2}') mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv
####################################################

############### Stratified by WB vs others ###########
## ukb200k
	# WB: CH_phenoCovar.WB.ukb200k_N161151.27cols.8Mar2023.tsv.gz; CH_phenoCovar.notWB.ukb200k_N32191.27cols.8Mar2023.tsv.gz 
      # while read lines; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.sh wb_ukb200k /Bulk/Imputation/Imputation_from_genotype_TOPmed/ /ukbb_CH/ch_pheno_8Mar2023/ /ukbb_CH/ch_gwas_2023/step2/ /ukbb_CH/ch_gwas_2023/step1/ ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen.bgi ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.sample CH_phenoCovar.WB.ukb200k_N161151.27cols.8Mar2023.tsv.gz hasCH,hasDNMT3A,hasTET2,hasDTA Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} Genetic_Sex,GenoBatch 400 0.01 40 0.3 $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print "chr"$1":"$2}') mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv
      # notWB
            # while read lines; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.sh notwb_ukb200k /Bulk/Imputation/Imputation_from_genotype_TOPmed/ /ukbb_CH/ch_pheno_8Mar2023/ /ukbb_CH/ch_gwas_2023/step2/ /ukbb_CH/ch_gwas_2023/step1/ ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen.bgi ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.sample CH_phenoCovar.notWB.ukb200k_N32191.27cols.8Mar2023.tsv.gz hasCH,hasDNMT3A,hasTET2,hasDTA Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} Genetic_Sex,GenoBatch 400 0.01 40 0.3 $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print "chr"$1":"$2}') mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv

      ## ukb250k CH_phenoCovar.WB.ukb250k_N204443.27cols.8Mar2023.tsv.gz CH_phenoCovar.notWB.ukb250k_N38907.27cols.8Mar2023.tsv.gz
      ## WB
      # while read lines; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.sh wb_ukb250k /Bulk/Imputation/Imputation_from_genotype_TOPmed/ /ukbb_CH/ch_pheno_8Mar2023/ /ukbb_CH/ch_gwas_2023/step2/ /ukbb_CH/ch_gwas_2023/step1/ ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen.bgi ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.sample CH_phenoCovar.WB.ukb250k_N204443.27cols.8Mar2023.tsv.gz hasCH,hasDNMT3A,hasTET2,hasDTA Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} Genetic_Sex,GenoBatch 400 0.01 40 0.3 $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print "chr"$1":"$2}') mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv

      ## notWB
      # while read lines; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step2.regenie_rap.sh notwb_ukb250k /Bulk/Imputation/Imputation_from_genotype_TOPmed/ /ukbb_CH/ch_pheno_8Mar2023/ /ukbb_CH/ch_gwas_2023/step2/ /ukbb_CH/ch_gwas_2023/step1/ ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.bgen.bgi ukb21007_c$(echo ${lines} | awk '{print $1}')_b0_v1.sample CH_phenoCovar.notWB.ukb250k_N38907.27cols.8Mar2023.tsv.gz hasCH,hasDNMT3A,hasTET2,hasDTA Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} Genetic_Sex,GenoBatch 400 0.01 40 0.3 $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print "chr"$1":"$2}') mmuddin2020/regenie:v3.2.5.2.gz step2.regenie_rap.sh mem1_ssd1_v2_x36 normal /home/unix/muddin/.local/bin/dx; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv

######################################################


### 
COHORT=${1} #"ukb200k"

GENO_PATH=${2} # /Bulk/Imputation/Imputation_from_genotype_TOPmed/
PHENO_DIR=${3} # /ukbb_CH/ch_pheno_8Mar2023/
OUT_DIR=${4}  # /ukbb_CH/ch_gwas_2023/step2/
STEP1_DIR=${5} # /ukbb_CH/ch_gwas_2023/step1/

BGEN=${6} # ukb21007_c${chr}_b0_v1.bgen
BGEN_INDEX=${7} # ukb21007_c${chr}_b0_v1.bgen.bgi
BGEN_SAMPLE=${8} # ukb21007_c${chr}_b0_v1.sample 

PHENO_FILE=${9} #"CH_phenoCovar.ukb200k_ukb250k_N421316.27cols.8Mar2023.tsv.gz"
PHENO_COL_LIST=${10} # "hasCH,hasDNMT3A,hasTET2,hasDTA"
COVAR_COL_LIST=${11} #"Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10}"
CAT_COVAR_COL_LIST=${12} # "Genetic_Sex,Ethnic_Background,GenoBatch"

GENO_BLOCK_SIZE=${13} # 400
P_THRESH=${14} # 0.01
MIN_MAC=${15} # 40
INFO=${16} # 0.3
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
## 
STEP1_LIST_FILE="NULL_MODEL.${COHORT}_pred.list"
##
${DX} run swiss-army-knife \
	-icmd="bash ${BASH_SCRIPT} ${COHORT} ${BGEN} ${BGEN_SAMPLE} ${PHENO_FILE} ${PHENO_COL_LIST} ${COVAR_COL_LIST} ${CAT_COVAR_COL_LIST} ${GENO_BLOCK_SIZE} ${P_THRESH} ${STEP1_LIST_FILE} ${MIN_MAC} ${INFO} ${CHR} ${CHR_RANGE} ${DOCKER_IMAGE} ${minCase}" \
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
	-iin="${PHENO_DIR}/${PHENO_FILE}" \
	-iin="${BASH_DIR}/${BASH_SCRIPT}" \
	--name="step2.${COHORT}" \
	--tag="step2" \
	--instance-type=${MEMORY_TYPE} \
	--priority=${PRIORITY} \
	--destination=${OUT_DIR} \
	--brief --yes

