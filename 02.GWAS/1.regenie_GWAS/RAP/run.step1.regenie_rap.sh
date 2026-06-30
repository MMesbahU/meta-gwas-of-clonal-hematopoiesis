#!/bin/bash

################## 2024 
### qc samples and SNP: script "CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/combine_ukb_BEDs.sh" "CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/qc_geno_ukb_v2.sh"
## upload filtered sample id: script used to 
# ~/.local/bin/dx upload /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/qc_pass_ukb_*_samples.v2.id --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/2024_pheno/
# upload qc SNPs
# ~/.local/bin/dx upload /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_pass_snps.*.chr1_22.maf01_hwe1e20_geno10pct.snplist --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/2024_pheno/

##################


# use Python-3.9
## upload files: /home/unix/muddin/.local/bin/dx upload /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/*regenie_rap.sh --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/ch_pheno_8Mar2023/
## /home/unix/muddin/.local/bin/dx upload /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.{WB,notWB}.{ukb200k,ukb250k}_N*.27cols.8Mar2023.tsv.gz --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/ch_pheno_8Mar2023/

### run 
	## Multi-ancestry
	## ukb200k
# bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step1.regenie_rap.sh ukb200k ukb_v2_GT_chr1_22 qc_pass_snps.ukb_v2.list qc_pass_ukb_all_samples.id CH_phenoCovar.ukb200k_N193342.27cols.8Mar2023.tsv.gz "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10}" "Genetic_Sex,Ethnic_Background,GenoBatch" 1000 mmuddin2020/regenie:v3.2.5.2.gz step1.regenie_rap_v1.sh mem3_ssd1_v2_x16 high

	## ukb250k
# bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step1.regenie_rap.sh ukb250k ukb_v2_GT_chr1_22 qc_pass_snps.ukb_v2.list qc_pass_ukb_all_samples.id CH_phenoCovar.ukb250k_N243350.27cols.8Mar2023.tsv.gz "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10}" "Genetic_Sex,Ethnic_Background,GenoBatch" 1000 mmuddin2020/regenie:v3.2.5.2.gz step1.regenie_rap_v2.sh mem1_ssd1_v2_x36 high
#############
	## Stratified
	        ## ukb200k
# bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step1.regenie_rap.sh wb_ukb200k ukb_v2_GT_chr1_22 qc_pass_snps.ukb_v2.list qc_pass_ukb_all_samples.id CH_phenoCovar.WB.ukb200k_N161151.27cols.8Mar2023.tsv.gz "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10}" "Genetic_Sex,GenoBatch" 1000 mmuddin2020/regenie:v3.2.5.2.gz step1.regenie_rap_v2.sh mem3_ssd1_v2_x16 high

# bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step1.regenie_rap.sh notwb_ukb200k ukb_v2_GT_chr1_22 qc_pass_snps.ukb_v2.list qc_pass_ukb_all_samples.id CH_phenoCovar.notWB.ukb200k_N32191.27cols.8Mar2023.tsv.gz "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10}" "Genetic_Sex,GenoBatch" 1000 mmuddin2020/regenie:v3.2.5.2.gz step1.regenie_rap_v2.sh mem3_ssd1_v2_x16 high

	       ## ukb250k
# bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step1.regenie_rap.sh wb_ukb250k ukb_v2_GT_chr1_22 qc_pass_snps.ukb_v2.list qc_pass_ukb_all_samples.id CH_phenoCovar.WB.ukb250k_N204443.27cols.8Mar2023.tsv.gz "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10}" "Genetic_Sex,GenoBatch" 1000 mmuddin2020/regenie:v3.2.5.2.gz step1.regenie_rap_v2.sh mem1_ssd1_v2_x36 high

# bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step1.regenie_rap.sh notwb_ukb250k ukb_v2_GT_chr1_22 qc_pass_snps.ukb_v2.list qc_pass_ukb_all_samples.id CH_phenoCovar.notWB.ukb250k_N38907.27cols.8Mar2023.tsv.gz "hasCH,hasDNMT3A,hasTET2,hasDTA" "Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10}" "Genetic_Sex,GenoBatch" 1000 mmuddin2020/regenie:v3.2.5.2.gz step1.regenie_rap_v2.sh mem1_ssd1_v2_x36 high


###
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
bash_script=${11} # step1.regenie_rap_v1.sh
memory_type=${12} # mem1_ssd1_v2_x36 70/640gb/0.8928/hr  mem2_ssd1_v2_x32 125/1116/.9056/h
priority=${13} # normal/ high/low
# 
#####
geno_dir="/ukbb_CH/ukb_qc_geno/"
pheno_dir="/ukbb_CH/ch_pheno_8Mar2023/"
out_dir="/ukbb_CH/ch_gwas_2023/step1/"
dx="/home/unix/muddin/.local/bin/dx"
##
${dx} run swiss-army-knife \
	-icmd="bash ${bash_script} ${cohort} ${PLINK_Prefix} ${KEEP_SNPs} ${KEEP_SAMPLE} ${PHENO_FILE} ${PHENO_COL_LIST} ${COVAR_COL_LIST} ${CAT_COVAR_COL_LIST} ${BIN_SIZE} ${docker_image}" \
	-iin="${geno_dir}/ukb_v2_GT_chr1_22.bed" \
	-iin="${geno_dir}/ukb_v2_GT_chr1_22.bim" \
	-iin="${geno_dir}/ukb_v2_GT_chr1_22.fam" \
	-iin="${geno_dir}/qc_pass_snps.ukb_v2.list" \
	-iin="${geno_dir}/qc_pass_ukb_all_samples.id" \
	-iin="${pheno_dir}/${PHENO_FILE}" \
	-iin="${pheno_dir}/${bash_script}" \
	--name="step1.${cohort}" \
	--tag="step1" \
	--instance-type=${memory_type} \
	--priority=${priority} \
	--destination="${out_dir}" \
	--brief --yes

