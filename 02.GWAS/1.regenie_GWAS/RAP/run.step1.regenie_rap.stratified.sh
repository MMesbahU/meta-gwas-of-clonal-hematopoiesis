#!/bin/bash


###
# use Python-3.9
###
### To load dx; we need python3
# use Anaconda3
## /home/unix/muddin/.local/bin/dx login 
# provide user name
# provide pass
# Use dx login --timeout to control the expiration date, 
# or dx logout to end this session
#################

################################################## 2024-09-26 
# rerun MultiAnc UKB GWAS
# qc_pass_ukb_MultiAnc_samples.v2.id qc_pass_snps.MultiAnc.chr1_22.maf01_hwe1e20_geno10pct.snplist CH_phenoCovar.ukb200k_N193342.28cols.03_05_2024.tsv.gz CH_phenoCovar.ukb250k_N243350.28cols.03_05_2024.tsv.gz

## Minimum case count of >=1

## ukb200k_N193342; ukb250k_N243350
# while read ukb_batch; do while read ANC; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step1.regenie_rap.stratified.sh ${ANC} ukb_v2_GT_chr1_22 qc_pass_snps.${ANC}.chr1_22.maf01_hwe1e20_geno10pct.snplist qc_pass_ukb_${ANC}_samples.v2.id CH_phenoCovar.${ukb_batch}.28cols.03_05_2024.tsv.gz "hasCH,hasCHvaf05,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_at_recruitment,sqrtAge_at_recruitment" "Genetic_Sex,knn,GenoBatch" 1000 mmuddin2020/regenie:v3.2.5.2.gz step1.regenie_rap_v1.sh mem3_ssd1_v2_x16 high "/ukbb_CH/ukb_qc_geno" "/ukbb_CH/2024/pheno" "/ukbb_CH/2024/gwas/step1/${ANC}" "/ukbb_CH/2024/pheno" "/ukbb_CH/2024/pheno" "/ukbb_CH/2024/scripts" "/home/unix/muddin/.local/bin/dx" 1; done < <(echo -e "MultiAnc"); done < <(echo -e "ukb200k_N193342\nukb250k_N243350")

####################################################


################## 2024 
### qc samples and SNP: script "CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/combine_ukb_BEDs.sh" "CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/qc_geno_ukb_v2.sh"
## upload filtered sample id: script used to 
# ~/.local/bin/dx upload /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/qc_pass_ukb_*_samples.v2.id --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/2024/pheno/
# upload qc SNPs
# ~/.local/bin/dx upload /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_pass_snps.*.chr1_22.maf01_hwe1e20_geno10pct.snplist --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/2024/pheno/
## bash scripts: 
# ~/.local/bin/dx upload /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/step1.regenie_rap_v1.sh  --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/2024/scripts/

## Upload Phenotype file:
## ~/.local/bin/dx upload /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.*.03_05_2024.tsv.gz --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/2024/pheno/
# ~/.local/bin/dx mkdir -p /ukbb_CH/2024/gwas/{step1,step2}
## Step 1 files from non EUR samples run in Braod HPC
# ~/.local/bin/dx upload -r /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step1/* --path Natarajan_Lab_WGS_WES_Project:/ukbb_CH/2024/gwas/step1/
#############

#### Run Step1 for all ancestries:
# while read ANC; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step1.regenie_rap.stratified.sh ${ANC} ukb_v2_GT_chr1_22 qc_pass_snps.${ANC}.chr1_22.maf01_hwe1e20_geno10pct.snplist qc_pass_ukb_${ANC}_samples.v2.id CH_phenoCovar.ukb200k_n_250k_${ANC}_N*.28cols.03_05_2024.tsv.gz "hasCH,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_at_recruitment,sqrtAge_at_recruitment" "Genetic_Sex,Batch,GenoBatch" 1000 mmuddin2020/regenie:v3.2.5.2.gz step1.regenie_rap_v1.sh mem3_ssd1_v2_x16 high "/ukbb_CH/ukb_qc_geno" "/ukbb_CH/2024/pheno" "/ukbb_CH/2024/gwas/step1/${ANC}" "/ukbb_CH/2024/pheno" "/ukbb_CH/2024/pheno" "/ukbb_CH/2024/scripts" "/home/unix/muddin/.local/bin/dx" 1; done < <(echo -e "AFR\nAMR\nSAS\nEAS\nEUR")
## Male vs Female
## qc samples and snp in "MultiAnc"
## have to adjust for genetic ancestry "knn
# ANC="MultiAnc"; while read SEX; do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/run.step1.regenie_rap.stratified.sh ${SEX} ukb_v2_GT_chr1_22 qc_pass_snps.${ANC}.chr1_22.maf01_hwe1e20_geno10pct.snplist qc_pass_ukb_${ANC}_samples.v2.id CH_phenoCovar.ukb200k_n_250k_${SEX}_N*.28cols.03_05_2024.tsv.gz "hasCH,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_at_recruitment,sqrtAge_at_recruitment" "knn,Batch,GenoBatch" 1000 mmuddin2020/regenie:v3.2.5.2.gz step1.regenie_rap_v1.sh mem3_ssd1_v2_x16 high "/ukbb_CH/ukb_qc_geno" "/ukbb_CH/2024/pheno" "/ukbb_CH/2024/gwas/step1/${SEX}" "/ukbb_CH/2024/pheno" "/ukbb_CH/2024/pheno" "/ukbb_CH/2024/scripts" "/home/unix/muddin/.local/bin/dx"; done < <(echo -e "Male\nFemale")

##################

################### 2023
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
# /ukbb_CH/ukb_qc_geno /ukbb_CH/2024/pheno /ukbb_CH/2024/gwas/step1/${AFR,} /ukbb_CH/2024/pheno /ukbb_CH/2024/pheno /ukbb_CH/2024/pheno
#####
geno_dir=${14} # "/ukbb_CH/ukb_qc_geno/"
pheno_dir=${15} #"/ukbb_CH/ch_pheno_8Mar2023/"
out_dir=${16} #"/ukbb_CH/ch_gwas_2023/step1/"
snp2keep_dir=${17}
sam2keep_dir=${18}
bash_script_dir=${19}
dx=${20} # "/home/unix/muddin/.local/bin/dx"

minCase=${21} ## minimum case counts
##
${dx} run swiss-army-knife \
	-icmd="bash ${bash_script} ${cohort} ${PLINK_Prefix} ${KEEP_SNPs} ${KEEP_SAMPLE} ${PHENO_FILE} ${PHENO_COL_LIST} ${COVAR_COL_LIST} ${CAT_COVAR_COL_LIST} ${BIN_SIZE} ${docker_image} ${minCase}" \
	-iin="${geno_dir}/${PLINK_Prefix}.bed" \
	-iin="${geno_dir}/${PLINK_Prefix}.bim" \
	-iin="${geno_dir}/${PLINK_Prefix}.fam" \
	-iin="${snp2keep_dir}/${KEEP_SNPs}" \
	-iin="${sam2keep_dir}/${KEEP_SAMPLE}" \
	-iin="${pheno_dir}/${PHENO_FILE}" \
	-iin="${bash_script_dir}/${bash_script}" \
	--name="step1.${cohort}" \
	--tag="step1" \
	--instance-type=${memory_type} \
	--priority=${priority} \
	--destination="${out_dir}" \
	--brief --yes

