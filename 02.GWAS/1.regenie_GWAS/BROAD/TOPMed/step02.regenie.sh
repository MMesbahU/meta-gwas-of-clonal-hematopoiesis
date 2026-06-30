#!/bin/bash


##### Note:
#### --ref-first	FLAG	
## Optional	Specify to use the first allele as the reference allele for BGEN or PLINK bed/bim/fam file input 
# [default is to use the last allele as the reference]
##

######
### TOPMed bgen has ref-allele first.
## have to use "--ref-first "
## for the GWAS run in April-May 2024; I did not use "--ref-first" for TOPMed data.
## Effect estimates and EAF are for Reference Allele in the resulting TOPMed GWAS summary 
########

### !!! Note
#### REF allele is the effect allele in MGBB53k GWAS with --ref-first flag
#### To make ALT allele as EA, omit --ref-first flag
### !!!

##########################################################################
source /broad/software/scripts/useuse
use Anaconda3
## to install
# conda create -n regenie_env -c conda-forge -c bioconda regenie
## Update from v3.1g to v3.1.3g
# conda update -n regenie_env -c conda-forge -c bioconda regenie
# conda update -n base conda
# source activate regenie_env
## Initialize Micromamba
eval "$(/home/unix/muddin/bin/micromamba shell hook --shell bash)"
# Activate the Environment
micromamba activate regenie_env
#########################################################################

######### April 2024

#Ancestry:
# while read ANC; do qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N step2.${ANC}.topmed74k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step02.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.${ANC}.04_23_2024.v02.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:20},AgeAtBloodDraw,Sqrd_AgeAtBloodDraw" "Genetic_Sex,STUDY,Sequencing_Center,TOPMed_Phase" 400 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step2 ${ANC}.TOPMed74k 0.01 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step1/NULL_MODEL_pred.list 1 0.3 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv 4 /medpop/esp2/projects/topmed/freeze.12c/bgen/freeze.12c pass_only 30 1; done < <(echo -e "EAS")

# while read ANC; do qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N step2.${ANC}.topmed74k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step02.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.${ANC}.04_23_2024.v02.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasSF,hasDDR" "PC{1:20},AgeAtBloodDraw,Sqrd_AgeAtBloodDraw" "Genetic_Sex,STUDY,Sequencing_Center,TOPMed_Phase" 400 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step2 ${ANC}.TOPMed74k 0.01 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step1/NULL_MODEL_pred.list 1 0.3 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv 4 /medpop/esp2/projects/topmed/freeze.12c/bgen/freeze.12c pass_only 30 1; done < <(echo -e "SAS")

## 1-10 | 11-20 | 21-90; 
# 91-120 | 121-159
# while read ANC; do qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N step2.${ANC}.topmed74k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step02.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.${ANC}.04_23_2024.v02.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:20},AgeAtBloodDraw,Sqrd_AgeAtBloodDraw" "Genetic_Sex,STUDY,Sequencing_Center,TOPMed_Phase" 400 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step2 ${ANC}.TOPMed74k 0.01 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step1/NULL_MODEL_pred.list 20 0.3 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv 4 /medpop/esp2/projects/topmed/freeze.12c/bgen/freeze.12c pass_only 30 5; done < <(echo -e "AFR\nAMR\nEUR")

## Sex: 
# while read ANC; do qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N step2.${ANC}.topmed74k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step02.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.${ANC}.04_23_2024.v02.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:20},AgeAtBloodDraw,Sqrd_AgeAtBloodDraw" "GenANC,STUDY,Sequencing_Center,TOPMed_Phase" 400 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step2 ${ANC}.TOPMed74k 0.01 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step1/NULL_MODEL_pred.list 20 0.3 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv 4 /medpop/esp2/projects/topmed/freeze.12c/bgen/freeze.12c pass_only 30 5; done < <(echo -e "Female\nMale")

# Multi
# while read ANC; do qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N step2.${ANC}.topmed74k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step02.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.${ANC}.04_23_2024.v02.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:20},AgeAtBloodDraw,Sqrd_AgeAtBloodDraw" "Genetic_Sex,GenANC,STUDY,Sequencing_Center,TOPMed_Phase" 400 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step2 ${ANC}.TOPMed74k 0.01 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step1/NULL_MODEL_pred.list 20 0.3 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv 4 /medpop/esp2/projects/topmed/freeze.12c/bgen/freeze.12c pass_only 30 5; done < <(echo -e "MultiANC")


####################################### END #################################################################

########################### Input files
# BGEN_PREFIX=${1} # /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged.chr{1.bgen}
# BGEN_SAMPLE_PREFIX=${2} # /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged.chr{1.sample}

PHENO_FILE=${1} # /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv

PHENO_COL_LIST=${2} # hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33

COVAR_COL_LIST=${3} # PC{1:10},Age_Genotyping,sqrAge_Genotyping

CAT_COVAR_COL_LIST=${4} # Sex,Ancestry_Self_cat,Batch_CHIP_call

GENO_BLOCK_SIZE=${5} # 400

OUTPUT_PATH=${6} # /broad/hptmp/mesbah/gwas/mgb53k

COHORT=${7}

P_THRESH=${8}

STEP1_LIST_FILE=${9} # /broad/hptmp/mesbah/gwas/mgbb53k/NULL_MODEL_pred.list

MIN_MAC=${10}

INFO=${11}

interval_file=${12}

cpus=${13}

bgen_prefix=${14}

bgen_suffix=${15}

max_CatLevels=${16}

minCase=${17}
#
chr=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${interval_file} )

CHR_RANGE=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1":"$2}' ${interval_file} ) # 1:1-62312655

BGEN=${bgen_prefix}.chr${chr}.${bgen_suffix}.bgen # /broad/hptmp/mesbah/dataset/ch_gwas/topmed/bgen/freeze.12c.chr13.pass_only.bgen

BGEN_SAMPLE=${bgen_prefix}.chr${chr}.${bgen_suffix}.sample # /broad/hptmp/mesbah/dataset/ch_gwas/topmed/bgen/freeze.12c.chr13.pass_only.sample

###
mkdir -p ${OUTPUT_PATH}
###
OUT_PREFIX=${OUTPUT_PATH}/${COHORT}.chr${CHR_RANGE}
#########################################################################
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
echo -e "##########################\n"
echo -e "Job interval: chr${chr}\t${CHR_RANGE}\n"
echo -e "##########################\n"
########
# For TOPMed bgen, have to use "--ref-first" to get effect estimates for Alternative Allele

## For MGBB, have to ommit "--ref-first", coz of the bgen allele format in MGBB topmed imputed data
## Run Regenie
regenie \
	--step 2 \
	--minCaseCount ${minCase} \
	--bgen ${BGEN} \
	--sample ${BGEN_SAMPLE} \
	--bt \
	--htp ${COHORT} \
	--phenoFile ${PHENO_FILE} \
	--phenoColList ${PHENO_COL_LIST} \
	--covarFile ${PHENO_FILE} \
	--covarColList ${COVAR_COL_LIST} \
	--catCovarList ${CAT_COVAR_COL_LIST} \
	--bsize ${GENO_BLOCK_SIZE} \
	--firth \
	--approx \
	--gz \
	--write-samples \
	--print-pheno \
	--no-split \
	--pThresh ${P_THRESH} \
	--pred ${STEP1_LIST_FILE} \
	--minMAC ${MIN_MAC} \
	--minINFO ${INFO} \
	--maxCatLevels ${max_CatLevels} \
	--range ${CHR_RANGE} \
	--threads ${cpus} \
	--out ${OUT_PREFIX}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

