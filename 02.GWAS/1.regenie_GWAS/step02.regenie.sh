#!/bin/bash

### !!! Note
#### REF allele is the effect allele in MGBB53k GWAS with --ref-first flag
#### To make ALT allele as EA, omit --ref-first flag
### !!!

##########################################################################
source /broad/software/scripts/useuse
use Anaconda
## to install
# conda create -n regenie_env -c conda-forge -c bioconda regenie
## Update from v3.1g to v3.1.3g
# conda update -n regenie_env -c conda-forge -c bioconda regenie
# conda update -n base conda
source activate regenie_env
#########################################################################

#############
## while read lines; do echo $(echo ${lines} | awk '{print $1}') $(echo ${lines} | awk '{print $2}');done < /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv
## 
# multi-ancestry GWAS:
# qsub -wd /broad/hptmp/mesbah/gwas/mgbb53k/step2 -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 10 -binding linear:10 -N step2.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/step02.regenie.sh /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv "hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,Ancestry_Self_cat,Batch_CHIP_call" 400 /broad/hptmp/mesbah/gwas/mgbb53k/step2 MGBB53k 0.01 /broad/hptmp/mesbah/gwas/mgbb53k/NULL_MODEL_pred.list 20 0.3 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv 10

## EUR only GWAS
# qsub -wd /broad/hptmp/mesbah/gwas/mgbb53k/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 8 -binding linear:8 -N EUR.step2.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/step02.regenie.sh /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k_WHITE.imp_new.noRel_sk.21Jul2022.tsv "hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,Batch_CHIP_call" 400 /broad/hptmp/mesbah/gwas/mgbb53k/step2_EUR MGBB53k_EUR 0.01 /broad/hptmp/mesbah/gwas/mgbb53k/null_EUR/NULL_MODEL_pred.list 20 0.3 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv 8

####################################################

########################### Input files
BGEN_PREFIX=${1} # /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged.chr{1.bgen}
BGEN_SAMPLE_PREFIX=${2} # /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged.chr{1.sample}
PHENO_FILE=${3} # /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv
PHENO_COL_LIST=${4} # hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33
COVAR_COL_LIST=${5} # PC{1:10},Age_Genotyping,sqrAge_Genotyping
CAT_COVAR_COL_LIST=${6} # Sex,Ancestry_Self_cat,Batch_CHIP_call
GENO_BLOCK_SIZE=${7} # 400
OUTPUT_PATH=${8} # /broad/hptmp/mesbah/gwas/mgb53k
COHORT=${9}
P_THRESH=${10}
STEP1_LIST_FILE=${11} # /broad/hptmp/mesbah/gwas/mgbb53k/NULL_MODEL_pred.list
MIN_MAC=${12}
INFO=${13}
interval_file=${14}
cpus=${15}
#
chr=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${interval_file} )
CHR_RANGE=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1":"$2}' ${interval_file} ) # 1:1-62312655
BGEN=${BGEN_PREFIX}.chr${chr}.bgen
BGEN_SAMPLE=${BGEN_SAMPLE_PREFIX}.chr${chr}.sample
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
#
## Run Regenie
regenie \
	--step 2 \
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
	--ref-first \
	--approx \
	--gz \
	--write-samples \
	--print-pheno \
	--no-split \
	--pThresh ${P_THRESH} \
	--pred ${STEP1_LIST_FILE} \
	--minMAC ${MIN_MAC} \
	--minINFO ${INFO} \
	--range ${CHR_RANGE} \
	--threads ${cpus} \
	--out ${OUT_PREFIX}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

