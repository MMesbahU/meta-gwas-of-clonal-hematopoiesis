#!/bin/bash

##########################################################################
source /broad/software/scripts/useuse
use Anaconda3
## to install
# conda create -n regenie_env -c conda-forge -c bioconda regenie
## Update from v3.1g to v3.1.3g
# conda update -n regenie_env -c conda-forge -c bioconda regenie
# conda update -n base conda
source activate regenie_env
#########################################################################

#############
# regenie --step 1 --loocv  --print-prs --bed /broad/hptmp/mesbah/gwas/mgbb53k/qced.mgbb53k.autosomes --phenoFile /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.v2.tsv --covarFile /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.v2.tsv --phenoColList hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33 --covarColList PC{1:10},Age_Genotyping,sqrAge_Genotyping --catCovarList Sex,Ancestry_Self_cat,Batch_CHIP_call --bt --bsize 1000 --lowmem --threads $(nproc --all) --out /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/NULL_Model 1>log.txt 2>err.txt &
	# Multi-Ancestry GWAS
# qsub -wd /broad/hptmp/mesbah/gwas/mgbb53k/tmpdir -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 8 -binding linear:8 -N step1.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/step01.regenie.sh /broad/hptmp/mesbah/gwas/mgbb53k/qced.mgbb53k.autosomes /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv "hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,Ancestry_Self_cat,Batch_CHIP_call" 1000 /broad/hptmp/mesbah/gwas/mgbb53k/null_v2 8
	# EUR only GWAS
# qsub -wd /broad/hptmp/mesbah/gwas/mgbb53k/tmpdir -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 8 -binding linear:8 -N step1.mgbb53k.eur /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/step01.regenie.sh /broad/hptmp/mesbah/gwas/mgbb53k/qced.mgbb53k.autosomes /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k_WHITE.imp_new.noRel_sk.21Jul2022.tsv "hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,Batch_CHIP_call" 1000 /broad/hptmp/mesbah/gwas/mgbb53k/null_EUR 8
####################################################

########################### Input files
PLINK_Prefix=${1} # /broad/hptmp/mesbah/gwas/mgbb53k/qced.mgbb53k.autosomes
PHENO_FILE=${2} # /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv
PHENO_COL_LIST=${3} # hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33
COVAR_COL_LIST=${4} # PC{1:10},Age_Genotyping,sqrAge_Genotyping
CAT_COVAR_COL_LIST=${5} # Sex,Ancestry_Self_cat,Batch_CHIP_call
BIN_SIZE=${6} # 1000
OUTPUT_PATH=${7} # /broad/hptmp/mesbah/gwas/mgb53k
cpus=${8}
#########################################################################
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
mkdir -p ${OUTPUT_PATH}
## Run Regenie
regenie \
	--step 1 \
	--loocv \
	--print-prs \
	--bed ${PLINK_Prefix} \
	--phenoFile ${PHENO_FILE} \
	--phenoColList ${PHENO_COL_LIST} \
	--covarFile ${PHENO_FILE} \
	--covarColList ${COVAR_COL_LIST} \
	--catCovarList ${CAT_COVAR_COL_LIST} \
	--bt \
	--bsize ${BIN_SIZE} \
	--lowmem \
	--threads ${cpus} \
	--out ${OUTPUT_PATH}/NULL_MODEL


######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

