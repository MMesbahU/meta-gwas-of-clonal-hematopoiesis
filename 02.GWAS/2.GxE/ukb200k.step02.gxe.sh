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
# WB only 450k GWAS
# qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg19_chrom_interval.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 6 -binding linear:6 -N eur_ukb200k_gxe /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/HF_GxCH/scripts/ukb200k.step02.gxe.sh /broad/ukbb/imputed_v3/ukb_imp_chr /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr3_v3_s487395.sample /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/eur_ukb450k.ch_hf.tsv has_disease Age,Age2,CHIP,PC{1:10} SEX,CHIP_Batch,GenoBatch 1000 0.01 CHIP /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg19_chrom_interval.tsv 6 /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/full_gxe HFxCH_ukb_eur
####################################################
########################### Input files
BGEN_PREFIX=${1} # /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen
BGEN_SAMPLE_PREFIX=${2} # /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr3_v3_s487395.sample 
PHENO_FILE=${3} # /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv
PHENO_COL=${4} # hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33
COVAR_COL_LIST=${5} # PC{1:10},Age_Genotyping,sqrAge_Genotyping
CAT_COVAR_COL_LIST=${6} # Sex,Ancestry_Self_cat,Batch_CHIP_call
GENO_BLOCK_SIZE=${7} # 400
P_THRESH=${8}
interaction_coVAR=${9}
interval_file=${10}
cpus=${11}
OUTPUT_PATH=${12}
COHORT=${13} 
#
chr=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${interval_file} )
CHR_RANGE=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1":"$2}' ${interval_file} ) # 1:1-62312655
BGEN=${BGEN_PREFIX}${chr}_v3.bgen # /broad/ukbb/imputed_v3/ukb_imp_chr1_v3.bgen
BGEN_SAMPLE=${BGEN_SAMPLE_PREFIX}
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

## Run Regenie
regenie \
	--step 2 \
	--ignore-pred \
	--bgen ${BGEN} \
	--sample ${BGEN_SAMPLE} \
	--bt \
	--htp ${COHORT} \
	--phenoFile ${PHENO_FILE} \
	--phenoCol ${PHENO_COL} \
	--covarFile ${PHENO_FILE} \
	--covarColList ${COVAR_COL_LIST} \
	--catCovarList ${CAT_COVAR_COL_LIST} \
	--interaction ${interaction_coVAR} \
	--bsize ${GENO_BLOCK_SIZE} \
	--firth \
	--ref-first \
	--approx \
	--gz \
	--write-samples \
	--print-pheno \
	--pThresh ${P_THRESH} \
	--range ${CHR_RANGE} \
	--threads ${cpus} \
	--minMAC 100 \
	--no-condtl \
	--out ${OUT_PREFIX}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

