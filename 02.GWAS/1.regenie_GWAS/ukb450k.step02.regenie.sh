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
# multi-ancestry 450k: 
# cat /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/450k_null/step1/NULL_MODEL_pred.list | sed 's:/mnt/data/output/gs/ukbb_v2/projects/muddin/chip_gwas/UKB450k_6Aug2022/step1/:/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/450k_null/step1/:g' > /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/450k_null/step1/NULL_MODEL_pred_hpc_modified.list
# qsub -wd /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg19_chrom_interval.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 10 -binding linear:10 -N step2.ukb450k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/ukb450k.step02.regenie.sh /broad/ukbb/imputed_v3/ukb_imp_chr /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr3_v3_s487395.sample /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipCallnomiss_CV_AB.ukb450k_all.tsv hasCHIP,hasDNMT3A,hasTET2,hasASXL1 Age,Age2,PC{1:10} Gender,Ethnicity,GenoBatch 400 /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/additive UKB450k 0.01 /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/450k_null/step1/NULL_MODEL_pred_hpc_modified.list 20 0.30 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg19_chrom_interval.tsv 10

# multi-ancestry 250k: 
# cat /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/250k_null/step1/NULL_MODEL_pred.list | sed 's:/mnt/data/output/gs/ukbb_v2/projects/muddin/chip_gwas/UKB250k_6Aug2022/step1/:/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/250k_null/step1/:g' > /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/250k_null/step1/NULL_MODEL_pred_hpc_modified.list
# qsub -wd /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg19_chrom_interval.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 10 -binding linear:10 -N step2.ukb250k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/ukb450k.step02.regenie.sh /broad/ukbb/imputed_v3/ukb_imp_chr /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr3_v3_s487395.sample /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipCallnomiss_CV_AB.ukb250k_all.tsv hasCHIP,hasDNMT3A,hasTET2,hasASXL1 Age,Age2,PC{1:10} Gender,Ethnicity,GenoBatch 400 /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/additive UKB250k 0.01 /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/250k_null/step1/NULL_MODEL_pred_hpc_modified.list 20 0.30 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg19_chrom_interval.tsv 10

# WB only 250k GWAS
# cat NULL_MODEL_pred.list | sed 's:/mnt/data/output/gs/ukbb_v2/projects/muddin/chip_gwas/UKB250k_WB_10Aug2022/step1/:/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/250k_WB_null/step1/:g' > NULL_MODEL_pred_hpc_modified.list
# qsub -wd /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg19_chrom_interval.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 10 -binding linear:10 -N WB.step2.ukb250k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/ukb450k.step02.regenie.sh /broad/ukbb/imputed_v3/ukb_imp_chr /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr3_v3_s487395.sample /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipCallnomiss_CV_AB.ukb250k_WB.tsv hasCHIP,hasDNMT3A,hasTET2 Age,Age2,PC{1:10} Gender,GenoBatch 400 /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_WB UKB250k_WB 0.01 /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/250k_WB_null/step1/NULL_MODEL_pred_hpc_modified.list 20 0.30 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg19_chrom_interval.tsv 10
####################################################

########################### Input files
BGEN_PREFIX=${1} # /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen
BGEN_SAMPLE_PREFIX=${2} # /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr3_v3_s487395.sample 
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

