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

############################ Job submission 
## --interaction VARNAME[BASE_LEVEL] (e.g. --interaction BMI[<25]).
## CH x Sex: 
# interection="CHxSex_refMale";ancestry="multi"; qsub -wd /broad/hptmp/mesbah/dataset/mgbb/gwas/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=20:00:00 -pe smp 8 -binding linear:8 -N step2.${interection}.${ancestry}.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/2.GxE/step02.GxE.regenie.sh /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,ancestry_pred,Batch_CHIP_call" 400 /broad/hptmp/mesbah/dataset/mgbb/gwas/MULTI/comb/${interection} ${interection}.${ancestry}.MGBB53k 0.01 /broad/hptmp/mesbah/dataset/mgbb/gwas/MULTI/comb/NULL_MODEL_pred.list 20 0.3 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv 8 "Sex[Male]"

## CH x Ancestry:
## interection="CHxAncestry_refEUR";ancestry="multi"; qsub -wd /broad/hptmp/mesbah/dataset/mgbb/gwas/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=20:00:00 -pe smp 8 -binding linear:8 -N step2.${interection}.${ancestry}.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/2.GxE/step02.GxE.regenie.sh /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,ancestry_pred,Batch_CHIP_call" 400 /broad/hptmp/mesbah/dataset/mgbb/gwas/MULTI/comb/${interection} ${interection}.${ancestry}.MGBB53k 0.01 /broad/hptmp/mesbah/dataset/mgbb/gwas/MULTI/comb/NULL_MODEL_pred.list 20 0.3 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv 8 "ancestry_pred[EUR]"
####################################### END 


####################################################
mkdir -p /broad/hptmp/mesbah/dataset/mgbb/gwas/MULTI/comb/${interection}
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
var_interaction=${16}
######

mkdir -p ${OUTPUT_PATH}

#####
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
# not using "--ref-first", coz of the bgen allele format in MGBB topmed imputed data
## Run Regenie
regenie \
	--step 2 \
	--bgen ${BGEN} \
	--sample ${BGEN_SAMPLE} \
	--bt \
	--htp ${COHORT} \
	--interaction ${var_interaction} \
	--no-condtl \
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
	--range ${CHR_RANGE} \
	--threads ${cpus} \
	--out ${OUT_PREFIX}

##################################

#####
####### Extract GxE summary
## zcat chr1.mac100.hfxch.ukb_v2_has_disease.regenie.gz | head -1 | cut -f1-21 |awk '{print $0"\tBETA\tSE\tMAC"}' | gzip -c > /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/chr1_22.HFxCH.tsv.gz

## for chr in {1..22}; do zcat chr${chr}.mac100.hfxch.ukb_v2_has_disease.regenie.gz | awk 'NR>1 {print $0}' | grep 'ADD-FIRTH-INT_SNPxCHIP' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:MAC=::g' | gzip -c >> /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/chr1_22.HFxCH.tsv.gz; done &

## 
# zcat /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/full_gxe/HFxCH_ukb_eur.chr1:1-24925062_has_disease.regenie.gz | head -1 | cut -f1-21 |awk '{print $0"\tBETA\tSE\tINFO\tMAC"}' | gzip -c > /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/chr1_22.ukb_all_eur.HFxCH.tsv.gz
# for files in $(ls -lhv /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/full_gxe/HFxCH_ukb_eur.chr*_has_disease.regenie.gz | awk '{print $NF}'); do zgrep 'ADD-FIRTH-INT_SNPxCHIP' ${files} | awk '$13>=0.01 && $13<=0.99' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | gzip -c >> /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/chr1_22.ukb_all_eur.HFxCH.tsv.gz; done &



#################################
######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

