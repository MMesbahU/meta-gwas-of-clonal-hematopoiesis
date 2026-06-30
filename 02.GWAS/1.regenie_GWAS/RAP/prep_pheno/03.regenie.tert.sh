#!/bin/bash

source /broad/software/scripts/useuse
use Anaconda
## to install
## conda create -n regenie_env -c conda-forge -c bioconda regenie
## Update from v3.1.3g v3.2.5
# conda update -n regenie_env -c conda-forge -c bioconda regenie
## conda update -n base conda
source activate regenie_env

## regenie --step 2 --ignore-pred --bgen /broad/ukbb/imputed_v3/ukb_imp_chr5_v3.bgen --sample /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr3_v3_s487395.sample --bt --htp ukb450k --phenoFile /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_ukb250k_N421316.27cols.8Mar2023.tsv.gz --phenoColList hasCH,hasCHvaf05,hasCHvaf10,hasDNMT3A,hasTET2 --covarFile /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_ukb250k_N421316.27cols.8Mar2023.tsv.gz --covarColList Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} --catCovarList Genetic_Sex,Ethnic_Background,Batch,GenoBatch --bsize 400 --firth --ref-first --approx --gz --write-samples --print-pheno --pThresh 0.01 --minMAC 40 --minINFO 0.60 --range 5:1253262-1295184 --threads 4 --out /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/tert/ukb450k.chr5_1253262_1295184
####
# qsub -wd /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/tert/tmpdir -R y -l h_vmem=10G -l h_rt=10:00:00 -pe smp 4 -binding linear:4 -N tert.ukb450k /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/03.regenie.tert.sh /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_ukb250k_N421316.27cols.8Mar2023.tsv.gz hasCH,hasCHvaf05,hasCHvaf10,hasDNMT3A,hasTET2,hasDTA Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} Genetic_Sex,Ethnic_Background,Batch,GenoBatch /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/tert ukb450k 4

## qsub -wd /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/tert/tmpdir -R y -l h_vmem=10G -l h_rt=10:00:00 -pe smp 4 -binding linear:4 -N tert.ukb250k /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/03.regenie.tert.sh /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb250k_N243350.27cols.8Mar2023.tsv.gz hasCH,hasCHvaf05,hasCHvaf10,hasDNMT3A,hasTET2,hasDTA Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} Genetic_Sex,Ethnic_Background,GenoBatch /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/tert ukb250k 4

## qsub -wd /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/tert/tmpdir -R y -l h_vmem=10G -l h_rt=10:00:00 -pe smp 4 -binding linear:4 -N tert.ukb200k /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/03.regenie.tert.sh /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_N193342.27cols.8Mar2023.tsv.gz hasCH,hasCHvaf05,hasCHvaf10,hasDNMT3A,hasTET2,hasDTA Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} Genetic_Sex,Ethnic_Background,GenoBatch /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/tert ukb200k 4

#############
#BGEN_PREFIX=${1} # /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen
BGEN_SAMPLE_PREFIX=/medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr3_v3_s487395.sample 
PHENO_FILE=${1} # /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv
PHENO_COL_LIST=${2} # hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33
COVAR_COL_LIST=${3} # PC{1:10},Age_Genotyping,sqrAge_Genotyping
CAT_COVAR_COL_LIST=${4} # Sex,Ancestry_Self_cat,Batch_CHIP_call
GENO_BLOCK_SIZE=400
OUTPUT_PATH=${5} # /broad/hptmp/mesbah/gwas/mgb53k
COHORT=${6} # "ukb450k"
P_THRESH=0.01
MIN_MAC=40
INFO=0.60
cpus=${7}
#
chr=5
CHR_RANGE="5:1253262-1295184"
BGEN=/broad/ukbb/imputed_v3/ukb_imp_chr5_v3.bgen
BGEN_SAMPLE=${BGEN_SAMPLE_PREFIX}
OUT_PREFIX=${OUTPUT_PATH}/${COHORT}.chr${CHR_RANGE}

#############
## Run Regenie
regenie \
	--step 2 \
	--ignore-pred \
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
	--minMAC ${MIN_MAC} \
	--minINFO ${INFO} \
	--range ${CHR_RANGE} \
	--threads ${cpus} \
	--out ${OUT_PREFIX}



##### RAP

## regenie --step 1 --loocv --bed ukb_v2_GT_chr1_22 --extract qc_pass_snps.ukb_v2.list --keep qc_pass_ukb_all_samples.id --phenoFile CH_phenoCovar.ukb200k_N193342.27cols.8Mar2023.tsv.gz --phenoColList hasCH,hasDNMT3A,hasTET2,hasDTA --covarFile CH_phenoCovar.ukb200k_N193342.27cols.8Mar2023.tsv.gz --covarColList Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} --catCovarList Genetic_Sex,Ethnic_Background,GenoBatch --bt --bsize 1000 --lowmem --threads 16 --out NULL_MODEL.ukb200k

## regenie --step 1 --loocv --bed ukb_v2_GT_chr1_22 --extract qc_pass_snps.ukb_v2.list --keep qc_pass_ukb_all_samples.id --phenoFile CH_phenoCovar.ukb200k_N193342.27cols.8Mar2023.tsv.gz --phenoColList hasCH,hasDNMT3A,hasTET2,hasDTA --covarFile CH_phenoCovar.ukb200k_N193342.27cols.8Mar2023.tsv.gz --covarColList Age_at_recruitment,sqrtAge_at_recruitment,PC{1:10} --catCovarList Genetic_Sex,Ethnic_Background,GenoBatch --bt --bsize 1000 --lowmem --threads 16 --out NULL_MODEL.ukb200k

