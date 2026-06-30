#!/bin/bash

source /broad/software/scripts/useuse

use Anaconda

source activate regenie_env

## for chr in {1..22}; do cp /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam /broad/hptmp/mesbah/ukb_chip/gxe/ukb_geno/ukb_gt_chr${chr}_v2.fam;done

# for chr in {1..22}; do cp /broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim /broad/hptmp/mesbah/ukb_chip/gxe/ukb_geno/ukb_gt_chr${chr}_v2.bim; done

# for chr in {1..22}; do cp /broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed /broad/hptmp/mesbah/ukb_chip/gxe/ukb_geno/ukb_gt_chr${chr}_v2.bed; done

############################ Job submission 

## CH x Sex: multi-ancestry:ancestry="multi"; qsub -wd /broad/hptmp/mesbah/dataset/mgbb/gwas/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv |awk '{print $1}') -R y -l h_vmem=5G -l h_rt=20:00:00 -pe smp 8 -binding linear:8 -N step2.${ancestry}.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/step02.regenie.sh /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,ancestry_pred,Batch_CHIP_call" 400 /broad/hptmp/mesbah/dataset/mgbb/gwas/MULTI/comb ${ancestry}.MGBB53k 0.01 /broad/hptmp/mesbah/dataset/mgbb/gwas/MULTI/comb/NULL_MODEL_pred.list 20 0.3 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv 8

## HF x CHIP:
## for chr in {1..22}; do qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/tmpdir -R y -l h_vmem=10G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N chr${chr}.ukb.hfxch /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/Step2.my_gxe.sh /broad/hptmp/mesbah/ukb_chip/gxe/ukb_geno/ukb_gt_chr${chr}_v2 /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/ukb450k_ukb200k_ukb250k.ch_hf.tsv has_disease "Age,Age2,CHIP,PC{1:10}" SEX,CHIP_Batch,Ethnicity,GenoBatch 1000 0.01 CHIP /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/chr${chr}.mac100.hfxch.ukb_v2 4; done
###########################

## Take input from the terminal
BED_FILE=${1}
PHENO_FILE=${2}
PHENO_COL=${3}
COVAR_COL_LIST=${4}
catCOVAR_COL_LIST=${5}
GENO_BLOCK_SIZE=${6}
P_THRESH=${7}
interaction_coVAR=${8}
outfile_prefix=${9}
N_threads=${10}
###### Run REGENIE Step2 w/o NULL model
regenie \
	--step 2 \
	--ignore-pred \
	--htp hfxch_ukb \
	--bed ${BED_FILE} \
	--phenoFile ${PHENO_FILE} \
	--phenoCol ${PHENO_COL} \
	--covarFile ${PHENO_FILE} \
	--covarColList ${COVAR_COL_LIST} \
	--catCovarList ${catCOVAR_COL_LIST} \
	--interaction ${interaction_coVAR} \
	--bsize ${GENO_BLOCK_SIZE} \
	--firth \
	--ref-first \
	--approx \
	--bt \
	--gz \
	--write-samples \
	--print-pheno \
	--pThresh ${P_THRESH} \
	--threads ${N_threads} \
	--minMAC 100 \
	--no-condtl \
	--out ${outfile_prefix}

####### Extract GxE summary
## zcat chr1.mac100.hfxch.ukb_v2_has_disease.regenie.gz | head -1 | cut -f1-21 |awk '{print $0"\tBETA\tSE\tMAC"}' | gzip -c > /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/chr1_22.HFxCH.tsv.gz

## for chr in {1..22}; do zcat chr${chr}.mac100.hfxch.ukb_v2_has_disease.regenie.gz | awk 'NR>1 {print $0}' | grep 'ADD-FIRTH-INT_SNPxCHIP' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:MAC=::g' | gzip -c >> /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/chr1_22.HFxCH.tsv.gz; done &

## 
zcat /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/full_gxe/HFxCH_ukb_eur.chr1:1-24925062_has_disease.regenie.gz | head -1 | cut -f1-21 |awk '{print $0"\tBETA\tSE\tINFO\tMAC"}' | gzip -c > /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/chr1_22.ukb_all_eur.HFxCH.tsv.gz
for files in $(ls -lhv /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/full_gxe/HFxCH_ukb_eur.chr*_has_disease.regenie.gz | awk '{print $NF}'); do zgrep 'ADD-FIRTH-INT_SNPxCHIP' ${files} | awk '$13>=0.01 && $13<=0.99' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | gzip -c >> /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/gxe.ch_hf/hf_ch/chr1_22.ukb_all_eur.HFxCH.tsv.gz; done &


