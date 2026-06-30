#!/bin/bash

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

## Genotype QC in "prep_MGBB.sh"
#### MGBB
# while read ancestry; do qsub -wd /broad/hptmp/mesbah/dataset/mgbb/gwas/tmpdir -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N step1.${ancestry}.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/step01.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/qced.mgbb53k.autosomes.${ancestry}.hew1e15maf01 /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.${ancestry}.Feb2024.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,Batch_CHIP_call" /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/qced.mgbb53k.autosomes.${ancestry}.hew1e15maf01.snplist /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/qced.mgbb53k.autosomes.${ancestry}.hew1e15maf01.id 1000 /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/${ancestry}/step1 4; done < <(echo -e "AFR\nEAS\nSAS")


####

### 2024 March: UKB
# mkdir -p /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step1/{SAS,AMR,AFR,EAS}
# FID	IID	hasCH	hasCHvaf05	hasCHvaf10	hasDNMT3A	hasTET2	hasASXL1	hasDTA	hasSF	hasDDR	Ethnic_Background	sqrtAge_at_recruitment	GenoBatch	Age_at_recruitment	Genetic_Sex	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10	Batch	knn
## SNPs to keep: /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_pass_snps.*.chr1_22.maf01_hwe1e20_geno10pct.snplist
## Samples to keep (only positive IDs): /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/qc_pass_ukb_AFR_samples.v2.id 
#Nodes=4; while read ancestry; do qsub -wd /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step1/${ancestry} -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp ${Nodes} -binding linear:${Nodes} -N step1.${ancestry}.ukb450k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/step01.regenie.sh /broad/hptmp/mesbah/dataset/ukb450/ukb_v2_GT_chr1_22 /medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_n_250k_${ancestry}_N*.28cols.03_05_2024.tsv.gz "hasCH,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_at_recruitment,sqrtAge_at_recruitment" "Genetic_Sex,Batch,GenoBatch" /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_pass_snps.${ancestry}.chr1_22.maf01_hwe1e20_geno10pct.snplist /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/qc_pass_ukb_${ancestry}_samples.v2.id 1000 /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step1/${ancestry} ${Nodes}; done < <(echo -e "AFR\nAMR\nEAS\nSAS")
#####


#### Feb 2024
########### Updated MGBB GWAS including U2AF1 variants
## new pheno file: /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz
## header: FID	IID	Biobank_Subject_ID	Sex	Age_Genotyping	Vital_Status	Incd_cad	Prev_cad	Has_cad	age_enroll	batch	ancestry_pred	Ancestry_Self_cat	Batch_CHIP_call	hasCHvaf02	hasCHvaf10	hasDNMT3A	hasTET2	hasASXL1	hasSF	hasDDR	sqrAge_Genotyping	PC1	PC2	PC3PC4	PC5	PC6	PC7	PC8	PC9	PC10	PC11	PC12	PC13	PC14	PC15	PC16	PC1PC18	PC19	PC20
### Create dirs
# mkdir -p /broad/hptmp/mesbah/dataset/mgbb/gwas/{AFR,AMR,EAS,EUR,SAS,MULTI}
# mkdir -p /broad/hptmp/mesbah/dataset/mgbb/gwas/MULTI/Sex/{XX,XY}

## multi-ancestry:ancestry="multi"; qsub -wd /broad/hptmp/mesbah/dataset/mgbb/gwas/tmpdir -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 8 -binding linear:8 -N step1.${ancestry}.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/step01.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/qced.mgbb53k.autosomes.hew1e20maf01mac100 /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,ancestry_pred,Batch_CHIP_call" 1000 /broad/hptmp/mesbah/dataset/mgbb/gwas/MULTI/comb 8

## Ancestry Specific:
## while read lines; do zcat /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz| awk -v ancestry=${lines} 'NR==1{print $0}(NR>1 && $12==ancestry){print $0}' | gzip -c > /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.${lines}.Feb2024.tsv.gz; done < <(echo -e "AFR\nAMR\nEAS\nEUR\nSAS")

## while read ancestry; do qsub -wd /broad/hptmp/mesbah/dataset/mgbb/gwas/tmpdir -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 8 -binding linear:8 -N step1.${ancestry}.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/step01.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/qced.mgbb53k.autosomes.hew1e20maf01mac100 /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.${ancestry}.Feb2024.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,Batch_CHIP_call" 1000 /broad/hptmp/mesbah/dataset/mgbb/gwas/${ancestry} 8; done < <(echo -e "AFR\nAMR\nEAS\nEUR\nSAS")

## XX/XY Stratified MultiAncestry GWAS:
## while read lines; do zcat /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz| awk -v sex=${lines} 'NR==1{print $0}(NR>1 && $4==sex){print $0}' | gzip -c > /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.${lines}.Feb2024.tsv.gz; done < <(echo -e "Male\nFemale")
## XX: qsub -wd /broad/hptmp/mesbah/dataset/mgbb/gwas/tmpdir -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 8 -binding linear:8 -N step1.XX.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/step01.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/qced.mgbb53k.autosomes.hew1e20maf01mac100 /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Female.Feb2024.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "ancestry_pred,Batch_CHIP_call" 1000 /broad/hptmp/mesbah/dataset/mgbb/gwas/MULTI/Sex/XX 8
## XY: qsub -wd /broad/hptmp/mesbah/dataset/mgbb/gwas/tmpdir -R y -l h_vmem=5G -l h_rt=40:00:00 -pe smp 8 -binding linear:8 -N step1.XY.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/step01.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/qced.mgbb53k.autosomes.hew1e20maf01mac100 /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Male.Feb2024.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "ancestry_pred,Batch_CHIP_call" 1000 /broad/hptmp/mesbah/dataset/mgbb/gwas/MULTI/Sex/XY 8

####################################### END 
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
SNPs_2_Keep=${6}
Samples_2_Keep=${7}
BIN_SIZE=${8} # 1000
OUTPUT_PATH=${9} # /broad/hptmp/mesbah/gwas/mgb53k
cpus=${10}
#########################################################################
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
# mkdir -p ${OUTPUT_PATH}
## Run Regenie
regenie \
	--step 1 \
	--minCaseCount 1 \
	--loocv \
	--print-prs \
	--bed ${PLINK_Prefix} \
	--extract ${SNPs_2_Keep} \
	--keep ${Samples_2_Keep} \
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

