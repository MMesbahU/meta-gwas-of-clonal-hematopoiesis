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


### TOPMed 
# FID	IID	NWD_ID	hasCHvaf02	hasCHvaf10	hasDNMT3A	hasTET2	hasASXL1	hasSF	hasDDR	AgeAtBloodDraw	Sqrd_AgeAtBloodDraw	Genetic_Sex	GenANC	STUDY	Sequencing_Center	TOPMed_Phase	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10	PC11	PC12	PC13	PC14	PC15	PC16	PC17	PC18	PC19	PC20
# MultiAncestry:
# mkdir -p /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir

# while read ANC; do qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N step1.${ANC}.TOPMed74k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step01.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/topmed_N74k.MultiANC.freeze.10b.chr1_22.plink_maf01hwe1e50geno10mind10 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.${ANC}.04_23_2024.v02.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:20},AgeAtBloodDraw,Sqrd_AgeAtBloodDraw" "Genetic_Sex,GenANC,STUDY,Sequencing_Center,TOPMed_Phase" /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/${ANC}/qcdGSAsnp.topmed.${ANC}.snplist /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/${ANC}/qcdGSAsnp.topmed.${ANC}.id 1000 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step1 4 30 5; done < <(echo -e "MultiANC")

## Sex Specific
# while read ANC; do qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N step1.${ANC}.TOPMed74k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step01.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/topmed_N74k.MultiANC.freeze.10b.chr1_22.plink_maf01hwe1e50geno10mind10 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.${ANC}.04_23_2024.v02.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:20},AgeAtBloodDraw,Sqrd_AgeAtBloodDraw" "GenANC,STUDY,Sequencing_Center,TOPMed_Phase" /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/${ANC}/qcdGSAsnp.topmed.${ANC}.snplist /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/${ANC}/qcdGSAsnp.topmed.${ANC}.id 1000 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step1 4 30 5; done < <(echo -e "Female\nMale")

## Ancestry Specific: 
## while read ANC; do qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N step1.${ANC}.TOPMed74k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step01.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/topmed_N74k.MultiANC.freeze.10b.chr1_22.plink_maf01hwe1e50geno10mind10 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.${ANC}.04_23_2024.v02.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:20},AgeAtBloodDraw,Sqrd_AgeAtBloodDraw" "Genetic_Sex,STUDY,Sequencing_Center,TOPMed_Phase" /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/${ANC}/qcdGSAsnp.topmed.${ANC}.snplist /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/${ANC}/qcdGSAsnp.topmed.${ANC}.id 1000 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step1 4 30 5; done < <(echo -e "AFR\nAMR\nEUR")

## while read ANC; do qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N step1.${ANC}.TOPMed74k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step01.regenie.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/topmed_N74k.MultiANC.freeze.10b.chr1_22.plink_maf01hwe1e50geno10mind10 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.${ANC}.04_23_2024.v02.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:20},AgeAtBloodDraw,Sqrd_AgeAtBloodDraw" "Genetic_Sex,STUDY,Sequencing_Center,TOPMed_Phase" /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/${ANC}/qcdGSAsnp.topmed.${ANC}.snplist /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/${ANC}/qcdGSAsnp.topmed.${ANC}.id 1000 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step1 4 30 1; done < <(echo -e "EAS\nSAS")

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

max_CatLevels=${11}

minCase=${12}

#########################################################################
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########

mkdir -p ${OUTPUT_PATH}

## Run Regenie
regenie \
	--step 1 \
	--minCaseCount ${minCase} \
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
	--maxCatLevels ${max_CatLevels} \
	--bt \
	--bsize ${BIN_SIZE} \
	--lowmem \
	--threads $(nproc --all) \
	--out ${OUTPUT_PATH}/NULL_MODEL


######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

