#!/bin/bash

## 2024 March 4:
# qsub -R y -wd /broad/hptmp/mesbah/dataset/ukb450 -N combine_plink_beds -l h_vmem=30G -l h_rt=5:00:00 -pe smp 6 -binding linear:6  /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/combine_ukb_BEDs.sh 1 6 /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/ukb_list_beds.txt /broad/hptmp/mesbah/dataset/ukb450
##

# qsub -wd /medpop/esp2/mesbah/CHIP_GWAS/ukb_v2 -N combine_plink_beds -l h_vmem=30G -l h_rt=5:00:00 -pe smp 6 -R y -binding linear:6  /medpop/esp2/mesbah/CHIP_GWAS/combine_ukb_BEDs.sh 1 6 /medpop/esp2/mesbah/CHIP_GWAS/ukb_list_beds.txt /medpop/esp2/mesbah/CHIP_GWAS/ukb_v2

# for chr in {2..22}; do echo "/broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed /broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam" >> ukb_list_beds.txt; done

source /broad/software/scripts/useuse

use .plink-1.90b

# input
chr=${1} # 1
numThreads=${2}
ukb_list_beds=${3} # ukb_list_beds.txt
outDir=${4}
plinkBEDfile=/broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed
plinkBIMfile=/broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim
plinkFAMfile=/medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam
outPlinkfile=${outDir}/ukb_v2_GT_chr1_22
# 
############################# Clock Time: Start
echo -e "Job started at: $(date)"

Job_START=$(date +%s)
############################

#plink2=/medpop/esp2/mesbah/tools/plink2
# ${plink2} 
plink \
	--bed ${plinkBEDfile} \
	--bim ${plinkBIMfile} \
	--fam ${plinkFAMfile} \
	--merge-list ${ukb_list_beds} \
	--threads ${numThreads} \
	--make-bed \
	--out ${outPlinkfile}

############################# Clock Time: END
echo "Job ended at: $(date)" 

Job_END=$(date +%s)


echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'

##############################


