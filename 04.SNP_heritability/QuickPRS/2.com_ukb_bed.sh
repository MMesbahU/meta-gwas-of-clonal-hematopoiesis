#!/bin/bash

## Nov 2024:
# qsub -R y -wd /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/tmpdir -N comb_plink_beds -l h_vmem=20G -l h_rt=5:00:00 -pe smp 6 -binding linear:6  /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/QuickPRS/com_ukb_bed.sh 1 6 /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3/ukb_list_beds.txt /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3
##
# for chr in {2..22}; do echo "/medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3/ukb_imp_chr${chr}_v3.bed /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3/ukb_imp_chr${chr}_v3.bim /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3/ukb_imp_chr${chr}_v3.fam" >> /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3/ukb_list_beds.txt; done

source /broad/software/scripts/useuse

use .plink-1.90b

# input
chr=${1} # 1

numThreads=${2}

ukb_list_beds=${3} # ukb_list_beds.txt

outDir=${4}

plinkBEDfile=/medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3/ukb_imp_chr${chr}_v3.bed

plinkBIMfile=/medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3/ukb_imp_chr${chr}_v3.bim

plinkFAMfile=/medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3/ukb_imp_chr${chr}_v3.fam

outPlinkfile=${outDir}/ukb_imp_chr1_22_v3
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


