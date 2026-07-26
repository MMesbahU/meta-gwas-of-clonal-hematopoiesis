#!/bin/bash -l
#$ -wd /medpop/esp2/lli/chip_gwas_rev/
#$ -N ukb_phewas
#$ -l h_vmem=50G
#$ -l h_rt=24:00:00
#$ -o /medpop/esp2/lli/chip_gwas_rev/logs/ukb_phewas_$TASK_ID.log
#$ -e /medpop/esp2/lli/chip_gwas_rev/logs/ukb_phewas_$TASK_ID.log
#$ -t 1-6

i=$(expr ${SGE_TASK_ID} - 1)

source /broad/software/scripts/useuse
use .r-4.0.0
use .gcc-13.2.0

# create looping parameters
phecode=($(cat /medpop/esp2/lli/chip_gwas_rev/Data/ukb_phewas/phecode_list.tsv))
phenCategory=${phecode[$i]}

Rscript /medpop/esp2/lli/chip_gwas_rev/Codes/phewas/ukb/b_ukbb_phewas.R --phenCategory ${phenCategory}
