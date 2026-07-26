#!/bin/bash -l
#$ -wd /medpop/esp2/lli/chip_gwas_rev/
#$ -N mgb_phewas
#$ -l h_vmem=20G
#$ -l h_rt=8:00:00
#$ -o /medpop/esp2/lli/chip_gwas_rev/logs/mgb_phewas_$TASK_ID.log
#$ -e /medpop/esp2/lli/chip_gwas_rev/logs/mgb_phewas_$TASK_ID.log
#$ -t 2-6

i=$(expr ${SGE_TASK_ID} - 1)

source /broad/software/scripts/useuse
use .r-4.0.0
use .gcc-13.2.0

# create looping parameters
phecode=($(cat /medpop/esp2/lli/chip_gwas_rev/Data/mgb_phewas/phecode_list.tsv))
phenCategory=${phecode[$i]}

Rscript /medpop/esp2/lli/chip_gwas_rev/Codes/phewas/mgbb/b_mgbb_phewas.R --phenCategory ${phenCategory}
