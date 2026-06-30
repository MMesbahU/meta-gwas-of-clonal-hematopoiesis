#!/bin/bash

source /broad/software/scripts/useuse

use R-4.1
use Anaconda3
use Tabix

############################################
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/PheWAS/prs_phewas/tmpdir -t 1-$( echo -e "CirculatoryRespiratory\nHematologicNeoplasmInfectious\nDermDigestGU\nMentalNeuroSensorySymptoms\nInjuriesPoisoningMSK\nPregnancyCongenitalEndocrine" | wc -l |awk '{print $1}') -pe smp 4 -binding linear:4 -l h_vmem=30G -l h_rt=20:00:00 -N chip_prs_phewas /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.PheWAS/run.PheWAS.CHIP_PRS.UKB_all_phecodes.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.PheWAS/PheWAS.CHIP_PRS.UKB_all_phecodes.R /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/PheWAS/prs_phewas/ukb
#############################################
PheWAS_Parallel_script=${1}

outdir=${2}

Phecode_Group=$(echo -e "CirculatoryRespiratory\nHematologicNeoplasmInfectious\nDermDigestGU\nMentalNeuroSensorySymptoms\nInjuriesPoisoningMSK\nPregnancyCongenitalEndocrine"| awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' )
###################################################################

############### Task: Worked
## 1. Run PheWAS-Analysis
###########################

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
##################################
Rscript ${PheWAS_Parallel_script} ${outdir} ${Phecode_Group}

##################################

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

