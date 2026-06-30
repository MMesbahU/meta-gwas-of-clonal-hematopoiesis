#!/bin/bash

## qsub -t 1-22 -R y -pe smp 2 -binding linear:2 -N add_hg38 -l h_rt=08:00:00 -l h_vmem=30G -wd /medpop/esp2/mesbah/Meta_GWAS/Aug2021/rsids/tmpdir /broad/hptmp/mesbah/dbSNP/run.add_rsids.GRCh38.sh /broad/hptmp/mesbah/dbSNP/add_rsids.GRCh38.R /medpop/esp2/mesbah/Meta_GWAS/Aug2021/eur_chip_maf1.29kTopmed_167kUKB_10kMGB_54kBioVU.tsv.gz 2 /medpop/esp2/mesbah/Meta_GWAS/Aug2021/rsids

source /broad/software/scripts/useuse

use R-4.0

myScript=${1}
GWAS=${2} # /medpop/esp2/mesbah/Meta_GWAS/Aug2021/eur_chip_maf1.29kTopmed_167kUKB_10kMGB_54kBioVU.tsv.gz
numThreads=${3}
outDir=${4}
lifoverlife=/broad/hptmp/mesbah/dbSNP/hg38/chr${SGE_TASK_ID}.dbSNP.hg38.tsv.gz
outFile=${outDir}/chr${SGE_TASK_ID}.$(basename ${GWAS})

##
Rscript ${myScript} ${GWAS} ${lifoverlife} ${numThreads} ${outFile}

