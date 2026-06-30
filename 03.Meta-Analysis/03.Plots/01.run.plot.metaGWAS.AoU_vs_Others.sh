#!/bin/bash

############################ August 2023 ######################################
## AoU vs Others
##
## while read trait; do bash /home/jupyter/myscripts/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/01.run.plot.metaGWAS.AoU_vs_Others.sh /home/jupyter/myscripts/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/01.plot.metaGWAS.AoU_vs_Others.R gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/gwas/MetaAnalyses/meta5/maf001n2.metagwas.has${trait}.AoU250k_topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023.tsv.gz /home/jupyter/ch_gwas/aou250k/meta5 ${trait} maf001_N700k.Aug2023 700000 0.001 0.999; done < <(echo -e "CH\nDNMT3A\nTET2") 

###############################################################################
## Run : qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige/locusZoom -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=20:00:00 -N plot_metagwas.21Aug2021 /broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige/scripts/run.plot.metaGWAS.sh /broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige/scripts/plot.metaGWAS.R

##
plot_script=${1}

gwas_summary=${2}

Outdir=${3}

myTrait=${4}

pic_prefix=${5}

Sample_N=${6}

minMAF=${7}

maxMAF=${8}

####
Rscript ${plot_script} ${gwas_summary} ${Outdir} ${myTrait} ${pic_prefix} ${Sample_N} ${minMAF} ${maxMAF}

##
