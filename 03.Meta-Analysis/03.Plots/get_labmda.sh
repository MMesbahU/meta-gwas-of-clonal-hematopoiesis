#!/bin/bash

source /broad/software/scripts/useuse

use R-4.1
# ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots/maf0.001*.tsv.gz | awk '{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/plots/list.of_summary.list 
# mkdir -p /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/plots/qqplots

## qsub -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/plots/list.of_summary.list | awk '{print $1}') -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/plots/tmpdir -N qqplots -l h_rt=2:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/get_labmda.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/get_labmda.R /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/plots/qqplots /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/plots/list.of_summary.list

##
rscript=${1}

outDir=${2} # path

gwas_summary_list=${3}
## 
gwas_summary=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${gwas_summary_list} )

qqplot_out_001=${outDir}/qqplot.eaf001.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${gwas_summary_list} ) ".tsv.gz").png

qqplot_out_002=${outDir}/qqplot.eaf002.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${gwas_summary_list} ) ".tsv.gz").png

qqplot_out_005=${outDir}/qqplot.eaf005.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${gwas_summary_list} ) ".tsv.gz").png

qqplot_out_01=${outDir}/qqplot.eaf01.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${gwas_summary_list} ) ".tsv.gz").png
## 
Rscript ${rscript} ${gwas_summary} ${qqplot_out_001} ${qqplot_out_002} ${qqplot_out_005} ${qqplot_out_01}

