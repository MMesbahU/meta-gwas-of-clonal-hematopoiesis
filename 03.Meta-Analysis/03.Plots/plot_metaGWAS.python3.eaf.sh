#!/bin/bash

source /broad/software/scripts/useuse

# use Bcftools

# use Tabix

use Anaconda3
#-----------------
# Jan 23, 2025
# main_list.gwama.maf001_p09_nstd2.txt
# qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/tmpdir -N plot_manhattan_eaf001 -l h_rt=2:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 -t 1-$(wc -l /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots/main_list.gwama.maf001_p09_nstd2.txt | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/plot_metaGWAS.python3.eaf.sh /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots/main_list.gwama.maf001_p09_nstd2.txt /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/plots/maf001 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/Miami_Plot.n_qc.eaf.py 0.001
#-----------------

### December 3, 2024
# qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/tmpdir -N plot_manhattan_eaf01 -l h_rt=2:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 -t 1-$(wc -l /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots/list.gwama.maf001_p09_nstd2.txt | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/plot_metaGWAS.python3.eaf.sh /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots/list.gwama.maf001_p09_nstd2.txt /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/plots/maf01 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/Miami_Plot.n_qc.eaf.py 0.01


# ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots/maf001_p09_nstd2.GWAMA.chr1_22.has*.tsv.gz | awk '{print $NF}' | sort -V > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots/list.gwama.maf001_p09_nstd2.txt
###################################
list_gwama_summary_file=${1}

outputDir=${2}

py_manhattan=${3}

EAF=${4}

### get GWAS file from t array file list
input_file=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_gwama_summary_file} )

output_plot_mean_3sd=${outputDir}/sd3.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_gwama_summary_file} ) ".tsv.gz").png

output_plot_median_3mad=${outputDir}/mad3.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_gwama_summary_file} ) ".tsv.gz").png

output_plot_iqr=${outputDir}/iqr.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_gwama_summary_file} ) ".tsv.gz").png

pheno_data=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_gwama_summary_file} | sed -e 's:\.:\t:g' -e 's:has::g'| awk '{print $5":"$4}' )

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)

###########
# Run the Python script with provided arguments
python ${py_manhattan} ${input_file} ${output_plot_mean_3sd} ${output_plot_median_3mad} ${output_plot_iqr} ${pheno_data} ${EAF}
#########

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########


