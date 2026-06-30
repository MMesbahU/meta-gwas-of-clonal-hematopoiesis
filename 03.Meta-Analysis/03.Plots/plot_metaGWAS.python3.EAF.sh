#!/bin/bash

source /broad/software/scripts/useuse

# use Bcftools

# use Tabix

use Anaconda3
#-----------------
# Jun, 2025
# mcps summary list
# ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots/maf0.01_p1.00_nsam500000_nstd2.GWAMA.chr1_22.has*MCPS136k.tsv.gz | awk '{print $NF}' > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mcps.multi_gwas.eaf01.list
# ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots/maf0.001_p1.00_nsam500000_nstd2.GWAMA.chr1_22.has*MCPS136k.tsv.gz | awk '{print $NF}' > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mcps.multi_gwas.eaf001.list

#################
# qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/tmpdir -N plot_manhattan_eaf01 -l h_rt=2:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 -t 1-$(wc -l /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mcps.multi_gwas.eaf01.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/plot_metaGWAS.python3.EAF.sh /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mcps.multi_gwas.eaf01.list /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/plots/maf01 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/Miami_Plot.EAF.py 0.01
###################
##########################
## EAF >=0.001
# EAF>= 1%:  qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/tmpdir -N plot_manhattan_eaf001_01 -l h_rt=2:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 -t 1-$(wc -l /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mcps.multi_gwas.eaf001.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/plot_metaGWAS.python3.EAF.sh /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mcps.multi_gwas.eaf001.list /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/plots/maf001 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/Miami_Plot.EAF.py 0.01

# EAF>=0.5%: qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/tmpdir -N plot_manhattan_eaf001_005 -l h_rt=2:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 -t 1-$(wc -l /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mcps.multi_gwas.eaf001.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/plot_metaGWAS.python3.EAF.sh /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mcps.multi_gwas.eaf001.list /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/plots/maf001 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/Miami_Plot.EAF.py 0.005

# EAF>=0.2%: qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/tmpdir -N plot_manhattan_eaf001_002 -l h_rt=2:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 -t 1-$(wc -l /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mcps.multi_gwas.eaf001.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/plot_metaGWAS.python3.EAF.sh /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mcps.multi_gwas.eaf001.list /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/plots/maf001 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/Miami_Plot.EAF.py 0.002
####
### new list of gwas summary
# mkdir -p /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/plots/tmpdir
# myOutDir="/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/plots/manhattan"
# myGWASlist="/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/plots/list.of_summary.list"
# EAF>= 1%: 

# while read eaf; do qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/plots/tmpdir -N plot_manhattan_eaf${eaf} -l h_rt=2:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 -t 1-$(wc -l ${myGWASlist} | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/plot_metaGWAS.python3.EAF.sh ${myGWASlist} ${myOutDir} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/Miami_Plot.EAF.py ${eaf}; done < <(echo -e "0.01\n0.005\n0.002\n0.001")

####
##########################
#-----------------
###################################
list_gwama_summary_file=${1}

outputDir=${2}

py_manhattan=${3}

EAF=${4}

### get GWAS file from t array file list
input_file=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_gwama_summary_file} )

output_plot=${outputDir}/${EAF}.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_gwama_summary_file} ) ".tsv.gz").png
# 
pheno_data=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_gwama_summary_file} | sed -e 's:\.:\t:g' -e 's:has::g'| awk '{print $7":"$6}' )

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)

###########
# Run the Python script with provided arguments
python ${py_manhattan} ${input_file} ${output_plot} ${pheno_data} ${EAF}
#########

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########


