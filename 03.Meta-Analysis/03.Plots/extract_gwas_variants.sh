#!/bin/bash

#######
### This script will grep SNV list from a GWAMA output file
######
source /broad/software/scripts/useuse
# MCPS 136k
# ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/*MCPS136k.out.gz | awk '{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/gwama.multi_afr_amr_vs_mcps136k.list
# mkdir -p /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/gwas_overlap/mcps 
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/gwas_overlap/tmpdir -N mcps.extract_snps -l h_rt=10:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/gwama.multi_afr_amr_vs_mcps136k.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/extract_gwas_variants.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/gwama.multi_afr_amr_vs_mcps136k.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/gwas_overlap/mcps /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/uniq.cojo.SNP_list.tsv

#------------
# mkdir -p /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/gwas_overlap/tmpdir
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/gwas_overlap/tmpdir -N extract_snps -l h_rt=20:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/n110_gwas.plus_stratified.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/extract_gwas_variants.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/n110_gwas.plus_stratified.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/gwas_overlap /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/uniq.cojo.SNP_list.tsv
#---
list_gwama_summary_file=${1} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/n110_gwas.plus_stratified.list
outDir=${2} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/gwas_overlap
query_snp_list=${3}
#--
gwas_file=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_gwama_summary_file} )

outFile=${outDir}/overlap_sr${SGE_TASK_ID}.$(basename ${gwas_file} ".gz").tsv

##
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)

###########
# Run
zgrep -f ${query_snp_list} ${gwas_file} >> ${outFile}
#########

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########
#---


