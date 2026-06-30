#!/bin/bash

source /broad/software/scripts/useuse

# use Tabix
################
# ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.gz |  awk '{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/list_of_113_gwama_outfiles.list

###### qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/list_of_113_gwama_outfiles.list |awk '{print $1}') -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=10:00:00 -N gwama_meta_summary /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/prep_GWAMA_Summary.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/list_of_113_gwama_outfiles.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary "pval" 0.0001

###########################
GWAMA_outfile_list=${1} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_96_gwama_inputfiles.list

outDir=${2}

## Run GWAMA
outprefix=${3} # pal

plevel=${4} # 0.0001
## input files
gwama_summary=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${GWAMA_outfile_list} )

output_file=${outDir}/${outprefix}_${plevel}.$(basename ${gwama_summary} ".out.gz").tsv
############################

#################################################################
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
##################################
# zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.has${traits}.21Aug2021_ukbEUR.tsv.gz | awk '(NR==1){print "SNPID\tN\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $2~/^[0-9]+/){print $1"\t"($14+$18)"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12}' | awk '!seen[$1]++'| gzip -c > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.ukb200k.hg37_eur_${traits}.tsv.gz
##
zcat ${gwama_summary} | awk -v min_P=${plevel} 'NR==1 {print "SNPID\tREF\tALT\tAAF\tBETA\tSE\tBETA_95L\tBETA_95U\tZ\tP\tminusLOG10_P\tHet_Q\tHet_P\tI2\tN_Studies\tN_Samples\tDirection"; next}
     NR > 1 && $10 < min_P {print $1 "\t" $3 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17}' | \
     { head -n 1; tail -n +2 | sort -V -k1,1; } > ${output_file}

###########


######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

