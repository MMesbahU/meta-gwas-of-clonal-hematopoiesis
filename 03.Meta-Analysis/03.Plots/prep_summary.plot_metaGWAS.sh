#!/bin/bash

source /broad/software/scripts/useuse

use Bcftools

use Tabix

# use Anaconda3
## Jun 2025
## modified: added filter
# qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/tmpdir -N prep_mcps_metasum -l h_rt=5:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/mcps_multi.meta6.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/prep_summary.plot_metaGWAS.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/mcps_multi.meta6.list 500000 2 0.01 1.00 /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots
#### end

##### multi-ancestry w/o mcps
# grep MultiANC.ukbb200k /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_113_gwama_output_files.list > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/multiancestry.n10.gwama_output_files.list
# qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/tmpdir -N prep_metasum_nst2_N400k_maf001 -l h_rt=1:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/multiancestry.n10.gwama_output_files.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/prep_summary.plot_metaGWAS.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/multiancestry.n10.gwama_output_files.list 400000 2 0.001 1.00 /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots

# qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/tmpdir -N prep_metasum_nst3_N400k_maf001 -l h_rt=1:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/multiancestry.n10.gwama_output_files.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/prep_summary.plot_metaGWAS.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/multiancestry.n10.gwama_output_files.list 400000 3 0.001 1.00 /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots
####################
# /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/Single_manhattan.py
# /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_113_gwama_output_files.list 

#
# qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/tmpdir -N prep_sum -l h_rt=5:00:00 -l h_vmem=30G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_113_gwama_output_files.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/prep_summary.plot_metaGWAS.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_113_gwama_output_files.list /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/plots

###################################
list_gwama_summary_file=${1}

n_sample_filter=${2} # >= N threshold

n_study_filter=${3} # >= number of studies has the snv

eaf_filter=${4} # >=0.01 i.e. >1%

p_filter=${5} # <0.99

outputDir=${6}

### get GWAS file from t array file list
gwama_summary_file=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_gwama_summary_file} )

summary_gwama_format=${outputDir}/maf${eaf_filter}_p${p_filter}_nsam${n_sample_filter}_nstd${n_study_filter}.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_gwama_summary_file} ) ".out.gz").tsv
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)

##### Prepare GWAMA summary for plot

## rs_number	reference_allele	other_allele	eaf	beta	se	beta_95L	beta_95U	z	p-value	_-log10_p-value	q_statistic	q_p-value	i2	n_studies	n_samples	effects
# chr1:10612:A:C	C	A	0.000051	1.574385	0.756204	0.092225	3.056545	2.081958	0.037366	1.427524	0.000000	1.000000	-nan	435506	+???

echo -e "SNPID\tEAF\tZ\tP\tN" > ${summary_gwama_format}
##
zcat ${gwama_summary_file} | awk -v n_sample=${n_sample_filter} -v n_study=${n_study_filter} -v eaf=${eaf_filter} -v p_filter=${p_filter} 'NR>1 && $4>=eaf && $4<=(1- eaf) && $10<p_filter && $15>=n_study && $16>=n_sample {print $1"\t"$4"\t"$9"\t"$10"\t"$16}' | sort -V -k1 >> ${summary_gwama_format} 

gzip -f ${summary_gwama_format}

#####
echo -e "N_SNV: $(zcat ${summary_gwama_format} | awk 'NR>1' | wc -l)\n"


## 
zcat ${summary_gwama_format} | awk -F '\t' 'NR==2 {min=max=$5} 
             NR>2 {
	        if ($5 < min) min = $5
		if ($5 > max) max = $5
	     } 
	     END {print "Min:", min, "Max:", max}'

#########

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########


