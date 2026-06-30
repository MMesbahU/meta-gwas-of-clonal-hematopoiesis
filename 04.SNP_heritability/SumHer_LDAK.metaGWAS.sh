#!/bin/bash

source /broad/software/scripts/useuse

use Bcftools

use Tabix

# Nov 5:
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N SumHer_EUR_notcons -l h_rt=5:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/gwas_list_113.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/SumHer_LDAK.sh /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.gbr.tagging 0.05 0.05 /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/gwas_list_113.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/input /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/output/nov2024/eur/notconstrained
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N SumHer_AFR -l h_rt=5:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/gwas_list_113.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/SumHer_LDAK.sh /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.afr.tagging 0.05 0.05 /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/gwas_list_113.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/input /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/output/nov2024/afr

## Nov 4
# while read ANC; do if [ "$ANC" == "AFR" ]; then TAGGING_FILE="/broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.afr.tagging"; elif [ "$ANC" == "EUR" ]; then TAGGING_FILE="/broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.gbr.tagging"; else echo "Unknown ANC: $ANC"; continue; fi; qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N SumHer${ANC} -l h_rt=10:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/gwas_list_${ANC}_v1.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/SumHer_LDAK.sh ${TAGGING_FILE} 0.05 0.05 /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/gwas_list_${ANC}_v1.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/input /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/output; done < <(echo -e "AFR\nEUR")

###################################
## Get the software
# wget https://dougspeed.com/wp-content/uploads/ldak5.2.linux_.zip
# unzip ldak5.2.linux_.zip
## HapMap tagging: bld.ldak.hapmap.gbr.tagging

#### http://dougspeed.com/pre-computed-tagging-files/
# Get precomputed Tagging SNPs
# LDAK-Thin Model: wget https://genetics.ghpc.au.dk/doug/ldak.thin.genotyped.gbr.tagging.gz
# BLD-LDAK-Lite-Alpha Model: wget https://genetics.ghpc.au.dk/doug/bld.ldak.lite.alpha.hapmap.gbr.tagging.gz
# BLD-LDAK Model: wget -P /broad/hptmp/mesbah/ukb_chip/ https://genetics.ghpc.au.dk/doug/bld.ldak.genotyped.gbr.tagging.gz

# ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/hg37/hg37.GWAMA.chr1_22.*.vcf.gz | awk 'NR<=30{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/gwas_list_30.list
# ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/hg37/hg37.GWAMA.chr1_22.*.vcf.gz | awk '{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/gwas_list_113.list

## Sumher: http://dougspeed.com/snp-heritability/
# /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.afr.tagging 0.05 0.05 /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/gwas_list_30.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/input /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/output

taggFile=${1}  # /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.afr.tagging

pop_prev=${2} # prevalence in general population

gwas_prev=${3} # prevalence in GWAS population

LDAK=${4} # /medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux

list_hg37_summary_file=${5} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/gwas_list_30.list

ldak_inpuDir=${6}

ldak_outputDir=${7}

### If the summary statistics come from analysing a binary phenotype, then you can use --prevalence <float> and --ascertainment <float> to specify the proportion of cases in the population and in the GWAS; LDAK will then also report estimates of variance explained on the liability scale.

### get GWAS file from t array file list
hg37_summary_file=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} )

summary_LDAK_format=${ldak_inpuDir}/ldak.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} ) ".vcf.gz").txt

out_file_prefix=${ldak_outputDir}/ldak_sumher.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} ) ".vcf.gz")

# pop_prev
# gwas_prev
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)

##### Prepare LDAK formated GWAS summary
# echo -e "Predictor A1 A2 Z n" > ${summary_LDAK_format}
#
# bcftools view -i 'EAF>=0.005 & EAF<=0.995 & N_Studies>=2' ${hg37_summary_file} | bcftools query -f '%CHROM:%POS %INFO/EffectAllele %INFO/OtherAllele %INFO/Z %INFO/N_Samples\n' | awk '{snp = $1; a1 = $2; a2 = $3; z = $4; n = $5} ($1 ~ /^[0-9]+/ && (a1 == "A" || a1 == "C" || a1 == "G" || a1 == "T") && (a2 == "A" || a2 == "C" || a2 == "G" || a2 == "T")) {print snp, a1, a2, z, n}' | awk '!seen[$1]++' >> ${summary_LDAK_format}
#########
${LDAK} \
	--summary ${summary_LDAK_format} \
	--tagfile ${taggFile} \
	--sum-hers ${out_file_prefix} \
	--prevalence ${pop_prev} \
	--ascertainment ${gwas_prev} \
	--cutoff 0.01 \
	--check-sums NO

######################

# gzip -f ${summary_LDAK_format}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########


