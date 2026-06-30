#!/bin/bash

source /broad/software/scripts/useuse

use Bcftools

use Tabix

## ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/hg37/hg37.GWAMA.chr1_22.has*.{AFR,EUR,MultiANC}.AoU250k_TOPMed72k_MGBB53k.noukbb.vcf.gz | awk '{print $NF}' | sort -V > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/list_hg37_aou_top_mgbb.noukb.txt

# Nov 5:
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N prep_ldak_sum -l h_rt=1:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/list_hg37_aou_top_mgbb.noukb.txt | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/QuickPRS/prep_LDAK_summary.sh /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/chr1_22.afr_or_gbr.ldak.var.hapmap.cors.bim /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/list_hg37_aou_top_mgbb.noukb.txt "ALT" "REF" /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs

##
hapmapFile=${1}  # /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.afr.tagging

list_hg37_summary_file=${2} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/gwas_list_30.list

EffectAllele=${3}

OtherAllele=${4}

ldak_inpuDir=${5}

### If the summary statistics come from analysing a binary phenotype, then you can use --prevalence <float> and --ascertainment <float> to specify the proportion of cases in the population and in the GWAS; LDAK will then also report estimates of variance explained on the liability scale.

### get GWAS file from t array file list
hg37_summary_file=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} )

tmp_summary_LDAK_format=${ldak_inpuDir}/tmp_ldak.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} ) ".vcf.gz").txt

summary_LDAK_format=${ldak_inpuDir}/ldak.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} ) ".vcf.gz").txt

# out_file_prefix=${ldak_outputDir}/ldak_sumher.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} ) ".vcf.gz")

# pop_prev
# gwas_prev
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)

##### Prepare LDAK formated GWAS summary
echo -e "varid Predictor A1 A2 Z n" > ${tmp_summary_LDAK_format}
#
bcftools view -i 'EAF>=0.005 & EAF<=0.995 & N_Studies>=2' ${hg37_summary_file} | \
	bcftools query -f "%CHROM:%POS:%REF:%ALT %CHROM:%POS %${EffectAllele} %${OtherAllele} %INFO/Z %INFO/N_Samples\n" | awk '{varid = $1; snp = $2; a1 = $3; a2 = $4; z = $5; n = $6} ($1 ~ /^[0-9]+/ && (a1 == "A" || a1 == "C" || a1 == "G" || a1 == "T") && (a2 == "A" || a2 == "C" || a2 == "G" || a2 == "T")) {print varid, snp, a1, a2, z, n}' | awk '!seen[$2]++' >> ${tmp_summary_LDAK_format}

##
awk 'NR==FNR { vars[$1]; next } FNR==1 || ($1 in vars) ' <(awk '{print $2":"$5":"$6}' ${hapmapFile} ) ${tmp_summary_LDAK_format} | cut -d ' ' -f2-6 > ${summary_LDAK_format}

echo -e "row $(wc -l ${tmp_summary_LDAK_format} ) \n"

echo -e "hapmap: $(wc -l ${hapmapFile} ) \n"

echo -e "qc summary: $(wc -l ${summary_LDAK_format} ) \n"

rm ${tmp_summary_LDAK_format}
#########

######################

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########


