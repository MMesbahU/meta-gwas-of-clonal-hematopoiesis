#!/bin/bash

source /broad/software/scripts/useuse

use Bcftools

use Tabix

## mkdir -p /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/output/nov2024/{ukbb,topmed,mgbb,aou,biovu}/{afr,eur}
## mkdir -p /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/input/{metagwas,ukbb,topmed,mgbb,aou,biovu}
### UKBB
# ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/vcf/hg37.*.vcf.gz | awk '{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/ukbb.gwas_list.list
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N SumHer_UKBB_afr -l h_rt=5:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/ukbb.gwas_list.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/SumHer_LDAK.GWAS.sh /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.afr.tagging 0.05 0.05 /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/ukbb.gwas_list.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/input/ukbb /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/output/nov2024/ukbb/afr
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N SumHer_UKBB_eur -l h_rt=5:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/ukbb.gwas_list.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/SumHer_LDAK.GWAS.sh /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.gbr.tagging 0.05 0.05 /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/ukbb.gwas_list.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/input/ukbb /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/output/nov2024/ukbb/eur

### TOPMed
# ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/vcf/hg37.*.vcf.gz | awk '{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/topmed.gwas_list.list
## afr: qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N SumHer_TOPMed_afr -l h_rt=5:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/topmed.gwas_list.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/SumHer_LDAK.GWAS.sh /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.afr.tagging 0.05 0.05 /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/topmed.gwas_list.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/input/topmed /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/output/nov2024/topmed/afr

# eur: 
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N SumHer_TOPMed_eur -l h_rt=5:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/topmed.gwas_list.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/SumHer_LDAK.GWAS.sh /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.gbr.tagging 0.05 0.05 /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/topmed.gwas_list.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/input/topmed /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/output/nov2024/topmed/eur
####
## MGBB and AOU
## while read ANC; do ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/${ANC}/vcf/hg37.*.vcf.gz | awk '{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/${ANC}.gwas_list.list; done < <(echo -e "mgbb\naou")
# AFR:
# while read ANC; do qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N SumHer_${ANC}_afr -l h_rt=5:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/${ANC}.gwas_list.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/SumHer_LDAK.GWAS.sh /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.afr.tagging 0.05 0.05 /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/${ANC}.gwas_list.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/input/${ANC} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/output/nov2024/${ANC}/afr; done < <(echo -e "mgbb\naou")

# EUR:
# while read ANC; do qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N SumHer_${ANC}_eur -l h_rt=5:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/${ANC}.gwas_list.list | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/SumHer_LDAK.GWAS.sh /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.hapmap.gbr.tagging 0.05 0.05 /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/${ANC}.gwas_list.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/input/${ANC} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/output/nov2024/${ANC}/eur; done < <(echo -e "mgbb\naou")

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
####
######################
# Define file paths
taggFile=${1}

pop_prev=${2} 

gwas_prev=${3} 

LDAK=${4} 

list_hg37_summary_file=${5} 

ldak_inpuDir=${6}

ldak_outputDir=${7}

# Get GWAS file from task array list
hg37_summary_file=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file})

summary_LDAK_format=${ldak_inpuDir}/ldak.$(basename ${hg37_summary_file} ".vcf.gz").txt

out_file_prefix=${ldak_outputDir}/ldak_sumher.$(basename ${hg37_summary_file} ".vcf.gz")

# Check if summary file exists and is .gz
if [ -f "${summary_LDAK_format}" ]; then
    
    echo -e "File ${summary_LDAK_format} exists, proceeding without creating a new file.\n"

elif [ -f "${hg37_summary_file}.gz" ]; then
    
    echo -e "Found compressed file ${hg37_summary_file}.gz, uncompressing...\n"
    
    gzip -d ${hg37_summary_file}.gz

else
    
    echo -e "File ${summary_LDAK_format} does not exist. Creating the file in LDAK format.\n"
		    
    # Prepare LDAK-formatted GWAS summary
    echo -e "Predictor A1 A2 Z n" > ${summary_LDAK_format}
    bcftools view -i 'EAF>=0.005 & EAF<=0.995' ${hg37_summary_file} | \
    bcftools query -f '%CHROM:%POS %INFO/EffectAllele %INFO/OtherAllele %INFO/BETA %INFO/SE %INFO/N_Samples\n' | \
    awk '{snp = $1; a1 = $2; a2 = $3; z = $4 / $5; n = $6} ($1 ~ /^[0-9]+/ && (a1 == "A" || a1 == "C" || a1 == "G" || a1 == "T") && (a2 == "A" || a2 == "C" || a2 == "G" || a2 == "T")) {print snp, a1, a2, z, n}' | \
    awk '!seen[$1]++' >> ${summary_LDAK_format}

fi

# Run LDAK command
${LDAK} \
	--summary ${summary_LDAK_format} \
	--tagfile ${taggFile} \
	--sum-hers ${out_file_prefix} \
	--prevalence ${pop_prev} \
	--ascertainment ${gwas_prev} \
	--cutoff 0.01 \
	--check-sums NO
######################

gzip -f ${summary_LDAK_format}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########


