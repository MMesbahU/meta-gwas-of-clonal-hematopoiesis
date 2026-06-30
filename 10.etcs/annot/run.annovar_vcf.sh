#!/bin/bash

### mkdir -p /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/mash/2025_run/tmpdir 
# qsub -wd /broad/hptmp/mesbah/dataset/GWAMA/out/annot/tmpdir -R y -l h_vmem=20G -l h_rt=2:00:00 -pe smp 1 -binding linear:1 -N annovar_hg19 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/11.etcs/annot/run.annovar_vcf.sh /broad/hptmp/mesbah/dataset/annovar hg19 "refGene,ensGene,cosmic96_coding,cosmic70,gnomad211_exome,dbnsfp41a" "g,f,f,f,f,f" /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/mash/2025_run /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/mash/2025_run/varID.vcf

# qsub -wd /broad/hptmp/mesbah/dataset/GWAMA/out/annot/tmpdir -R y -l h_vmem=20G -l h_rt=2:00:00 -pe smp 1 -binding linear:1 -N annovar_hg19 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/11.etcs/annot/run.annovar_vcf.sh /broad/hptmp/mesbah/dataset/annovar hg19 "refGene,ensGene,exac03,cytoBand,avsnp150,cosmic96_coding,cosmic70,gnomad211_exome,gnomad40_exome,gnomad40_genome,dbnsfp41a,clinvar_20200316" "g,f,f,f,f,f,f,f,f,f,f,f" /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/hg38.vcf.list


###############################################
source /broad/software/scripts/useuse

use Tabix

use Bcftools

use R-4.0

use .gcc-7.3.0
#################
path_annovar=${1} # /broad/hptmp/mesbah/dataset/annovar 

ref_build=${2} # hg38

protocol=${3} # "refGene,cytoBand,avsnp150,cosmic96_coding,cosmic70,gnomad211_exome,dbnsfp41a,clinvar_20200316"

operations=${4} # "g,f,f,f,f,f,f,f"

outDir=${5}

in_vcf_file=${6}
#list_hg37_summary_file=${6}

### get GWAS file from t array file list
# in_vcf_file=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} )

# out_prefix=${outDir}/annot.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} ) ".vcf.gz")

out_prefix=${outDir}/annot_${ref_build}.$(basename ${in_vcf_file} ".vcf")

##
############################################### Annovar ################################################
## Run annovar
${path_annovar}/table_annovar.pl ${in_vcf_file} ${path_annovar}/humandb/ \
	-buildver ${ref_build} \
	--out ${out_prefix} \
	-remove \
	-protocol ${protocol} \
	-operation ${operations} \
	-nastring . \
	-vcfinput \
	-polish

# 
rm ${out_prefix}.avinput

bgzip ${out_prefix}.${ref_build}_multianno.txt

tabix -s 1 -b 33 -e 33 ${out_prefix}.${ref_build}_multianno.txt.gz

bgzip ${out_prefix}.${ref_build}_multianno.vcf

tabix -p vcf ${out_prefix}.${ref_build}_multianno.vcf.gz

########################################################################################################



