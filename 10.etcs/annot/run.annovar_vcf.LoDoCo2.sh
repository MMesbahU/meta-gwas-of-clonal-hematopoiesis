#!/bin/bash

### July 2025
# LoDoCo2:
# qsub -wd /medpop/esp2/mesbah/others/for_honigberg_lab/LoDoCo2 -R y -l h_vmem=20G -l h_rt=1:00:00 -pe smp 1 -binding linear:1 -N LoDoCo2_hg38 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/11.etcs/annot/run.annovar_vcf.LoDoCo2.sh /broad/hptmp/mesbah/dataset/annovar hg38 "refGene,ensGene,exac03,cytoBand,avsnp150,cosmic96_coding,cosmic70,gnomad211_exome,gnomad40_exome,gnomad40_genome,dbnsfp41a,clinvar_20200316" "g,f,f,f,f,f,f,f,f,f,f,f" /medpop/esp2/mesbah/others/for_honigberg_lab/LoDoCo2 /medpop/esp2/mesbah/others/for_honigberg_lab/LoDoCo2/LoDoCo2.CHIP_var.annot.vcf.gz
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

out_prefix=${outDir}/annovar_${ref_build}.$(basename ${in_vcf_file} ".vcf.gz")

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

# tabix -s 1 -b 33 -e 33 ${out_prefix}.${ref_build}_multianno.txt.gz

bgzip ${out_prefix}.${ref_build}_multianno.vcf

tabix -p vcf ${out_prefix}.${ref_build}_multianno.vcf.gz

########################################################################################################



