#!/bin/bash

### Top SNVs; P<0.0001
# ls  /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/unique_varID_hg38_p0001.vcf.gz | awk '{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/hg38.vcf.list
# qsub -wd /broad/hptmp/mesbah/dataset/GWAMA/out/annot/tmpdir -R y -t 1-1 -l h_vmem=20G -l h_rt=5:00:00 -pe smp 1 -binding linear:1 -N annovar_hg38_p0001 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/11.etcs/annot/run.annovar.vcf_list.sh /broad/hptmp/mesbah/dataset/annovar hg38 "refGene,ensGene,exac03,cytoBand,avsnp150,cosmic96_coding,cosmic70,gnomad211_exome,gnomad40_exome,gnomad40_genome,dbnsfp41a,clinvar_20200316" "g,f,f,f,f,f,f,f,f,f,f,f" /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/hg38.vcf.list

### Dec 2024
# mkdir -p /broad/hptmp/mesbah/dataset/GWAMA/out/annot/{tmpdir,out}
# ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/hg37/hg37.GWAMA.chr1_22.has*.vcf.gz | awk '{print $NF}' > /broad/hptmp/mesbah/dataset/GWAMA/out/annot/hg37.vcf.list
# qsub -wd /broad/hptmp/mesbah/dataset/GWAMA/out/annot/tmpdir -R y -t 1-$(wc -l /broad/hptmp/mesbah/dataset/GWAMA/out/annot/hg37.vcf.list | awk '{print $1}') -tc 30 -l h_vmem=20G -l h_rt=5:00:00 -pe smp 1 -binding linear:1 -N annovar_gwas_sum /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/11.etcs/annot/run.annovar.vcf_list.sh /broad/hptmp/mesbah/dataset/annovar hg19 "refGene,cosmic96_coding,cosmic70,gnomad211_exome,dbnsfp41a" "g,f,f,f,f" /broad/hptmp/mesbah/dataset/GWAMA/out/annot/out /broad/hptmp/mesbah/dataset/GWAMA/out/annot/hg37.vcf.list

####

## while read pheno; do bash /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/run.annovar.sh /broad/hptmp/mesbah/dataset/annovar /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/p5e6eaf001nstd2.metagwas.has${pheno}.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023.vcf.gz /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/annovar.p5e6eaf001nstd2.metagwas.has${pheno}.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023 1>>annovar.log 2>>annovar.err; done < <(echo -e "CH\nDNMT3A\nTET2") &

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

list_hg37_summary_file=${6}

### get GWAS file from t array file list
in_vcf_file=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} )

out_prefix=${outDir}/annot.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} ) ".vcf.gz")

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



