#!/bin/bash

# qsub -wd /broad/hptmp/mesbah/gwas/cojo/tmpdir -R y -l h_vmem=20G -l h_rt=20:00:00 -pe smp 1 -binding linear:1 -N sum_cojo.mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/06.COJO/01.prepMetaGWAS.ref_mgbb53k.sh mgbb53k_ref.summary_4_cojo /broad/hptmp/mesbah/gwas/cojo/input

source /broad/software/scripts/useuse

use Bcftools

# 
# MAF 0.1%
# n Studies>=2
cojo_input=${1}

cojoINdir=${2}

mkdir -p ${cojoINdir}

for files in $(ls -lhv /medpop/esp/mesbah/GWAMA_VCF/hg38/*hg38_dbSNP.vcf.gz | awk '{print $NF}')
do
	gwas_summary=${cojoINdir}/${cojo_input}.$(basename ${files} ".vcf.gz").tsv

#### COJO input from VCF file
# COJO input fields: 
	echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" > ${gwas_summary}
	
	bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2' ${files} | bcftools query -f '%INFO/varID_hg38\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/Pvalue\t%INFO/N_METAL\n' | awk '!seen[$1]++' >> ${gwas_summary}

done 


