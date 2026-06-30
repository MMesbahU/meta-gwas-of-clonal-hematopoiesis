#!/bin/bash

## qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -N meta4_summary.Aug2023 -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=20G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/04.gwama_Summary_Stats.sh /broad/hptmp/mesbah/dataset/ch_gwas/gwama_vcf_list.list 

##### Extract Variants from VCF files to Meta-analyse with AoU GWAS summary
## Need: SNPID   REF     ALT     AAF     BETA    SE      N       P

####
source /broad/software/scripts/useuse
use Bcftools
use Tabix

### Aug 2023
## ls -lhrt /broad/hptmp/mesbah/dataset/ch_gwas/metaGWAS.has*.topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023.vcf.gz | awk '{print $NF}' > /broad/hptmp/mesbah/dataset/ch_gwas/gwama_vcf_list.list
gwama_vcf_list=${1} # 



### Pvalue contains '0' in many instances; have to use 10^-log10_Pvalue
## convert log10(pvalue) to pvalue
while read lines
do 
	# SNPID   REF     ALT     AAF     BETA    SE      N       P
	out_file="$(dirname ${lines})"/"$(basename ${lines} ".vcf.gz").eaf001.tsv"

	echo -e "SNPID\tREF\tALT\tAAF\tBETA\tSE\tN\tP\tDirection_Effect" > ${out_file}
	# hg38 assembly: chromosome in chr1 format
	bcftools view  -i 'EAF>=0.001' ${lines} | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/log10_Pvalue\t%INFO/N_METAL\t%INFO/Effect_Direction\n' | awk '$1~/chr[0-9]+/{print $3"\t"$7"\t"$6"\t"$8"\t"$9"\t"$10"\t"$12"\t"10^-$11"\t"$13}' >> ${out_file}

	gzip -f ${out_file}

done <${gwama_vcf_list}


