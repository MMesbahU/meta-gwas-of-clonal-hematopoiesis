#!/bin/bash

## allUKB: while read pheno; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -N getinfo_allukb${pheno} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=20G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/01.PrepSummary/get_info.sh /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${pheno}.11Aug2021_ukb200k.vcf.gz /broad/hptmp/mesbah/ukb_chip/allUKB/tmp_info.chr1_22.has${pheno}.11Aug2021_ukb200k /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${pheno}.11Aug2021_ukb200k.INFO.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## qsub -R y -wd /broad/hptmp/mesbah/tmpdir -N getinfo_eurukb${pheno} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=30G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/01.PrepSummary/get_info.sh /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.hasCHIP.21Aug2021_ukbEUR.vcf.gz /broad/hptmp/mesbah/ukb_chip/eurUKB/tmp_info.chr1_22.hasCHIP.11Aug2021_ukb200k /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.hasCHIP.21Aug2021_ukbEUR.INFO.tsv

source /broad/software/scripts/useuse

use Tabix
use VCFtools

outVCF=${1}
out_info_file=${2}
final_output=${3}

##
tabix -h ${outVCF} $(echo chr{1..22}) | vcftools \
	--vcf - \
	--get-INFO AAF \
	--get-INFO REGENIE_BETA \
	--get-INFO REGENIE_SE \
	--get-INFO PVAL \
	--get-INFO N \
	--get-INFO old_hg19_id \
	--out ${out_info_file}

## Add CHROM:POS:REF:ALT 
head -1 ${out_info_file}.INFO | awk '{print "SNPID\t"$0}' > ${final_output}

awk 'NR>1{print $1":"$2":"$3":"$4"\t"$0}' ${out_info_file}.INFO | sort -k1 -V >> ${final_output}

rm ${out_info_file}.INFO

# compress 
bgzip -f ${final_output}

tabix -f -s 2 -b 3 -e 3 ${final_output}.gz
