#!/bin/bash

### October 2023
##### Extract Variants from VCF files
## EAF>=0.1% 
## N Studies>1
## 

## /home/jupyter/ch_gwas/aou250k/metaAnalyses/Oct2023/gwama_meta/metagwas*.vcf.gz
### Pvalue contains '0' in many instances; have to use 10^-log10_Pvalue
## convert log10(pvalue) to pvalue
##
## zgrep -v '^##' /home/jupyter/ch_gwas/aou250k/metaAnalyses/Oct2023/gwama_meta/metagwas.hasCH.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023.vcf.gz | head -2
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
# chr1    10416   chr1:10416:CCCTAA:C     CCCTAA  C       .       .       varID_hg38=chr1:10416:CCCTAA:C;EffectAllele=C;OtherAllele=CCCTAA;EAF=0.007152;BETA=0.070201;SE=0.110246;BETA_95L=-0.145881;BETA_95U=0.286283;Z=0.636767;Pvalue=0.524296;log10_Pvalue=0.280424;Q_stat=0;Q_Pvalue=1;I_sqr=-nan;N_Studies=1;N_Samples=166046;Effect_Direction=+?????;P_METAL=0.5243;N_METAL=166046;EAF_METAL=0.9928;BETA_METAL=-0.070201;SE_METAL=0.110246

### bash ~/myscripts/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/05.Summary_Stats.metaGWAS_N900k.sh 1>>log_file.txt 2>>err_file.txt &
##
while read trait; do echo -e "CHR\tPOS\tSNP\tREF\tATL\tEffectAllele\tOtherAllele\tEAF\tBETA\tSE\tlogP\tN\tDirection\tHet_P" > /home/jupyter/ch_gwas/aou250k/metaAnalyses/Oct2023/gwama_meta/has${trait}.multiAncestry.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023.eaf001_min2Studies.tsv && bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>1' /home/jupyter/ch_gwas/aou250k/metaAnalyses/Oct2023/gwama_meta/metagwas.has${trait}.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/log10_Pvalue\t%INFO/N_METAL\t%INFO/Effect_Direction\t%INFO/Q_Pvalue\n'  >> /home/jupyter/ch_gwas/aou250k/metaAnalyses/Oct2023/gwama_meta/has${trait}.multiAncestry.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023.eaf001_min2Studies.tsv; done < <(echo -e "CH\nDNMT3A\nTET2")





